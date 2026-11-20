"""Sub-score calculators (V / M / S) for MAES v1.0.

All sub-scores are floats in ``[0, 1]`` where 1 is best.

V  validity
    ``V = V_PB ** (1 - alpha) * V_strain ** alpha`` with ``V_strain =
    exp(-dE_strain / kT)`` (RDKit UFF/MMFF relative to a relaxed conformer) and
    ``V_PB`` the PoseBusters pass rate.  PoseBusters is an *optional* backend:
    if the package is missing (offline sandbox) a built-in RDKit sanity checker
    is used and the backend name is recorded.

M  membrane consistency
    ``M = M_loc ** w1 * M_hyd ** w2 * M_ori ** w3`` using the projection of the
    ligand atoms onto the reference membrane plane (``membrane_map.MembranePlane``).
    ``M_hyd`` uses a *simplified atomic hydrophobicity proxy* (documented as
    such — it is not the published Wimley-White residue scale; pass an explicit
    ``hydrophobicity`` mapping to plug in your own scale).
    When no membrane plane is available M is neutral (1.0) and flagged.

S  state consistency
    ``S = exp(-d(phi_pred, phi_ref) / sigma)`` on conformational-state
    descriptors (ABC transporter: NBD interface distance, TM-helix angle,
    pocket volume, gating-residue distance).  A GMM classifier on reference
    descriptors is provided for the discrete IF/OF/intermediate label.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from .config import ScoreConfig
from .membrane_map import MembraneMap, MembranePlane

__all__ = [
    "ATOMIC_HYDROPHOBICITY",
    "CONFIDENCE_WEIGHTS",
    "DEFAULT_STATE_FIELDS",
    "MembraneResult",
    "StateResult",
    "ValidityResult",
    "build_subscores",
    "burial_factor",
    "confidence_weight",
    "fallback_validity",
    "gating_distance",
    "helix_axis_angle_deg",
    "insertion_energy",
    "ligand_principal_axis",
    "ligand_strain_energy",
    "membrane_consistency",
    "nbd_interface_distance",
    "posebusters_available",
    "posebusters_validity",
    "pocket_volume_proxy",
    "principal_axis_angle_deg",
    "state_consistency_score",
    "strain_term",
    "validity_score",
    "GMMStateClassifier",
]

# Boltzmann constant in kcal/mol/K (kT = 0.5925 kcal/mol at 298.15 K)
KB_KCAL = 0.00198720425864083

# Simplified atomic hydrophobicity proxy (arbitrary, documented as a proxy).
ATOMIC_HYDROPHOBICITY: Dict[str, float] = {
    "C": 1.0, "S": 0.6, "P": 0.4, "F": -0.2, "Cl": 0.5, "Br": 0.5,
    "I": 0.5, "N": -0.6, "O": -0.7, "H": 0.0, "Se": 0.4,
}

DEFAULT_STATE_FIELDS = ("nbd_distance", "tm_angle_deg", "pocket_volume", "gating_distance")

#: C-line ``membrane_confidence`` -> weight applied to the M sub-score
CONFIDENCE_WEIGHTS = {"high": 1.0, "medium": 0.5, "none": 0.0}


def confidence_weight(plane: Optional[MembranePlane]) -> float:
    """Weight for a membrane plane based on its C-line confidence grade."""
    if plane is None:
        return 0.0
    grade = str(getattr(plane, "confidence", "high")).strip().lower()
    return float(CONFIDENCE_WEIGHTS.get(grade, 1.0))


def _clip01(value: float) -> float:
    return float(min(1.0, max(0.0, float(value))))


# --------------------------------------------------------------------------- V
@dataclass
class ValidityResult:
    pass_rate: float
    n_checks: int
    n_passed: int
    backend: str
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pass_rate": float(self.pass_rate),
            "n_checks": int(self.n_checks),
            "n_passed": int(self.n_passed),
            "backend": self.backend,
            "details": dict(self.details),
        }


def posebusters_available() -> bool:
    """True when the optional PoseBusters backend can be imported."""
    try:  # pragma: no cover - depends on the environment
        import posebusters  # noqa: F401

        return True
    except Exception:
        return False


def _as_mol(mol_or_block: Any, remove_hs: bool = False):
    """Accept an RDKit Mol, a mol block or a SMILES string."""
    if mol_or_block is None:
        return None
    from rdkit import Chem

    if isinstance(mol_or_block, Chem.Mol):
        work = Chem.Mol(mol_or_block)
    else:
        text = str(mol_or_block).strip()
        work = Chem.MolFromMolBlock(text, removeHs=remove_hs, sanitize=True)
        if work is None:
            work = Chem.MolFromSmiles(text)
            if work is not None:
                from rdkit.Chem import AllChem

                work = Chem.AddHs(work)
                if AllChem.EmbedMolecule(work, randomSeed=42) != 0:
                    return None
    return work


def fallback_validity(mol_pred: Any) -> ValidityResult:
    """Cheap RDKit-only sanity checks used when PoseBusters is unavailable.

    Checks: (1) readable mol, (2) 3D conformer present, (3) no severe atom clash,
    (4) all covalent bonds have a plausible length, (5) connected molecule.
    This is intentionally a *lower bound* on the PoseBusters protocol.
    """
    mol = _as_mol(mol_pred)
    checks: Dict[str, bool] = {}
    if mol is None:
        return ValidityResult(0.0, 5, 0, "fallback", {"readable": False})
    checks["readable"] = True
    if mol.GetNumConformers() == 0:
        return ValidityResult(1 / 5, 5, 1, "fallback", checks)
    conf = mol.GetConformer()
    coords = np.asarray([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())], dtype=float)
    checks["has_3d_coords"] = bool(np.ptp(coords, axis=0).max() > 0.5)

    # severe clash: any pair of non-bonded heavy atoms closer than 0.4 A
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    clash = False
    if len(heavy) > 1:
        pos = coords[heavy]
        dist = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
        np.fill_diagonal(dist, np.inf)
        from rdkit import Chem

        bonded = set()
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if i in heavy and j in heavy:
                bonded.add((i, j))
                bonded.add((j, i))
        for i, ii in enumerate(heavy):
            for j, jj in enumerate(heavy):
                if i >= j:
                    continue
                if (ii, jj) in bonded:
                    continue
                if dist[i, j] < 0.4:
                    clash = True
    checks["no_atom_clash"] = not clash

    # bond lengths plausible
    bond_ok = True
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        d = float(np.linalg.norm(coords[i] - coords[j]))
        if not 0.5 <= d <= 2.6:
            bond_ok = False
            break
    checks["plausible_bonds"] = bond_ok
    checks["connected"] = len(_fragment_sizes(mol)) == 1
    n_passed = int(sum(1 for ok in checks.values() if ok))
    return ValidityResult(n_passed / 5.0, 5, n_passed, "fallback", checks)


def _fragment_sizes(mol) -> list:
    from rdkit import Chem

    return [len(frag) for frag in Chem.GetMolFrags(mol)]


def _select_posebusters_config(
    pdb_path: Optional[str],
    mol_true: Any,
    config: Optional[str],
) -> str:
    """Choose a PoseBusters config compatible with the available inputs."""
    if config is not None:
        return str(config)
    if pdb_path is not None and mol_true is not None:
        return "redock"
    if pdb_path is not None:
        return "dock"
    if mol_true is not None:
        return "redock"
    return "mol"


def posebusters_validity(
    mol_pred: Any,
    mol_true: Any = None,
    pdb_path: Optional[str] = None,
    config: Optional[str] = None,
) -> ValidityResult:
    """Run PoseBusters (0.6.x ``bust`` API, with 0.3.x ``check`` fallback).

    PoseBusters 0.6.5 (installed 2026-11-16) exposes
    ``PoseBusters(config).bust(mol_pred, mol_true, mol_cond)``; older releases
    used ``check(..., pdb=...)``.  Both are supported.  When the call raises,
    the built-in :func:`fallback_validity` is returned with the error recorded
    in ``details`` so the paper can see which backend produced V.
    """
    mol = _as_mol(mol_pred)
    if mol is None:
        return ValidityResult(0.0, 1, 0, "none", {"readable": False})
    cfg = _select_posebusters_config(pdb_path, mol_true, config)
    try:
        from posebusters import PoseBusters

        buster = PoseBusters(config=cfg)
        kwargs: Dict[str, Any] = {"mol_pred": mol}
        if mol_true is not None:
            kwargs["mol_true"] = _as_mol(mol_true)
        if pdb_path is not None:
            kwargs["mol_cond"] = str(pdb_path)
        if hasattr(buster, "bust"):  # posebusters >= 0.4
            table = buster.bust(**kwargs)
        else:  # posebusters 0.3.x
            legacy: Dict[str, Any] = {"mol_pred": mol}
            if mol_true is not None:
                legacy["mol_true"] = _as_mol(mol_true)
            if pdb_path is not None:
                legacy["pdb"] = str(pdb_path)
            table = buster.check(**legacy)  # type: ignore[attr-defined]
        bool_cols = [c for c in table.columns if pd.api.types.is_bool_dtype(table[c])]
        if not bool_cols:
            raise RuntimeError("PoseBusters returned no boolean check columns")
        row = table[bool_cols].iloc[0]
        n_passed = int(row.sum())
        n_checks = int(len(bool_cols))
        details = {str(k): bool(v) for k, v in row.items()}
        try:
            import importlib.metadata as _md

            details["posebusters_version"] = _md.version("posebusters")
        except Exception:  # pragma: no cover
            details["posebusters_version"] = "unknown"
        details["config"] = cfg
        return ValidityResult(n_passed / n_checks, n_checks, n_passed, "posebusters", details)
    except Exception as exc:  # missing package OR runtime failure
        result = fallback_validity(mol)
        result.details["posebusters_error"] = f"{type(exc).__name__}: {exc}"
        result.backend = "fallback (posebusters error)"
        return result


def ligand_strain_energy(
    mol_pred: Any,
    n_confs: int = 10,
    forcefield: str = "UFF",
    seed: int = 42,
    max_iters: int = 500,
) -> Optional[float]:
    """RDKit strain energy of a pose relative to its relaxed conformer (kcal/mol).

    Returns ``None`` when the ligand cannot be parameterised.  Values below zero
    are clipped to 0 (the pose was already at/below the reference minimum).
    """
    mol = _as_mol(mol_pred)
    if mol is None or mol.GetNumConformers() == 0:
        return None
    try:  # pragma: no cover - requires rdkit force-field machinery
        from rdkit import Chem
        from rdkit.Chem import AllChem

        work = Chem.Mol(mol)
        if forcefield.upper() == "MMFF":
            if AllChem.MMFFGetMoleculeProperties(work) is None:
                return None
            pose_ff = AllChem.MMFFGetMoleculeForceField(work, AllChem.MMFFGetMoleculeProperties(work), confId=0)
        else:
            pose_ff = AllChem.UFFGetMoleculeForceField(work, confId=0)
        if pose_ff is None:
            return None
        e_pose = float(pose_ff.CalcEnergy())

        ref = Chem.Mol(mol)
        ref.RemoveAllConformers()
        if all(atom.GetAtomicNum() != 1 for atom in ref.GetAtoms()):
            ref = Chem.AddHs(ref)
        params = AllChem.ETKDGv3()
        params.randomSeed = int(seed)
        conf_ids = list(AllChem.EmbedMultipleConfs(ref, numConfs=int(n_confs), params=params))
        if not conf_ids:
            return None
        if forcefield.upper() == "MMFF":
            energies = [e for e, _ in AllChem.MMFFOptimizeMoleculeConfs(ref, maxIters=int(max_iters))]
        else:
            energies = [e for e, _ in AllChem.UFFOptimizeMoleculeConfs(ref, maxIters=int(max_iters))]
        if not energies:
            return None
        return float(max(0.0, e_pose - min(energies)))
    except Exception:
        return None


def strain_term(
    strain_kcal: Optional[float],
    temperature_k: float = 298.15,
    mode: str = "exp",
    sigma_kcal: float = 3.0,
) -> float:
    """Map a strain energy to ``[0, 1]`` (1 = strain-free)."""
    if strain_kcal is None or not np.isfinite(float(strain_kcal)):
        return 1.0
    strain = max(0.0, float(strain_kcal))
    if mode == "exp":
        kt = KB_KCAL * float(temperature_k)
        if kt <= 0:
            raise ValueError("temperature_k must be positive")
        return _clip01(np.exp(-strain / kt))
    if mode == "squared":
        if sigma_kcal <= 0:
            raise ValueError("sigma_kcal must be positive")
        return _clip01(np.exp(-0.5 * (strain / sigma_kcal) ** 2))
    raise ValueError(f"unknown strain mode: {mode!r}")


def validity_score(
    pb_pass_rate: float,
    strain_kcal: Optional[float] = None,
    alpha: float = 0.2,
    temperature_k: float = 298.15,
    strain_mode: str = "exp",
    sigma_kcal: float = 3.0,
) -> float:
    """Combine PoseBusters pass rate and ligand strain into V in ``[0, 1]``."""
    if not 0.0 <= float(alpha) <= 1.0:
        raise ValueError("alpha must be in [0, 1]")
    v_pb = _clip01(pb_pass_rate)
    v_strain = strain_term(strain_kcal, temperature_k=temperature_k, mode=strain_mode, sigma_kcal=sigma_kcal)
    v = (v_pb ** (1.0 - alpha)) * (v_strain ** alpha) if alpha > 0 else v_pb
    return _clip01(v)


# --------------------------------------------------------------------------- M
@dataclass
class MembraneResult:
    M: float
    M_loc: float
    M_hyd: float
    M_ori: float
    neutral: bool = False
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "M": float(self.M),
            "M_loc": float(self.M_loc),
            "M_hyd": float(self.M_hyd),
            "M_ori": float(self.M_ori),
            "neutral": bool(self.neutral),
            "details": dict(self.details),
        }


def _neutral_membrane(reason: str) -> MembraneResult:
    return MembraneResult(1.0, 1.0, 1.0, 1.0, neutral=True, details={"reason": reason})


def burial_factor(z: Any, half_thickness: float, interface_width: float = 8.0) -> np.ndarray:
    """Fraction of membrane burial in ``[0, 1]`` for signed distances ``z``."""
    if half_thickness <= 0:
        raise ValueError("half_thickness must be positive")
    z_abs = np.abs(np.asarray(z, dtype=float))
    factor = np.zeros_like(z_abs)
    core = z_abs < half_thickness
    factor[core] = 1.0
    if interface_width > 0:
        iface = (z_abs >= half_thickness) & (z_abs < half_thickness + interface_width)
        factor[iface] = 1.0 - (z_abs[iface] - half_thickness) / interface_width
    return np.clip(factor, 0.0, 1.0)


def insertion_energy(
    z: Any,
    elements: Sequence[str],
    half_thickness: float,
    interface_width: float = 8.0,
    hydrophobicity: Optional[Mapping[str, float]] = None,
) -> float:
    """Simplified membrane insertion free energy (proxy, not the WW residue scale)."""
    z_arr = np.asarray(z, dtype=float).reshape(-1)
    elements = [str(e) for e in elements]
    if z_arr.size != len(elements):
        raise ValueError("z and elements must have the same length")
    scale = dict(hydrophobicity) if hydrophobicity is not None else dict(ATOMIC_HYDROPHOBICITY)
    hydro = np.array([float(scale.get(e.strip().capitalize(), 0.0)) for e in elements])
    return float(np.sum(hydro * burial_factor(z_arr, half_thickness, interface_width)))


def ligand_principal_axis(coords: Any) -> np.ndarray:
    """First principal component of the ligand heavy-atom coordinates."""
    xyz = np.asarray(coords, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or xyz.shape[0] < 2:
        raise ValueError("coords must be an (n>=2, 3) array")
    centred = xyz - xyz.mean(axis=0)
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    axis = np.asarray(vt[0], dtype=float)
    norm = float(np.linalg.norm(axis))
    if norm == 0:
        raise ValueError("degenerate ligand coordinates")
    return axis / norm


def principal_axis_angle_deg(coords: Any, normal: Any) -> float:
    """Angle (deg, in ``[0, 90]``) between the ligand principal axis and the membrane normal."""
    axis = ligand_principal_axis(coords)
    n = np.asarray(normal, dtype=float).reshape(3)
    n = n / np.linalg.norm(n)
    cos = float(abs(np.dot(axis, n)))
    return float(np.degrees(np.arccos(np.clip(cos, 0.0, 1.0))))


def membrane_consistency(
    cfg: ScoreConfig,
    z_lig: Optional[Any] = None,
    z_ref: Optional[Any] = None,
    plane: Optional[MembranePlane] = None,
    angle_lig_deg: Optional[float] = None,
    angle_ref_deg: Optional[float] = None,
    dG_lig: Optional[float] = None,
    dG_ref: Optional[float] = None,
    interface_width: Optional[float] = None,
    centroid_lig: Optional[Any] = None,
    centroid_ref: Optional[Any] = None,
    confidence_weight: float = 1.0,
) -> MembraneResult:
    """M sub-score: location + hydrophobic match + orientation consistency.

    ``centroid_lig`` / ``centroid_ref`` (optional 3-vectors) add an in-plane
    lateral-displacement penalty to ``M_loc``; without them ``M_loc`` only sees
    the membrane-normal projection (design doc v0.1 behaviour).

    ``confidence_weight`` shrinks M toward the neutral value 1.0 according to
    the C-line ``membrane_confidence`` grade: high = 1.0 (full weight),
    medium = 0.5, none/unavailable = 0.0 (neutral + flagged).  This keeps a
    low-confidence geometric plane from dominating the score.
    """
    if plane is None or z_lig is None or z_ref is None:
        return _neutral_membrane("no membrane plane / no reference projection")
    z_lig_arr = np.asarray(z_lig, dtype=float).reshape(-1)
    z_ref_arr = np.asarray(z_ref, dtype=float).reshape(-1)
    if z_lig_arr.size == 0 or z_ref_arr.size == 0:
        return _neutral_membrane("empty projection")
    width = float(interface_width) if interface_width is not None else float(cfg.interface_width_angstrom)

    # --- location: soft |dz| tolerance with a hard compartment-consistency penalty
    zbar_lig = float(np.mean(z_lig_arr))
    zbar_ref = float(np.mean(z_ref_arr))
    soft = float(np.exp(-abs(zbar_lig - zbar_ref) / float(cfg.membrane_tau_angstrom)))
    comp_lig = str(plane.compartment(zbar_lig, interface_width=width))
    comp_ref = str(plane.compartment(zbar_ref, interface_width=width))
    m_loc = soft * (1.0 if comp_lig == comp_ref else 0.5)

    lateral = None
    if centroid_lig is not None and centroid_ref is not None:
        c_lig = np.asarray(centroid_lig, dtype=float).reshape(-1)
        c_ref = np.asarray(centroid_ref, dtype=float).reshape(-1)
        delta = c_lig - c_ref
        normal = plane.normal_array
        perp = delta - float(np.dot(delta, normal)) * normal
        lateral = float(np.linalg.norm(perp))
        m_loc *= float(np.exp(-lateral / float(cfg.membrane_tau_angstrom)))

    # --- hydrophobic match
    if dG_lig is not None and dG_ref is not None:
        m_hyd = float(np.exp(-abs(float(dG_lig) - float(dG_ref)) / float(cfg.membrane_lambda_kcal)))
    else:
        m_hyd = 1.0

    # --- orientation
    if angle_lig_deg is not None and angle_ref_deg is not None:
        m_ori = float(np.exp(-abs(float(angle_lig_deg) - float(angle_ref_deg)) / float(cfg.membrane_angle_tau_deg)))
    else:
        m_ori = 1.0

    weights = cfg.membrane_weights
    m = (
        _clip01(m_loc) ** float(weights.get("loc", 1 / 3))
        * _clip01(m_hyd) ** float(weights.get("hyd", 1 / 3))
        * _clip01(m_ori) ** float(weights.get("ori", 1 / 3))
    )
    w = float(max(0.0, min(1.0, confidence_weight)))
    m = 1.0 - w + w * _clip01(m)  # shrink toward neutral for low-confidence planes
    return MembraneResult(
        M=_clip01(m),
        M_loc=_clip01(m_loc),
        M_hyd=_clip01(m_hyd),
        M_ori=_clip01(m_ori),
        neutral=bool(w <= 0.0),
        details={
            "z_lig_mean": zbar_lig,
            "z_ref_mean": zbar_ref,
            "compartment_lig": comp_lig,
            "compartment_ref": comp_ref,
            "lateral_deviation_angstrom": lateral,
            "confidence_weight": w,
            "confidence": str(getattr(plane, "confidence", "high")),
            "method": str(getattr(plane, "method", "")),
        },
    )


# --------------------------------------------------------------------------- S
@dataclass
class StateResult:
    S: float
    distance: float
    neutral: bool = False
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "S": float(self.S),
            "distance": float(self.distance),
            "neutral": bool(self.neutral),
            "details": dict(self.details),
        }


def nbd_interface_distance(a_coords: Any, b_coords: Any) -> float:
    """Minimum interatomic distance between two NBD half-transporters (Angstrom)."""
    a = np.asarray(a_coords, dtype=float).reshape(-1, 3)
    b = np.asarray(b_coords, dtype=float).reshape(-1, 3)
    if a.size == 0 or b.size == 0:
        raise ValueError("empty coordinate sets")
    dist = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)
    return float(dist.min())


def helix_axis_angle_deg(axis_a: Any, axis_b: Any) -> float:
    """Angle between two (unit or unnormalised) helix axes, in ``[0, 90]``."""
    a = np.asarray(axis_a, dtype=float).reshape(3)
    b = np.asarray(axis_b, dtype=float).reshape(3)
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    return float(np.degrees(np.arccos(np.clip(abs(float(np.dot(a, b))), 0.0, 1.0))))


def pocket_volume_proxy(pocket_coords: Any) -> float:
    """Convex-hull volume of pocket atoms (Angstrom^3) as a cheap pocket-size proxy."""
    xyz = np.asarray(pocket_coords, dtype=float).reshape(-1, 3)
    if xyz.shape[0] < 4:
        raise ValueError("need at least 4 pocket atoms")
    try:
        from scipy.spatial import ConvexHull

        return float(ConvexHull(xyz).volume)
    except Exception:
        # degraded proxy: axis-aligned bounding-box volume of the pocket atoms
        span = np.ptp(xyz, axis=0)
        return float(np.prod(span))


def gating_distance(res_a_coords: Any, res_b_coords: Any) -> float:
    """Minimum distance between two gating-residue atom sets (Angstrom)."""
    return nbd_interface_distance(res_a_coords, res_b_coords)


def _phi_vector(phi: Any, fields: Sequence[str] = DEFAULT_STATE_FIELDS) -> np.ndarray:
    if isinstance(phi, Mapping):
        return np.asarray([float(phi.get(f, np.nan)) for f in fields], dtype=float)
    if isinstance(phi, str):
        return _phi_vector(json.loads(phi), fields=fields)
    arr = np.asarray(phi, dtype=float).reshape(-1)
    return arr


def state_consistency_score(
    phi_pred: Any,
    phi_ref: Any,
    sigma: float = 1.0,
    scale: Optional[np.ndarray] = None,
) -> StateResult:
    """Continuous state consistency ``S = exp(-d(phi_pred, phi_ref) / sigma)``."""
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    if phi_pred is None or phi_ref is None:
        return StateResult(1.0, 0.0, neutral=True, details={"reason": "no state descriptors"})
    a = _phi_vector(phi_pred)
    b = _phi_vector(phi_ref)
    if a.shape != b.shape:
        raise ValueError("phi_pred and phi_ref must have the same length")
    diff = np.abs(a - b)
    if scale is not None:
        scale = np.asarray(scale, dtype=float).reshape(-1)
        if scale.shape != diff.shape:
            raise ValueError("scale must have the same length as phi")
        diff = diff / np.where(scale > 0, scale, 1.0)
    distance = float(np.nanmean(diff))
    return StateResult(float(np.exp(-distance / float(sigma))), distance)


class GMMStateClassifier:
    """Gaussian mixture model over reference state descriptors (skeleton)."""

    def __init__(self, n_components: int = 2, seed: int = 42):
        if n_components < 1:
            raise ValueError("n_components must be >= 1")
        self.n_components = int(n_components)
        self.seed = int(seed)
        self._gmm = None
        self.state_names: Dict[int, str] = {}

    # ------------------------------------------------------------------ fitting
    def fit(self, phi_ref: Any, state_names: Optional[Sequence[str]] = None) -> "GMMStateClassifier":
        x = np.asarray(phi_ref, dtype=float)
        if x.ndim != 2:
            raise ValueError("phi_ref must be an (n_samples, n_features) array")
        if x.shape[0] < self.n_components:
            raise ValueError(
                f"need >= n_components samples to fit GMM ({x.shape[0]} < {self.n_components})"
            )
        from sklearn.mixture import GaussianMixture

        self._gmm = GaussianMixture(
            n_components=self.n_components,
            covariance_type="full",
            random_state=self.seed,
            n_init=5,
            reg_covar=1e-4,
        ).fit(x)
        # order clusters by the first descriptor so that labels are deterministic
        order = np.argsort(self._gmm.means_[:, 0])
        self.state_names = {}
        for rank, cluster in enumerate(order):
            if state_names is not None and rank < len(state_names):
                self.state_names[int(cluster)] = str(state_names[rank])
            else:
                self.state_names[int(cluster)] = f"state_{rank}"
        return self

    def predict(self, phi: Any) -> np.ndarray:
        if self._gmm is None:
            raise RuntimeError("GMMStateClassifier is not fitted")
        x = np.asarray(phi, dtype=float)
        if x.ndim == 1:
            x = x.reshape(1, -1)
        clusters = self._gmm.predict(x)
        names = np.array([self.state_names.get(int(c), f"cluster_{c}") for c in clusters], dtype=object)
        return names

    def predict_proba(self, phi: Any) -> np.ndarray:
        if self._gmm is None:
            raise RuntimeError("GMMStateClassifier is not fitted")
        x = np.asarray(phi, dtype=float)
        if x.ndim == 1:
            x = x.reshape(1, -1)
        return self._gmm.predict_proba(x)

    def to_dict(self) -> Dict[str, Any]:
        if self._gmm is None:
            raise RuntimeError("GMMStateClassifier is not fitted")
        return {
            "n_components": self.n_components,
            "seed": self.seed,
            "weights": self._gmm.weights_.tolist(),
            "means": self._gmm.means_.tolist(),
            "covariances": self._gmm.covariances_.tolist(),
            "state_names": {str(k): v for k, v in self.state_names.items()},
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "GMMStateClassifier":
        model = cls(n_components=int(data["n_components"]), seed=int(data.get("seed", 42)))
        from sklearn.mixture import GaussianMixture

        gmm = GaussianMixture(n_components=model.n_components, covariance_type="full")
        gmm.weights_ = np.asarray(data["weights"], dtype=float)
        gmm.means_ = np.asarray(data["means"], dtype=float)
        gmm.covariances_ = np.asarray(data["covariances"], dtype=float)
        gmm.precisions_cholesky_ = np.linalg.cholesky(np.linalg.inv(gmm.covariances_))
        gmm.n_features_in_ = gmm.means_.shape[1]
        gmm.converged_ = True
        model._gmm = gmm
        model.state_names = {int(k): v for k, v in data.get("state_names", {}).items()}
        return model


# ------------------------------------------------------------ table-level glue
def build_subscores(
    df: pd.DataFrame,
    cfg: ScoreConfig,
    membrane_map: Optional[MembraneMap] = None,
) -> pd.DataFrame:
    """Compute the ``V`` / ``M`` / ``S`` columns from a Track-B pose table.

    Expected primitive columns (all optional except ``pos`` identity)::

        pose_id, system_id
        pb_pass_rate  or posebusters_pass     (V)
        strain_kcal                           (V, optional)
        z_lig_mean, z_ref_mean                (M)
        angle_lig_deg, angle_ref_deg          (M, optional)
        dG_lig, dG_ref                        (M, optional)
        phi_lig, phi_ref                      (S, optional; sequences or JSON)
        V, M, S                               (pre-computed; used if primitives absent)

    Missing membrane/state information yields neutral (1.0) terms flagged in
    ``m_neutral`` / ``s_neutral`` so the paper can report coverage.
    """
    out = df.copy()
    n = len(out)
    V = np.ones(n, dtype=float)
    M = np.ones(n, dtype=float)
    S = np.ones(n, dtype=float)
    m_neutral = np.zeros(n, dtype=bool)
    m_confidence = np.array(["none"] * n, dtype=object)
    m_frame_mismatch = np.zeros(n, dtype=bool)
    s_neutral = np.zeros(n, dtype=bool)
    v_backend = np.array(["precomputed"] * n, dtype=object)

    for i, (_, row) in enumerate(out.iterrows()):
        # --- V
        if "V" in out.columns and not pd.isna(row.get("V")):
            V[i] = _clip01(row["V"])
            v_backend[i] = "precomputed"
        else:
            if "pb_pass_rate" in out.columns and not pd.isna(row.get("pb_pass_rate")):
                pb = float(row["pb_pass_rate"])
            elif "posebusters_pass" in out.columns and not pd.isna(row.get("posebusters_pass")):
                pb = 1.0 if bool(row["posebusters_pass"]) else 0.0
            else:
                raise KeyError("need 'V', 'pb_pass_rate' or 'posebusters_pass' for the V sub-score")
            strain = row.get("strain_kcal") if "strain_kcal" in out.columns else None
            strain = None if strain is None or pd.isna(strain) else float(strain)
            V[i] = validity_score(
                pb,
                strain,
                alpha=cfg.alpha,
                temperature_k=cfg.strain_temperature_k,
                strain_mode=cfg.strain_mode,
                sigma_kcal=cfg.strain_sigma_kcal,
            )
            v_backend[i] = "posebusters+UFF" if strain is not None else "posebusters"

        # --- M
        sid = str(row.get("system_id", ""))
        if "membrane_key" in out.columns and not pd.isna(row.get("membrane_key")):
            sid = str(row["membrane_key"])
        plane = membrane_map.get(sid) if membrane_map is not None else None
        # frame consistency: a pose expressed in RCSB coordinates must not be
        # scored against an OPM-frame plane (and vice versa); mismatch is
        # reported as neutral M instead of silently comparing frames
        pose_frame = None
        if "coord_frame" in out.columns and not pd.isna(row.get("coord_frame")):
            raw = str(row["coord_frame"]).strip().lower()
            pose_frame = "opm" if raw in {"opm", "opm_oriented", "opm-oriented"} else (
                "rcsb" if raw in {"rcsb", "rcsb_original", "rcsb-original"} else None
            )
        frame_mismatch = bool(
            plane is not None and pose_frame is not None
            and str(getattr(plane, "coord_frame", "opm")) != pose_frame
        )
        z_lig = row.get("z_lig_mean") if "z_lig_mean" in out.columns else None
        z_ref = row.get("z_ref_mean") if "z_ref_mean" in out.columns else None
        if frame_mismatch:
            M[i] = 1.0
            m_neutral[i] = True
            m_confidence[i] = "none"
            m_frame_mismatch[i] = True
        elif plane is not None and z_lig is not None and z_ref is not None and not (
            pd.isna(z_lig) or pd.isna(z_ref)
        ):
            weight = confidence_weight(plane)
            res = membrane_consistency(
                cfg,
                z_lig=[float(z_lig)],
                z_ref=[float(z_ref)],
                plane=plane,
                angle_lig_deg=row.get("angle_lig_deg") if "angle_lig_deg" in out.columns else None,
                angle_ref_deg=row.get("angle_ref_deg") if "angle_ref_deg" in out.columns else None,
                dG_lig=row.get("dG_lig") if "dG_lig" in out.columns else None,
                dG_ref=row.get("dG_ref") if "dG_ref" in out.columns else None,
                confidence_weight=weight,
            )
            M[i] = res.M
            m_neutral[i] = bool(res.neutral)
            m_confidence[i] = str(getattr(plane, "confidence", "high"))
        elif "M" in out.columns and not pd.isna(row.get("M")):
            M[i] = _clip01(row["M"])
        else:
            M[i] = 1.0
            m_neutral[i] = True

        # --- S
        phi_lig = row.get("phi_lig") if "phi_lig" in out.columns else None
        phi_ref = row.get("phi_ref") if "phi_ref" in out.columns else None
        if phi_lig is not None and phi_ref is not None and not (
            (isinstance(phi_lig, float) and pd.isna(phi_lig)) or (isinstance(phi_ref, float) and pd.isna(phi_ref))
        ):
            res_s = state_consistency_score(phi_lig, phi_ref, sigma=cfg.state_sigma)
            S[i] = res_s.S
            s_neutral[i] = res_s.neutral
        elif "S" in out.columns and not pd.isna(row.get("S")):
            S[i] = _clip01(row["S"])
        else:
            S[i] = 1.0
            s_neutral[i] = True

    out["V"] = np.clip(V, 0.0, 1.0)
    out["M"] = np.clip(M, 0.0, 1.0)
    out["S"] = np.clip(S, 0.0, 1.0)
    out["m_neutral"] = m_neutral
    out["m_confidence"] = m_confidence
    out["m_frame_mismatch"] = m_frame_mismatch
    out["s_neutral"] = s_neutral
    out["v_backend"] = v_backend
    return out
