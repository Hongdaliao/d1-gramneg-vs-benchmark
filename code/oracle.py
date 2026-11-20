"""Oracle harness: validate the MAES pipeline on the 68 reference (co-crystal) poses.

The B-line predicted poses are not available yet, so this module uses the
**reference ligands themselves as "oracle predictions"** and adds controlled
rigid-body perturbations with known RMSD.  This makes it possible to

1. sanity-check the V/M/S *numerics* on real crystal geometries and real
   membrane planes (OPM-oriented frames, ``membrane_planes.json`` W2);
2. exercise the full calibration -> evaluation chain (AUROC / Brier / ECE /
   rho / divergence cases / ablation) before B-line data arrive;
3. check that the M term behaves correctly: in-pocket perturbations keep M high,
   lateral "wrong pocket" and out-of-membrane translations drive M down.

It is explicitly **not** a scientific result: V is constant under rigid-body
perturbation (PoseBusters checks are internal to the ligand), and S is neutral
because a rigid ligand transform does not change the protein state.  S will be
exercised once B-line provides predicted receptor/state arms.

Run::

    python -c "from maes.oracle import run_oracle; run_oracle()"
"""

from __future__ import annotations

import gzip
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from .calibration import MAESModel, label_fingerprint
from .config import MAESConfig, default_config
from .evaluation import (
    ablation_study,
    auroc,
    brier_score,
    ece,
    reliability_curve,
    spearman_rho,
)
from .membrane_map import MembranePlane, TableMembraneMap, load_c_line_planes
from .subscores import (
    confidence_weight,
    insertion_energy,
    ligand_principal_axis,
    ligand_strain_energy,
    membrane_consistency,
    posebusters_validity,
    principal_axis_angle_deg,
    validity_score,
)

__all__ = ["OracleReference", "load_reference_poses", "perturb_pose", "run_oracle"]

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_STRUCT_DIR = PROJECT_ROOT / "data" / "trackB" / "cache" / "struct"
DEFAULT_MANIFEST = PROJECT_ROOT / "data" / "trackB" / "pose_manifest_curated.csv"
DEFAULT_PLANES = PROJECT_ROOT / "data" / "trackB" / "membrane_planes.json"
DEFAULT_STATES = PROJECT_ROOT / "data" / "trackB" / "state_descriptors.csv"


# --------------------------------------------------------------------- helpers
def _read_lines(path: Path) -> List[str]:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", errors="ignore") as fh:
            return fh.read().splitlines()
    return path.read_text(encoding="utf-8", errors="ignore").splitlines()


def _coord_path(pdb_id: str, struct_dir: Path) -> Optional[Path]:
    for suffix in (".opm.pdb", ".pdb", ".pdb.gz"):
        candidate = struct_dir / f"{pdb_id}{suffix}"
        if candidate.exists():
            return candidate
    return None


def extract_ligand(lines: Sequence[str], ligand_id: str, chain: Optional[str]) -> Optional[Dict[str, Any]]:
    """Extract ligand HETATM records and build an RDKit mol from that PDB block."""
    from rdkit import Chem

    lig = str(ligand_id)
    records = []
    for line in lines:
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip()
        if resname != lig:
            continue
        chain_id = line[21:22].strip()
        altloc = line[16:17]
        if chain is not None and str(chain).strip() and chain_id != str(chain).strip():
            continue
        if altloc not in (" ", "", "A"):
            continue
        records.append(line)
    if len(records) < 5:
        return None
    block = "\n".join(records + ["END"])
    mol = Chem.MolFromPDBBlock(block, sanitize=True, removeHs=False, proximityBonding=True)
    if mol is None:
        return None
    coords = np.asarray([list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())], dtype=float)
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    return {"mol": mol, "coords": coords, "elements": elements, "n_atoms": len(records), "chain": records[0][21:22].strip()}


def _with_coords(mol, coords: np.ndarray):
    from rdkit import Chem

    new = Chem.Mol(mol)
    new.RemoveAllConformers()
    conf = Chem.Conformer(new.GetNumAtoms())
    for i in range(new.GetNumAtoms()):
        conf.SetAtomPosition(i, [float(x) for x in coords[i]])
    new.AddConformer(conf, assignId=True)
    return new


@dataclass
class OracleReference:
    pose_id: str
    pdb_id: str
    target: str
    ligand_class: str
    chain: str
    coords: np.ndarray
    elements: List[str]
    mol: Any
    plane: Optional[MembranePlane]
    plane_source: str
    z_ref: np.ndarray
    angle_ref_deg: float
    dG_ref: float
    centroid_ref: np.ndarray
    buried_frac: float
    state_phi: Optional[List[float]] = None
    pb_pass_rate: float = 1.0
    pb_backend: str = "unknown"
    pb_n_checks: int = 0
    strain_kcal: Optional[float] = None
    V: float = 1.0


def load_reference_poses(
    manifest_path: str | Path = DEFAULT_MANIFEST,
    planes_path: str | Path = DEFAULT_PLANES,
    struct_dir: str | Path = DEFAULT_STRUCT_DIR,
    states_path: Optional[str | Path] = DEFAULT_STATES,
    subset: str = "primary_main",
    cfg: Optional[MAESConfig] = None,
    strain_confs: int = 5,
) -> tuple:
    """Load usable reference poses (coords + mol + membrane plane + state descriptors)."""
    from .api import ANALYSIS_SUBSETS, select_analysis_set

    cfg = (cfg or default_config()).validate()
    manifest = pd.read_csv(Path(manifest_path))
    if "is_primary" not in manifest.columns:
        from .api import attach_pose_ids

        manifest = attach_pose_ids(manifest)
    subset_frame = select_analysis_set(manifest, subset)
    planes = load_c_line_planes(planes_path) if Path(planes_path).exists() else TableMembraneMap({})
    state = pd.read_csv(Path(states_path)) if states_path and Path(states_path).exists() else pd.DataFrame()

    refs: List[OracleReference] = []
    skipped: List[Dict[str, str]] = []
    for _, row in subset_frame.iterrows():
        pdb_id = str(row["pdb_id"])
        path = _coord_path(pdb_id, Path(struct_dir))
        if path is None:
            skipped.append({"pose_id": str(row["pose_id"]), "reason": "no coord file"})
            continue
        lig = extract_ligand(_read_lines(path), str(row["ligand_id"]), row.get("chain"))
        if lig is None:
            skipped.append({"pose_id": str(row["pose_id"]), "reason": "ligand not parsed"})
            continue
        coords = lig["coords"]
        plane = planes.get(pdb_id)
        # Missing plane is NOT a reason to drop the pose: production behaviour is
        # to keep it with a neutral M (m_neutral=True), so the oracle harness does
        # the same and reports coverage.
        if plane is not None:
            z_ref = plane.project(coords)
            angle_ref = principal_axis_angle_deg(coords, plane.normal_array)
            dG_ref = insertion_energy(
                z_ref, lig["elements"], plane.half_thickness,
                interface_width=cfg.scores.interface_width_angstrom,
            )
            plane_source = str(plane.source)
        else:
            z_ref = np.zeros(len(coords), dtype=float)
            angle_ref = float("nan")
            dG_ref = float("nan")
            plane_source = "unavailable"
        buried = 0.5
        phi = None
        if len(state):
            match = state[
                (state["pdb_id"].astype(str) == pdb_id)
                & (state["ligand_id"].astype(str) == str(row["ligand_id"]))
            ]
            if len(match):
                buried = float(pd.to_numeric(match["ligand_buried_frac"], errors="coerce").fillna(0.5).iloc[0])
                phi_values = [
                    pd.to_numeric(match.get(column), errors="coerce").iloc[0]
                    for column in ("tm_angle_mean_deg", "pocket_void_A3", "chain_interface_min_A")
                ]
                if all(np.isfinite(v) for v in phi_values):
                    phi = [float(v) for v in phi_values]
        pb = posebusters_validity(lig["mol"], config="mol")
        strain = ligand_strain_energy(lig["mol"], n_confs=int(strain_confs))
        v = validity_score(
            pb.pass_rate, strain,
            alpha=cfg.scores.alpha,
            temperature_k=cfg.scores.strain_temperature_k,
            strain_mode=cfg.scores.strain_mode,
            sigma_kcal=cfg.scores.strain_sigma_kcal,
        )
        refs.append(
            OracleReference(
                pose_id=str(row["pose_id"]), pdb_id=pdb_id, target=str(row["target"]),
                ligand_class=str(row["ligand_class"]), chain=str(row.get("chain", "")),
                coords=coords, elements=lig["elements"], mol=lig["mol"], plane=plane,
                plane_source=plane_source, z_ref=z_ref, angle_ref_deg=angle_ref,
                dG_ref=dG_ref, centroid_ref=coords.mean(axis=0), buried_frac=buried,
                state_phi=phi, pb_pass_rate=float(pb.pass_rate), pb_backend=pb.backend,
                pb_n_checks=int(pb.n_checks), strain_kcal=(None if strain is None else float(strain)), V=float(v),
            )
        )
    return refs, skipped


# --------------------------------------------------------------- perturbations
PERTURBATION_KINDS = ("in_pocket", "near_miss", "wrong_pocket", "out_of_membrane")


def perturb_pose(
    ref: OracleReference,
    kind: str,
    rng: np.random.Generator,
    interface_width: float = 8.0,
) -> Dict[str, Any]:
    """Rigid-body perturbation with known RMSD (same receptor frame)."""
    from scipy.spatial.transform import Rotation

    coords = ref.coords.copy()
    centroid = coords.mean(axis=0)
    normal = ref.plane.normal_array if ref.plane is not None else np.array([0.0, 0.0, 1.0])
    half = float(ref.plane.half_thickness) if ref.plane is not None else 15.0

    if kind == "in_pocket":
        shift = rng.normal(size=3)
        shift = shift / (np.linalg.norm(shift) + 1e-9) * rng.uniform(0.0, 1.2)
        angle_deg = rng.uniform(0.0, 8.0)
        axis = rng.normal(size=3)
    elif kind == "near_miss":
        shift = rng.normal(size=3)
        shift = shift / (np.linalg.norm(shift) + 1e-9) * rng.uniform(1.2, 3.5)
        angle_deg = rng.uniform(8.0, 25.0)
        axis = rng.normal(size=3)
    elif kind == "wrong_pocket":
        lateral = rng.normal(size=3)
        lateral = lateral - np.dot(lateral, normal) * normal  # in-plane
        lateral = lateral / (np.linalg.norm(lateral) + 1e-9) * rng.uniform(4.0, 10.0)
        shift = lateral
        angle_deg = rng.uniform(10.0, 40.0)
        axis = rng.normal(size=3)
    elif kind == "out_of_membrane":
        sign = rng.choice([-1.0, 1.0])
        shift = sign * normal * (half + interface_width + rng.uniform(1.0, 8.0))
        angle_deg = rng.uniform(10.0, 50.0)
        axis = rng.normal(size=3)
    else:
        raise ValueError(f"unknown perturbation kind: {kind!r}")

    rotation = Rotation.from_rotvec(axis / (np.linalg.norm(axis) + 1e-9) * np.deg2rad(angle_deg))
    perturbed = rotation.apply(coords - centroid) + centroid + shift
    rmsd = float(np.sqrt(np.mean(np.sum((perturbed - ref.coords) ** 2, axis=1))))
    return {"coords": perturbed, "rmsd": rmsd, "kind": kind, "shift": shift, "angle_deg": float(angle_deg)}


# ------------------------------------------------------------------ scoring
def _score_perturbed(ref: OracleReference, perturbed: Dict[str, Any], cfg: MAESConfig) -> Dict[str, Any]:
    coords = perturbed["coords"]
    z_lig = ref.plane.project(coords) if ref.plane is not None else np.zeros(len(coords))
    angle = principal_axis_angle_deg(coords, ref.plane.normal_array) if ref.plane is not None else ref.angle_ref_deg
    dG = insertion_energy(
        z_lig, ref.elements,
        ref.plane.half_thickness if ref.plane is not None else 15.0,
        interface_width=cfg.scores.interface_width_angstrom,
    )
    centroid = coords.mean(axis=0)
    if ref.plane is not None:
        weight = confidence_weight(ref.plane)
        membrane = membrane_consistency(
            cfg.scores,
            z_lig=z_lig, z_ref=ref.z_ref, plane=ref.plane,
            angle_lig_deg=angle, angle_ref_deg=ref.angle_ref_deg,
            dG_lig=dG, dG_ref=ref.dG_ref,
            centroid_lig=centroid, centroid_ref=ref.centroid_ref,
            confidence_weight=weight,
        )
        M, M_loc, M_hyd, M_ori, m_neutral = membrane.M, membrane.M_loc, membrane.M_hyd, membrane.M_ori, membrane.neutral
        m_details = membrane.details
        m_confidence = str(getattr(ref.plane, "confidence", "high"))
    else:
        M = M_loc = M_hyd = M_ori = 1.0
        m_neutral = True
        m_details = {}
        m_confidence = "none"
    rmsd = perturbed["rmsd"]
    plif_proxy = float(np.clip(1.0 - rmsd / 6.0, 0.0, 1.0) * ref.buried_frac)
    return {
        "pose_id": f"{ref.pose_id}__{perturbed['kind']}",
        "ref_pose_id": ref.pose_id,
        "pdb_id": ref.pdb_id,
        "system_id": ref.target,
        "target": ref.target,
        "ligand_class": ref.ligand_class,
        "perturbation": perturbed["kind"],
        "rmsd": float(rmsd),
        "plif_recovery": plif_proxy,
        "posebusters_pass": bool(ref.pb_pass_rate >= 1.0),
        "pb_pass_rate": float(ref.pb_pass_rate),
        "strain_kcal": ref.strain_kcal,
        "V": float(ref.V),
        "M": float(M),
        "M_loc": float(M_loc),
        "M_hyd": float(M_hyd),
        "M_ori": float(M_ori),
        "S": 1.0,
        "m_neutral": bool(m_neutral),
        "m_confidence": m_confidence,
        "s_neutral": True,
        "z_lig_mean": float(np.mean(z_lig)),
        "z_ref_mean": float(np.mean(ref.z_ref)),
        "angle_lig_deg": float(angle),
        "angle_ref_deg": float(ref.angle_ref_deg),
        "dG_lig": float(dG),
        "dG_ref": float(ref.dG_ref),
        "lateral_deviation_angstrom": m_details.get("lateral_deviation_angstrom"),
        "compartment_lig": m_details.get("compartment_lig"),
        "compartment_ref": m_details.get("compartment_ref"),
        "plane_source": ref.plane_source,
        "pb_backend": ref.pb_backend,
    }


def _reference_row(ref: OracleReference) -> Dict[str, Any]:
    return {
        "pose_id": ref.pose_id + "__reference",
        "ref_pose_id": ref.pose_id,
        "pdb_id": ref.pdb_id,
        "system_id": ref.target,
        "target": ref.target,
        "ligand_class": ref.ligand_class,
        "perturbation": "reference",
        "rmsd": 0.0,
        "plif_recovery": float(ref.buried_frac),
        "posebusters_pass": bool(ref.pb_pass_rate >= 1.0),
        "pb_pass_rate": float(ref.pb_pass_rate),
        "strain_kcal": ref.strain_kcal,
        "V": float(ref.V),
        "M": 1.0,
        "M_loc": 1.0,
        "M_hyd": 1.0,
        "M_ori": 1.0,
        "S": 1.0,
        "m_neutral": ref.plane is None,
        "m_confidence": (str(getattr(ref.plane, "confidence", "high")) if ref.plane is not None else "none"),
        "s_neutral": True,
        "z_lig_mean": float(np.mean(ref.z_ref)),
        "z_ref_mean": float(np.mean(ref.z_ref)),
        "angle_lig_deg": float(ref.angle_ref_deg),
        "angle_ref_deg": float(ref.angle_ref_deg),
        "dG_lig": float(ref.dG_ref),
        "dG_ref": float(ref.dG_ref),
        "lateral_deviation_angstrom": 0.0,
        "compartment_lig": ref.plane.compartment(float(np.mean(ref.z_ref))) if ref.plane else "unknown",
        "compartment_ref": ref.plane.compartment(float(np.mean(ref.z_ref))) if ref.plane else "unknown",
        "plane_source": ref.plane_source,
        "pb_backend": ref.pb_backend,
    }


def _m_sanity(frame: pd.DataFrame, refs: Sequence["OracleReference"], cfg: MAESConfig) -> Dict[str, Any]:
    """Reference membrane occupancy + perturbation monotonicity checks."""
    interface = cfg.scores.interface_width_angstrom
    compartments: Dict[str, int] = {}
    in_membrane: List[bool] = []
    sources: Dict[str, int] = {}
    confidences: Dict[str, int] = {}
    for ref in refs:
        sources[ref.plane_source] = sources.get(ref.plane_source, 0) + 1
        if ref.plane is None:
            compartments["unavailable"] = compartments.get("unavailable", 0) + 1
            confidences["none"] = confidences.get("none", 0) + 1
            continue
        confidences[str(getattr(ref.plane, "confidence", "high"))] = (
            confidences.get(str(getattr(ref.plane, "confidence", "high")), 0) + 1
        )
        z_mean = float(np.mean(ref.z_ref))
        compartment = str(ref.plane.compartment(z_mean, interface_width=interface))
        compartments[compartment] = compartments.get(compartment, 0) + 1
        in_membrane.append(bool(abs(z_mean) <= ref.plane.half_thickness + interface))
    by_kind = frame.groupby("perturbation")["M_loc"].agg(["mean", "median", "count"]).round(4).to_dict(orient="index")
    order = ["in_pocket", "near_miss", "wrong_pocket", "out_of_membrane"]
    means = [by_kind.get(k, {}).get("mean", float("nan")) for k in order]
    finite = [m for m in means if np.isfinite(m)]
    monotone = bool(all(finite[i] >= finite[i + 1] - 0.02 for i in range(len(finite) - 1)))
    oom = frame[frame["perturbation"] == "out_of_membrane"]
    oom_valid = oom[~oom["m_neutral"].astype(bool)] if "m_neutral" in oom.columns else oom
    m_by_conf = (
        frame[frame["perturbation"] != "reference"]
        .groupby("m_confidence")["M"]
        .agg(["mean", "median", "count"])
        .round(4)
        .to_dict(orient="index")
        if "m_confidence" in frame.columns
        else {}
    )
    return {
        "n_reference": int(len(refs)),
        "reference_compartments": compartments,
        "reference_in_membrane_fraction": float(np.mean(in_membrane)) if in_membrane else float("nan"),
        "M_loc_by_perturbation": by_kind,
        "M_by_confidence_perturbed": m_by_conf,
        "perturbation_order": order,
        "M_loc_monotone": monotone,
        "out_of_membrane_M_loc_lt_0.2_fraction": (
            float((oom_valid["M_loc"] < 0.2).mean()) if len(oom_valid) else float("nan")
        ),
        "out_of_membrane_n_evaluated": int(len(oom_valid)),
        "plane_sources": sources,
        "plane_confidences": confidences,
        "confidence_grading": "high=1.0, medium=0.5, none=0.0 (M shrunk toward neutral)",
        "m_neutral_fraction_reference": float(
            np.mean([ref.plane is None for ref in refs])
        ) if refs else float("nan"),
    }


# ------------------------------------------------------- state descriptor GMM
MAIN_STATE_FEATURES = ("tm_angle_mean_deg", "pocket_void_A3", "chain_interface_min_A")
NBD_FEATURE = "nbd_centroid_dist_A"


def state_cluster_report(
    manifest: pd.DataFrame,
    states_path: str | Path = DEFAULT_STATES,
    n_components: int = 2,
    seed: int = 42,
) -> Dict[str, Any]:
    """GMM clustering on TM angle + pocket volume + chain interface (C-line W3).

    NBD centroid distance is reported separately for the 6 ABC-transporter rows
    (MsbA 5 + LptB 1) and NaN is NEVER imputed to 0 (that would fabricate a state).
    """
    from sklearn.mixture import GaussianMixture

    states_path = Path(states_path)
    if not states_path.exists():
        return {"available": False}
    states = pd.read_csv(states_path)
    main = manifest.copy()
    if "system_id" not in main.columns:
        from .api import attach_pose_ids

        main = attach_pose_ids(main)
    needed = ["pdb_id", "ligand_id"] + list(MAIN_STATE_FEATURES) + [NBD_FEATURE, "ligand_buried_frac"]
    # merge on pdb_id + ligand_id only: the same ligand can appear in several
    # chains/copies while the descriptor is a per-pose quantity
    merged = main.merge(states[needed], on=["pdb_id", "ligand_id"], how="left", suffixes=("", "_state"))
    features = merged[list(MAIN_STATE_FEATURES)].apply(pd.to_numeric, errors="coerce")
    coverage = {c: int(features[c].notna().sum()) for c in MAIN_STATE_FEATURES}
    imputed = features.fillna(features.median())
    mean = imputed.mean()
    std = imputed.std(ddof=0).replace(0, 1.0)
    x = ((imputed - mean) / std).to_numpy(dtype=float)
    k = int(min(max(2, n_components), max(2, min(4, len(merged)))))
    gmm = GaussianMixture(n_components=k, covariance_type="full", random_state=seed, n_init=10).fit(x)
    labels = gmm.predict(x)
    sizes = pd.Series(labels).value_counts().sort_index().to_dict()
    centers = pd.DataFrame(gmm.means_, columns=list(MAIN_STATE_FEATURES))
    nbd = merged.loc[merged[NBD_FEATURE].notna(), ["target", "pdb_id", "ligand_id", NBD_FEATURE]]
    return {
        "available": True,
        "feature_set": list(MAIN_STATE_FEATURES),
        "coverage_non_na": coverage,
        "n_rows": int(len(merged)),
        "n_components": k,
        "cluster_sizes": {str(kk): int(vv) for kk, vv in sizes.items()},
        "cluster_centers_z": centers.round(3).to_dict(orient="records"),
        "nbd_abc_sublayer": {
            "n_rows": int(len(nbd)),
            "note": "NBD centroid only for ABC transporters (MsbA 5 + LptB 1); NaN not imputed to 0",
            "values": nbd.to_dict(orient="records"),
        },
        "imputation": "median imputation for GMM only; raw coverage reported above",
    }


# --------------------------------------------------------------------- n_eff
def n_eff_block(frame: pd.DataFrame, manifest: pd.DataFrame, cfg: MAESConfig) -> Dict[str, Any]:
    """target x CF cluster count + Kish + ICC-corrected n_eff (Lead W3 rule)."""
    from .api import attach_pose_ids, cluster_keys, select_analysis_set
    from .evaluation import cluster_bootstrap_mean, design_effect_n_eff

    attach = attach_pose_ids(manifest)
    primary = select_analysis_set(attach, "primary_main")
    keys_primary = cluster_keys(primary, "target_crystal_form")
    counts = keys_primary.value_counts()
    kish_primary = float((counts.sum() ** 2) / np.square(counts.to_numpy(dtype=float)).sum())

    pdb_to_cf = dict(zip(primary["pdb_id"].astype(str), primary["cluster_target_crystal_form"].astype(str)))
    cf = frame["pdb_id"].astype(str).map(pdb_to_cf).fillna("unknown").to_numpy()
    if {"rmsd", "posebusters_pass"} <= set(frame.columns):
        accepted = (
            (pd.to_numeric(frame["rmsd"], errors="coerce").to_numpy(dtype=float) < 2.0)
            & frame["posebusters_pass"].fillna(False).astype(bool).to_numpy()
        ).astype(float)
    else:
        accepted = (frame["MAES"] >= 0.5).astype(float).to_numpy()
    icc_binary = design_effect_n_eff(accepted, cf)
    icc_continuous = design_effect_n_eff(frame["MAES"].to_numpy(dtype=float), cf)
    boot = cluster_bootstrap_mean(accepted, cf, n_boot=2000, seed=cfg.scores.random_state)

    # pose-level aggregation (the correct unit for the 40-pose primary set):
    # one value per reference pose, clustered by target x CF
    pose = (
        frame.assign(_accepted=accepted)
        .groupby("ref_pose_id")
        .agg(mean_maes=("MAES", "mean"), accept_frac=("_accepted", "mean"), pdb_id=("pdb_id", "first"))
        .reset_index()
    )
    pose_cf = pose["pdb_id"].astype(str).map(pdb_to_cf).fillna("unknown").to_numpy()
    icc_pose_maes = design_effect_n_eff(pose["mean_maes"].to_numpy(dtype=float), pose_cf)
    icc_pose_accept = design_effect_n_eff(pose["accept_frac"].to_numpy(dtype=float), pose_cf)
    icc_pose_maes["n_eff_40_pose_projection"] = (
        40.0 / icc_pose_maes["design_effect"] if icc_pose_maes["design_effect"] else float("nan")
    )
    icc_pose_accept["n_eff_40_pose_projection"] = (
        40.0 / icc_pose_accept["design_effect"] if icc_pose_accept["design_effect"] else float("nan")
    )

    reference = frame[frame["perturbation"] == "reference"]
    ref_cf = reference["pdb_id"].astype(str).map(pdb_to_cf).fillna("unknown").to_numpy()
    ref_icc = design_effect_n_eff(np.ones(len(reference)), ref_cf)
    return {
        "n_obs_primary": int(len(primary)),
        "n_clusters_primary": int(keys_primary.nunique()),
        "cluster_definition": "target_x_crystal_form (space_group | cell)",
        "kish_n_eff_primary": kish_primary,
        "pose_level_icc_maes": icc_pose_maes,
        "pose_level_icc_accept_frac": icc_pose_accept,
        "oracle_harness_icc_binary_acceptable": icc_binary,
        "oracle_harness_icc_continuous_maes": icc_continuous,
        "oracle_harness_acceptance_bootstrap": boot,
        "reference_only_all_acceptable": ref_icc,
        "headline": (
            f"n_obs=40, n_clusters={int(keys_primary.nunique())}, "
            f"n_eff~{icc_pose_maes['n_eff_40_pose_projection']:.1f} "
            "(target x CF cluster bootstrap, ICC from oracle harness)"
        ),
        "note": (
            "n_eff = n_obs / (1 + (m_bar - 1) * ICC). ICC is estimated on the oracle harness "
            "at the pose level (mean MAES per reference pose, 40 poses x 11 CF) and MUST be "
            "recomputed on B-line predictions. CF01 (17 of 40 inhibitor poses) dominates the "
            "design effect; the cluster count 11 is the upper bound, the ICC projection is the "
            "number to quote in the text."
        ),
    }


# ---------------------------------------------------------------------- run
def run_oracle(
    manifest_path: str | Path = DEFAULT_MANIFEST,
    planes_path: str | Path = DEFAULT_PLANES,
    struct_dir: str | Path = DEFAULT_STRUCT_DIR,
    states_path: Optional[str | Path] = DEFAULT_STATES,
    subset: str = "primary_main",
    n_perturb: int = 6,
    seed: int = 42,
    figure_path: Optional[str | Path] = None,
    report_path: Optional[str | Path] = None,
) -> Dict[str, Any]:
    """Run the oracle validation and return a JSON-serialisable report."""
    cfg = default_config()
    rng = np.random.default_rng(seed)
    try:  # keep RDKit force-field chatter out of the report logs
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")
    except Exception:  # pragma: no cover
        pass
    refs, skipped = load_reference_poses(
        manifest_path, planes_path, struct_dir, states_path, subset=subset, cfg=cfg
    )
    if not refs:
        raise RuntimeError("no usable reference poses (check coord files / membrane planes)")

    rows: List[Dict[str, Any]] = []
    for ref in refs:
        rows.append(_reference_row(ref))
        kinds = [PERTURBATION_KINDS[i % len(PERTURBATION_KINDS)] for i in range(n_perturb)]
        for kind in kinds:
            perturbed = perturb_pose(ref, kind, rng, interface_width=cfg.scores.interface_width_angstrom)
            rows.append(_score_perturbed(ref, perturbed, cfg))
    frame = pd.DataFrame(rows)

    model = MAESModel(cfg).fit(frame)
    scores = model.predict(frame)
    eval_primary = (frame["rmsd"].to_numpy(dtype=float) < cfg.labels.eval_rmsd).astype(int)
    calib = ((frame["rmsd"].to_numpy(dtype=float) < cfg.labels.calib_rmsd)
             & frame["posebusters_pass"].to_numpy(dtype=bool)).astype(int)
    metrics = {
        "n_rows": int(len(frame)),
        "n_references": int(len(refs)),
        "n_positive_calibration": int(calib.sum()),
        "n_positive_eval_primary": int(eval_primary.sum()),
        "auroc_primary": auroc(eval_primary, scores),
        "brier_primary": brier_score(eval_primary, scores),
        "ece_primary": ece(eval_primary, scores),
        "rho_maes_rmsd": spearman_rho(scores, frame["rmsd"].to_numpy(dtype=float)),
        "maes_range": [float(scores.min()), float(scores.max())],
    }
    try:
        ablation = ablation_study(frame, cfg).to_dict(orient="records")
    except Exception as exc:  # pragma: no cover
        ablation = [{"error": f"{type(exc).__name__}: {exc}"}]

    # divergence cases (checklist C5)
    divergence = frame.assign(MAES=scores)
    high_maes_bad_rmsd = divergence[(divergence["MAES"] >= 0.5) & (divergence["rmsd"] > 2.0)]
    low_maes_good_rmsd = divergence[(divergence["MAES"] < 0.5) & (divergence["rmsd"] <= 1.0)]
    examples = []
    for label, sub in (("MAES_high_RMSD_bad", high_maes_bad_rmsd), ("MAES_low_RMSD_good", low_maes_good_rmsd)):
        for _, row in sub.head(5).iterrows():
            examples.append(
                {
                    "class": label,
                    "pose_id": row["ref_pose_id"],
                    "perturbation": row["perturbation"],
                    "rmsd": float(row["rmsd"]),
                    "MAES": float(row["MAES"]),
                    "M_loc": float(row["M_loc"]),
                    "M_ori": float(row["M_ori"]),
                }
            )
    m_sanity = _m_sanity(frame, refs, cfg)
    manifest_frame = pd.read_csv(Path(manifest_path))
    state_clusters = state_cluster_report(manifest_frame, states_path or DEFAULT_STATES, seed=seed)
    try:
        neff = n_eff_block(frame.assign(MAES=scores), manifest_frame, cfg)
    except Exception as exc:  # pragma: no cover - defensive
        neff = {"error": f"{type(exc).__name__}: {exc}"}
    report: Dict[str, Any] = {
        "maes_version": cfg.version,
        "mode": "oracle-reference-pose + rigid-body perturbations (pipeline validation only)",
        "subset": subset,
        "skipped": skipped,
        "metrics": metrics,
        "ablation": ablation,
        "divergence_counts": {
            "MAES_high_RMSD_bad": int(len(high_maes_bad_rmsd)),
            "MAES_low_RMSD_good": int(len(low_maes_good_rmsd)),
        },
        "divergence_examples": examples,
        "m_sanity": m_sanity,
        "state_clusters": state_clusters,
        "n_eff": neff,
        "harness_notes": [
            "V is constant under rigid-body perturbation (PoseBusters checks are internal "
            "to the ligand), so the oracle harness validates the pipeline, not scientific "
            "discrimination.",
            "S is neutral (1.0) here: a rigid ligand transform does not change the protein "
            "state; S will be exercised once B-line provides predicted receptor/state arms.",
            "|rho(MAES, RMSD)| > 0.95 is expected in this harness because M is driven by the "
            "same translation magnitude as RMSD; on real predictions the band 0.5-0.8 is the "
            "target and rho_rmsd_proxy_warning will flag a proxy-like MAES.",
            "Reference poses are M=1 by construction (M compares against the reference); the "
            "M sanity claim is the monotone decrease across perturbation classes.",
        ],
        "model_metadata": model.metadata,
        "label_fingerprints": {
            "calibration": label_fingerprint(calib),
            "evaluation_primary": label_fingerprint(eval_primary),
        },
    }
    if figure_path is not None:
        path = Path(figure_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        _oracle_figure(frame, scores, cfg, path)
        report["figure"] = str(path)
    if report_path is not None:
        path = Path(report_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
        report["report_path"] = str(path)
    return report


def _oracle_figure(frame: pd.DataFrame, scores: np.ndarray, cfg: MAESConfig, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frame = frame.assign(MAES=scores)
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8), dpi=200)
    order = ["reference", "in_pocket", "near_miss", "wrong_pocket", "out_of_membrane"]
    order = [k for k in order if k in set(frame["perturbation"])]
    data = [frame.loc[frame["perturbation"] == k, "M_loc"].to_numpy() for k in order]
    axes[0].boxplot(data, labels=order, showfliers=False)
    axes[0].set_ylabel("M_loc")
    axes[0].set_title("M sanity: membrane location")
    axes[0].tick_params(axis="x", rotation=30, labelsize=7)

    colors = {"reference": "black", "in_pocket": "tab:green", "near_miss": "tab:olive",
              "wrong_pocket": "tab:orange", "out_of_membrane": "tab:red"}
    for kind in order:
        sub = frame[frame["perturbation"] == kind]
        axes[1].scatter(sub["rmsd"], sub["MAES"], s=16, alpha=0.8, label=kind, c=colors.get(kind))
    axes[1].axvline(cfg.labels.eval_rmsd, ls=":", c="grey", lw=1)
    axes[1].axvline(cfg.labels.calib_rmsd, ls=":", c="black", lw=1)
    axes[1].set_xlabel("RMSD (A)")
    axes[1].set_ylabel("MAES")
    axes[1].set_title("oracle harness: MAES vs RMSD")
    axes[1].legend(fontsize=6, frameon=False)

    y_primary = (frame["rmsd"].to_numpy(dtype=float) < cfg.labels.eval_rmsd).astype(int)
    curve = reliability_curve(y_primary, scores, n_bins=8)
    valid = curve["count"] > 0
    axes[2].plot([0, 1], [0, 1], ls="--", lw=1, c="grey")
    axes[2].plot(curve.loc[valid, "mean_pred"], curve.loc[valid, "frac_pos"], "o-", ms=4)
    axes[2].set_xlabel("predicted P(acceptable)")
    axes[2].set_ylabel("observed fraction")
    axes[2].set_title("reliability (oracle)")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
