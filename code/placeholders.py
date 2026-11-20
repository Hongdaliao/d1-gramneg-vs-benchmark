"""Three synthetic placeholder systems used to smoke-test MAES end to end.

The Track B pose table (co-crystal ligands, ~120 poses) is not ready yet, so
v1.0 uses three *placeholder bodies* with known generative structure:

===========  =====================================================  ==================
system       regime                                               membrane info
===========  =====================================================  ==================
PL1_MsbA     docking-like, medium difficulty                       OPM-like plane
PL2_LptB     co-folding-like, harder (more failures)               OPM-like plane
PL3_soluble  soluble-target control (M must degrade to neutral)     none
===========  =====================================================  ==================

Each placeholder row carries the *primitive* columns consumed by
``subscores.build_subscores`` plus ground-truth labels (``rmsd``,
``plif_recovery``, ``posebusters_pass``), so the same code path is exercised as
with real Track B data.  Swap in the real table with
``api.build_features`` when Track B lands; only the column contract matters.
"""

from __future__ import annotations

import json
import zlib
from typing import Dict, List

import numpy as np
import pandas as pd

from .membrane_map import MembranePlane, TableMembraneMap

__all__ = [
    "PLACEHOLDER_SPECS",
    "load_placeholders",
    "make_placeholder",
    "placeholder_membrane_map",
]

PLACEHOLDER_SPECS: Dict[str, Dict[str, object]] = {
    "PL1_MsbA": {
        "description": "docking-like, medium difficulty, ABC transporter pocket (membrane available)",
        "membrane_available": True,
        "q_beta": (2.0, 2.0),
        "method_cycle": ("vina", "boltz2", "drugclip"),
        "phi_ref": (1.20, 4.50, 3.00, 0.80),
    },
    "PL2_LptB": {
        "description": "co-folding-like, harder, LptB2FG-like pocket (membrane available)",
        "membrane_available": True,
        "q_beta": (1.4, 3.0),
        "method_cycle": ("boltz2", "af3", "vina"),
        "phi_ref": (1.05, 5.80, 2.40, 1.10),
    },
    "PL3_soluble": {
        "description": "soluble-target control: no membrane plane -> M must be neutral",
        "membrane_available": False,
        "q_beta": (2.0, 2.0),
        "method_cycle": ("vina", "drugclip", "gnina"),
        "phi_ref": (0.90, 3.10, 4.20, 0.50),
    },
}


def _clamp01(x: np.ndarray) -> np.ndarray:
    return np.clip(x, 0.0, 1.0)


def make_placeholder(
    system_id: str,
    n_poses: int = 24,
    seed: int = 42,
) -> pd.DataFrame:
    """Generate one placeholder body with primitive columns + ground-truth labels."""
    if system_id not in PLACEHOLDER_SPECS:
        raise KeyError(f"unknown placeholder {system_id!r}; choose from {sorted(PLACEHOLDER_SPECS)}")
    spec = PLACEHOLDER_SPECS[system_id]
    # deterministic per-system seed (built-in hash() is salted per process)
    rng = np.random.default_rng(seed + int(zlib.crc32(system_id.encode("utf-8")) % 10_000))
    n = int(n_poses)
    q = rng.beta(float(spec["q_beta"][0]), float(spec["q_beta"][1]), size=n)  # quality in [0,1]
    failures = 1.0 - q

    # --- ground truth labels
    rmsd = np.clip(0.35 + 6.5 * failures ** 2 + rng.normal(0, 0.25, n), 0.05, None)
    plif_recovery = _clamp01(0.15 + 0.9 * q + rng.normal(0, 0.12, n))

    # --- V primitives
    pb_pass_rate = _clamp01(0.50 + 0.50 * q + rng.normal(0, 0.10, n))
    posebusters_pass = rng.random(n) < pb_pass_rate
    strain_kcal = np.clip(0.4 + 11.0 * failures ** 1.5 + rng.normal(0, 1.2, n), 0.0, None)

    # --- M primitives (only meaningful when a plane exists)
    membrane_available = bool(spec["membrane_available"])
    z_ref = rng.uniform(-11.0, 11.0, n)          # reference ligand centre inside the bilayer
    z_lig = z_ref + failures * 7.5 + rng.normal(0, 0.8, n)
    angle_ref = rng.uniform(20.0, 70.0, n)
    angle_lig = angle_ref + failures * 45.0 + rng.normal(0, 4.0, n)
    dG_ref = rng.uniform(-2.0, 2.0, n)
    dG_lig = dG_ref - 6.0 * failures + rng.normal(0, 0.7, n)
    if not membrane_available:
        z_ref = np.full(n, np.nan)
        z_lig = np.full(n, np.nan)
        angle_ref = np.full(n, np.nan)
        angle_lig = np.full(n, np.nan)
        dG_ref = np.full(n, np.nan)
        dG_lig = np.full(n, np.nan)

    # --- S primitives (normalised descriptors; Track B should z-score over refs)
    phi_ref = np.asarray(spec["phi_ref"], dtype=float)
    phi_lig = [
        (phi_ref + failures[i] * np.array([0.60, 1.20, 0.90, 0.40]) + rng.normal(0, 0.08, 4)).tolist()
        for i in range(n)
    ]

    method_cycle: tuple = spec["method_cycle"]
    methods = [method_cycle[i % len(method_cycle)] for i in range(n)]

    frame = pd.DataFrame(
        {
            "pose_id": [f"{system_id}_pose{i:03d}" for i in range(n)],
            "system_id": system_id,
            "method": methods,
            "ligand_id": [f"LIG{i % 7:02d}" for i in range(n)],
            "membrane_available": membrane_available,
            "pb_pass_rate": pb_pass_rate,
            "posebusters_pass": posebusters_pass,
            "strain_kcal": strain_kcal,
            "z_lig_mean": z_lig,
            "z_ref_mean": z_ref,
            "angle_lig_deg": angle_lig,
            "angle_ref_deg": angle_ref,
            "dG_lig": dG_lig,
            "dG_ref": dG_ref,
            "phi_lig": [json.dumps(vec) for vec in phi_lig],
            "phi_ref": [json.dumps(phi_ref.tolist()) for _ in range(n)],
            "rmsd": rmsd,
            "plif_recovery": plif_recovery,
        }
    )
    return frame


def load_placeholders(n_poses: int = 24, seed: int = 42) -> pd.DataFrame:
    """Concatenate the three placeholder bodies (deterministic given the seed)."""
    frames: List[pd.DataFrame] = []
    for offset, system_id in enumerate(sorted(PLACEHOLDER_SPECS)):
        frames.append(make_placeholder(system_id, n_poses=n_poses, seed=seed + offset))
    return pd.concat(frames, ignore_index=True)


def placeholder_membrane_map() -> TableMembraneMap:
    """A plane table that covers PL1/PL2 only (PL3 exercises the neutral path)."""
    planes = {
        "PL1_MsbA": MembranePlane(
            normal=(0.0, 0.0, 1.0), center=(0.0, 0.0, 0.0), thickness=30.0, source="placeholder-opm"
        ),
        "PL2_LptB": MembranePlane(
            normal=(0.05, 0.0, 0.9975), center=(0.0, 0.0, 1.5), thickness=28.0, source="placeholder-opm"
        ),
    }
    return TableMembraneMap(planes)
