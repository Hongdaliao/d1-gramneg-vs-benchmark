"""W5 task-35 — three-axis transferability tables (R6).

Uses only existing Track A/B results (no new docking):
1. Temporal: full `temporal_split == test` actives vs novel-scaffold subset,
   positives vs the same decoys; AUC/EF5 deltas with paired cluster bootstrap.
2. Cross-species: LpxC_Ec vs LpxC_Pa paired on shared scored decoys; AcrB/MexB
   descriptive placeholder (insufficient matched poses).
3. Structure source: Track B holo vs AFDB predicted paired poses (n=34),
   pocket-aligned / ligand-superposition caveats.

Outputs under ``results/trackA/axes/`` plus a summary JSON.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd

from .evaluation import auroc, enrichment_factor
from .tracka import _load_decoy_cluster_map, _normalize_ligand_key, load_score_table

PROJECT_ROOT = Path(__file__).resolve().parents[2]
TRACKA = PROJECT_ROOT / "results" / "trackA" / "full"
SPLITS = PROJECT_ROOT / "data" / "trackA" / "splits"
OUT = PROJECT_ROOT / "results" / "trackA" / "axes"
TRACKB = PROJECT_ROOT / "results" / "trackB"

TARGETS = [
    "LpxC_Ec", "LpxC_Pa", "DHFR_Ec", "Gyrase_Ec", "GyrA_Ec",
    "GyrB_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec",
]


def _accepted(target: str) -> pd.DataFrame:
    frame = load_score_table(TRACKA / target / f"{target}_scores.csv")
    if "box_clipped" in frame.columns:
        clipped = frame["box_clipped"].map(
            lambda value: str(value).strip().lower() in {"true", "1", "1.0", "yes", "y"}
        )
        frame = frame[~clipped]
    frame = frame[np.isfinite(frame["score"].astype(float))]
    return frame


def _novel_ids(target: str) -> set[str]:
    path = SPLITS / f"{target}_temporal_novel.csv"
    if not path.exists():
        return set()
    frame = pd.read_csv(path)
    column = "compound_chembl_id" if "compound_chembl_id" in frame.columns else frame.columns[0]
    return {_normalize_ligand_key(value) for value in frame[column].astype(str)}


def _clusters(frame: pd.DataFrame, target: str) -> np.ndarray:
    pool = _load_decoy_cluster_map(target)
    active_to_inchikey = pool["active_to_inchikey"] if pool else {}
    by_decoy_inchikey = pool.get("by_decoy_inchikey", {}) if pool else {}
    by_decoy = pool["by_decoy"] if pool else {}
    keys: List[str] = []
    for raw_id, norm, label in zip(frame["ligand_id"], frame["ligand_id_norm"], frame["label"]):
        norm = str(norm)
        if int(label) == 1:
            inchikey = active_to_inchikey.get(_normalize_ligand_key(norm), "")
            keys.append("active:" + str(inchikey) if inchikey else "active:" + norm)
        else:
            raw = str(raw_id)
            decoy_inchikey = raw[4:] if raw.lower().startswith("dec_") else raw
            info = by_decoy_inchikey.get(decoy_inchikey) or by_decoy.get(norm) or {}
            inchikey = info.get("inchikey", "")
            keys.append("active:" + str(inchikey) if inchikey else "decoy:" + norm)
    return np.array(keys, dtype=object)


def _metrics(labels: np.ndarray, scores: np.ndarray) -> Dict[str, float]:
    y = np.asarray(labels, dtype=int)
    s = np.asarray(scores, dtype=float)
    valid = np.isfinite(s)
    y, s = y[valid], s[valid]
    if len(y) == 0 or y.sum() == 0 or y.sum() == len(y):
        return {"auc": float("nan"), "ef5": float("nan"), "n_act": int(y.sum()), "n_dec": int((y == 0).sum())}
    return {
        "auc": float(auroc(y, s)),
        "ef5": float(enrichment_factor(y, s, 0.05)),
        "n_act": int(y.sum()),
        "n_dec": int((y == 0).sum()),
    }


def _paired_bootstrap_temporal(
    frame: pd.DataFrame,
    novel: set[str],
    clusters: np.ndarray,
    n_boot: int = 2000,
    seed: int = 42,
) -> Dict[str, Any]:
    labels = frame["label"].astype(int).to_numpy()
    scores = frame["score"].astype(float).to_numpy()
    norms = frame["ligand_id_norm"].astype(str).to_numpy()
    is_novel = (labels == 1) & np.array([norm in novel for norm in norms])
    unique, inverse = np.unique(clusters, return_inverse=True)
    groups = [np.where(inverse == i)[0] for i in range(len(unique))]
    rng = np.random.default_rng(int(seed))
    draws = {"auc_full": [], "ef5_full": [], "auc_novel": [], "ef5_novel": [], "delta_auc": [], "delta_ef5": []}
    for _ in range(int(n_boot)):
        chosen = rng.integers(0, len(unique), size=len(unique))
        idx = np.concatenate([groups[j] for j in chosen])
        full = _metrics(labels[idx], scores[idx])
        nov_mask = is_novel[idx]
        # keep all decoys + novel actives
        keep = (labels[idx] == 0) | nov_mask
        if keep.sum() < 3:
            continue
        nov = _metrics(labels[idx][keep], scores[idx][keep])
        if not all(np.isfinite([full["auc"], full["ef5"], nov["auc"], nov["ef5"]])):
            continue
        draws["auc_full"].append(full["auc"]); draws["ef5_full"].append(full["ef5"])
        draws["auc_novel"].append(nov["auc"]); draws["ef5_novel"].append(nov["ef5"])
        draws["delta_auc"].append(nov["auc"] - full["auc"])
        draws["delta_ef5"].append(nov["ef5"] - full["ef5"])
    def ci(key: str) -> List[float]:
        values = draws[key]
        if len(values) < 20:
            return [float("nan"), float("nan")]
        return [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]
    return {
        "ci": {key: ci(key) for key in draws},
        "n_boot_nonempty": int(len(draws["auc_full"])),
    }


def temporal_axis() -> tuple[pd.DataFrame, Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for target in TARGETS:
        frame = _accepted(target)
        frame = frame[frame["temporal_split"].astype(str).str.lower().isin({"train", "test", "decoy"})].copy()
        frame = frame[(frame["label"] == 0) | (frame["label"] == 1)].copy()
        test = frame[(frame["label"] == 1) & (frame["temporal_split"].astype(str).str.lower() == "test")]
        decoys = frame[frame["label"] == 0]
        novel = _novel_ids(target) & set(test["ligand_id_norm"].astype(str))
        subset = pd.concat([test, decoys], ignore_index=True)
        clusters = _clusters(subset, target)
        full = _metrics(subset["label"].to_numpy(), subset["score"].to_numpy())
        nov_frame = subset[(subset["label"] == 0) | subset["ligand_id_norm"].astype(str).isin(novel)]
        nov = _metrics(nov_frame["label"].to_numpy(), nov_frame["score"].to_numpy())
        boot = _paired_bootstrap_temporal(subset, novel, clusters)
        ci = boot["ci"]
        rows.append({
            "target": target,
            "n_test": int(len(test)),
            "n_novel": int(len(novel)),
            "n_decoy": int(len(decoys)),
            "auc_full": full["auc"], "ef5_full": full["ef5"],
            "auc_novel": nov["auc"], "ef5_novel": nov["ef5"],
            "delta_auc": nov["auc"] - full["auc"],
            "delta_ef5": nov["ef5"] - full["ef5"],
            "auc_full_ci_lo": ci["auc_full"][0], "auc_full_ci_hi": ci["auc_full"][1],
            "auc_novel_ci_lo": ci["auc_novel"][0], "auc_novel_ci_hi": ci["auc_novel"][1],
            "delta_auc_ci_lo": ci["delta_auc"][0], "delta_auc_ci_hi": ci["delta_auc"][1],
            "delta_ef5_ci_lo": ci["delta_ef5"][0], "delta_ef5_ci_hi": ci["delta_ef5"][1],
            "n_boot_nonempty": boot["n_boot_nonempty"],
            "note": "novel-scaffold subset; decoys shared; cluster bootstrap resamples active/decoys clusters; n_novel<5 descriptive only" if len(novel) < 5 else "novel-scaffold subset; decoys shared; cluster bootstrap",
        })
    out = pd.DataFrame(rows)
    summary = {
        "n_targets": int(len(out)),
        "n_targets_n_novel_ge5": int((out["n_novel"] >= 5).sum()),
        "median_delta_auc": float(out["delta_auc"].median()),
        "median_delta_ef5": float(out["delta_ef5"].median()),
        "n_positive_delta_auc": int((out["delta_auc"] > 0).sum()),
        "n_targets_with_novel": int((out["n_novel"] > 0).sum()),
        "caveat": "n_novel<5 targets (LepB n=1) descriptive only; no confirmatory claim",
    }
    return out, summary


def cross_species_axis() -> tuple[pd.DataFrame, Dict[str, Any]]:
    ec = _accepted("LpxC_Ec")
    pa = _accepted("LpxC_Pa")
    ec_dec = ec[(ec["label"] == 0)].copy()
    pa_dec = pa[(pa["label"] == 0)].copy()
    shared = sorted(set(ec_dec["ligand_id_norm"].astype(str)) & set(pa_dec["ligand_id_norm"].astype(str)))
    ec_shared = ec_dec[ec_dec["ligand_id_norm"].astype(str).isin(shared)]
    pa_shared = pa_dec[pa_dec["ligand_id_norm"].astype(str).isin(shared)]
    ec_act = ec[ec["label"] == 1]
    pa_act = pa[pa["label"] == 1]
    ec_auc = float(auroc(np.r_[np.ones(len(ec_act)), np.zeros(len(ec_shared))], np.r_[ec_act["score"], ec_shared["score"]]))
    pa_auc = float(auroc(np.r_[np.ones(len(pa_act)), np.zeros(len(pa_shared))], np.r_[pa_act["score"], pa_shared["score"]]))
    ec_ef5 = float(enrichment_factor(np.r_[np.ones(len(ec_act)), np.zeros(len(ec_shared))], np.r_[ec_act["score"], ec_shared["score"]], 0.05))
    pa_ef5 = float(enrichment_factor(np.r_[np.ones(len(pa_act)), np.zeros(len(pa_shared))], np.r_[pa_act["score"], pa_shared["score"]], 0.05))
    rng = np.random.default_rng(42)
    deltas, ec_draws, pa_draws = [], [], []
    for _ in range(2000):
        ec_idx = rng.integers(0, len(ec_act), len(ec_act))
        pa_idx = rng.integers(0, len(pa_act), len(pa_act))
        dec_idx = rng.integers(0, len(shared), len(shared))
        dec_keys = [shared[i] for i in dec_idx]
        ec_d = ec_shared[ec_shared["ligand_id_norm"].astype(str).isin(dec_keys)]
        pa_d = pa_shared[pa_shared["ligand_id_norm"].astype(str).isin(dec_keys)]
        if len(ec_d) < 2 or len(pa_d) < 2:
            continue
        y_ec = np.r_[np.ones(len(ec_act)), np.zeros(len(ec_d))]
        s_ec = np.r_[ec_act["score"].to_numpy(), ec_d["score"].to_numpy()]
        y_pa = np.r_[np.ones(len(pa_act)), np.zeros(len(pa_d))]
        s_pa = np.r_[pa_act["score"].to_numpy(), pa_d["score"].to_numpy()]
        a_ec = auroc(y_ec, s_ec); a_pa = auroc(y_pa, s_pa)
        if np.isfinite(a_ec) and np.isfinite(a_pa):
            ec_draws.append(a_ec); pa_draws.append(a_pa); deltas.append(a_pa - a_ec)
    def ci(values: List[float]) -> List[float]:
        return [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))] if len(values) >= 20 else [float("nan"), float("nan")]
    rows = [{
        "axis": "cross_species",
        "comparison": "LpxC_Ec vs LpxC_Pa (shared scored decoys)",
        "available": True,
        "n_shared_scored_decoys": len(shared),
        "n_ec_actives": len(ec_act), "n_pa_actives": len(pa_act),
        "auc_ec": ec_auc, "auc_pa": pa_auc, "delta_auc_pa_minus_ec": pa_auc - ec_auc,
        "ef5_ec": ec_ef5, "ef5_pa": pa_ef5, "delta_ef5_pa_minus_ec": pa_ef5 - ec_ef5,
        "auc_ec_ci_lo": ci(ec_draws)[0], "auc_ec_ci_hi": ci(ec_draws)[1],
        "auc_pa_ci_lo": ci(pa_draws)[0], "auc_pa_ci_hi": ci(pa_draws)[1],
        "delta_auc_ci_lo": ci(deltas)[0], "delta_auc_ci_hi": ci(deltas)[1],
        "note": "same shared-decoy set for both targets; decoys resampled at decoy level; not a full cross-species docking benchmark",
    }, {
        "axis": "cross_species",
        "comparison": "AcrB_Ec vs MexB_Pa",
        "available": False,
        "n_shared_scored_decoys": None,
        "n_ec_actives": None, "n_pa_actives": None,
        "auc_ec": None, "auc_pa": None, "delta_auc_pa_minus_ec": None,
        "ef5_ec": None, "ef5_pa": None, "delta_ef5_pa_minus_ec": None,
        "auc_ec_ci_lo": None, "auc_ec_ci_hi": None,
        "auc_pa_ci_lo": None, "auc_pa_ci_hi": None,
        "delta_auc_ci_lo": None, "delta_auc_ci_hi": None,
        "note": "AcrB/MexB: only Track B cognate poses (AcrB 34, MexB 1); no matched cross-species panel; descriptive/not available",
    }]
    summary = {
        "lpxc_shared_scored_decoys": len(shared),
        "auc_ec": ec_auc, "auc_pa": pa_auc,
        "delta_auc_pa_minus_ec": pa_auc - ec_auc,
        "delta_auc_ci": ci(deltas),
        "acrb_mexb": "descriptive only / not available",
    }
    return pd.DataFrame(rows), summary


def structure_axis() -> tuple[pd.DataFrame, Dict[str, Any]]:
    paired = pd.read_csv(TRACKB / "afdb_dock" / "paired_holo_vs_afdb.csv")
    for column in ("holo_rmsd", "afdb_rmsd", "holo_plif", "afdb_plif", "tm_score", "pocket_rmsd"):
        if column in paired.columns:
            paired[column] = pd.to_numeric(paired[column], errors="coerce")
    valid = paired[paired["afdb_rmsd"].notna() & paired["holo_rmsd"].notna()].copy()
    valid["delta_rmsd_afdb_minus_holo"] = valid["afdb_rmsd"] - valid["holo_rmsd"]
    valid["delta_plif_afdb_minus_holo"] = valid["afdb_plif"] - valid["holo_plif"]
    out = valid[[
        "pose_key", "target", "pdb_id", "ligand_id", "afdb_alignment_mode",
        "excluded_predicted", "tm_score", "pocket_rmsd", "holo_rmsd", "afdb_rmsd",
        "delta_rmsd_afdb_minus_holo", "holo_plif", "afdb_plif", "delta_plif_afdb_minus_holo",
        "holo_pb", "afdb_pb",
    ]].copy()
    rng = np.random.default_rng(42)
    deltas_rmsd, deltas_plif = [], []
    for _ in range(2000):
        idx = rng.integers(0, len(valid), len(valid))
        sample = valid.iloc[idx]
        deltas_rmsd.append(float(sample["afdb_rmsd"].mean() - sample["holo_rmsd"].mean()))
        deltas_plif.append(float(sample["afdb_plif"].mean() - sample["holo_plif"].mean()))
    def ci(values: List[float]) -> List[float]:
        return [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]
    summary = {
        "n_paired": int(len(valid)),
        "holo_rmsd_mean": float(valid["holo_rmsd"].mean()),
        "afdb_rmsd_mean": float(valid["afdb_rmsd"].mean()),
        "delta_rmsd_mean": float(valid["afdb_rmsd"].mean() - valid["holo_rmsd"].mean()),
        "delta_rmsd_ci95": ci(deltas_rmsd),
        "holo_plif_mean": float(valid["holo_plif"].mean()),
        "afdb_plif_mean": float(valid["afdb_plif"].mean()),
        "delta_plif_mean": float(valid["afdb_plif"].mean() - valid["holo_plif"].mean()),
        "delta_plif_ci95": ci(deltas_plif),
        "holo_pb_counts": valid["holo_pb"].astype(str).value_counts().to_dict(),
        "afdb_pb_counts": valid["afdb_pb"].astype(str).value_counts().to_dict(),
        "alignment": "pocket-aligned (12 Å pocket superposition), not blind prediction",
        "rmsd_definition": "ligand superposition RMSD",
        "status": "H2 exploratory (Amendment #1 trigger: 28/40 excluded)",
    }
    return out, summary


def build() -> Dict[str, Any]:
    OUT.mkdir(parents=True, exist_ok=True)
    temporal, temporal_summary = temporal_axis()
    cross, cross_summary = cross_species_axis()
    structure, structure_summary = structure_axis()
    temporal.round(6).to_csv(OUT / "temporal.csv", index=False)
    cross.round(6).to_csv(OUT / "cross_species.csv", index=False)
    structure.round(6).to_csv(OUT / "structure_source.csv", index=False)
    summary = {
        "task": "W5 task-35 three-axis transferability",
        "temporal": temporal_summary,
        "cross_species": cross_summary,
        "structure_source": structure_summary,
        "caveats": [
            "Temporal novel-scaffold subset: n_novel<5 targets descriptive only; no confirmatory claim.",
            "Cross-species: only LpxC Ec/Pa has a matched shared-decoy comparison; AcrB/MexB descriptive/not available.",
            "Structure source: AFDB arm is pocket-aligned (not blind prediction), ligand-superposition RMSD; H2 exploratory after 28/40 exclusions.",
        ],
    }
    (OUT / "axes_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    return {"out_dir": str(OUT), "summary": summary, "n_temporal": len(temporal), "n_structure": len(structure)}


if __name__ == "__main__":
    print(json.dumps(build(), indent=2, default=str))
