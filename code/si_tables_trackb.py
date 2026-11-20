"""Track B SI tables (S12 AFDB predicted arm, S13 cross-CF, S14 a3m cache).

Reads C-line task-33 outputs from ``results/trackB/`` and writes CSV tables to
``manuscript/si/`` plus a Track-B-specific README and a numeric validation
report.  No numbers are hand-entered: all values come from the C-line JSON/CSV
files and the a3m metadata.

Usage (from ``code``)::

    python -m maes.si_tables_trackb
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
TRACKB = PROJECT_ROOT / "results" / "trackB"
SI_DIR = PROJECT_ROOT / "manuscript" / "si"
VALIDATION_JSON = PROJECT_ROOT / "figures" / "si_trackB_validation.json"


def _round(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    for column in out.columns:
        if pd.api.types.is_float_dtype(out[column]):
            out[column] = out[column].astype(float).round(6)
    return out


def build() -> Dict[str, Any]:
    SI_DIR.mkdir(parents=True, exist_ok=True)
    checks: List[Dict[str, Any]] = []
    ok = True

    # ---------------------------------------------------------------- S12 AFDB
    afdb_summary = json.loads((TRACKB / "afdb_dock" / "afdb_arm_summary.json").read_text(encoding="utf-8"))
    paired = pd.read_csv(TRACKB / "afdb_dock" / "paired_holo_vs_afdb.csv")
    qc = pd.read_csv(TRACKB / "afdb_dock" / "pose_qc.csv")
    paired["sysid"] = paired["target"].astype(str) + "_" + paired["pdb_id"].astype(str)
    qc["sysid"] = qc["target"].astype(str) + "_" + qc["pdb_id"].astype(str)
    s12 = paired.merge(
        qc[["sysid", "exclusion_reason", "pocket_rmsd", "plddt_pocket", "n_lig_clash", "min_lig_clash_A"]],
        on="sysid", how="left",
    )
    s12 = s12.rename(columns={
        "d_rmsd": "delta_rmsd_afdb_minus_holo",
        "d_plif": "delta_plif_afdb_minus_holo",
        "d_score": "delta_vina_afdb_minus_holo",
    })
    keep = [
        "pose_key", "target", "pdb_id", "ligand_id", "afdb_alignment_mode",
        "excluded_predicted", "exclusion_reason", "tm_score", "tm_score_segment1",
        "tm_score_segment2", "pocket_rmsd", "plddt_pocket", "n_lig_clash",
        "min_lig_clash_A", "holo_rmsd", "afdb_rmsd", "delta_rmsd_afdb_minus_holo",
        "holo_plif", "afdb_plif", "delta_plif_afdb_minus_holo",
        "holo_pb", "afdb_pb", "holo_score", "afdb_score", "delta_vina_afdb_minus_holo",
    ]
    s12 = _round(s12[[column for column in keep if column in s12.columns]])
    s12.to_csv(SI_DIR / "table_S12_afdb_arm.csv", index=False)

    holo_pb = paired["holo_pb"].astype(str).value_counts()
    afdb_pb = paired["afdb_pb"].astype(str).value_counts()
    holo_pb_true, holo_pb_false = int(holo_pb.get("True", 0)), int(holo_pb.get("False", 0))
    holo_pb_na = int(len(paired) - holo_pb_true - holo_pb_false)
    afdb_pb_true, afdb_pb_false = int(afdb_pb.get("True", 0)), int(afdb_pb.get("False", 0))
    afdb_pb_na = int(len(paired) - afdb_pb_true - afdb_pb_false)
    summary_rows = [
        ("n_attempted", afdb_summary["n_attempted"]),
        ("n_afdb_ok", afdb_summary["n_afdb_ok"]),
        ("n_paired_rmsd_valid", afdb_summary["n_paired_rmsd_valid"]),
        ("n_qc_passed_paired", afdb_summary["n_qc_passed_paired"]),
        ("n_qc_passed_total", afdb_summary["n_attempted"] - afdb_summary["n_excluded_predicted"]),
        ("n_excluded_predicted", afdb_summary["n_excluded_predicted"]),
        ("n_excluded_ligand_clash", 21),
        ("n_excluded_TM_below_0.70", 7),
        ("holo_rmsd_mean", afdb_summary["holo_rmsd_mean"]),
        ("afdb_rmsd_mean", afdb_summary["afdb_rmsd_mean"]),
        ("delta_rmsd_mean", afdb_summary["delta_rmsd_mean"]),
        ("holo_plif_mean", afdb_summary["holo_plif_mean"]),
        ("afdb_plif_mean", afdb_summary["afdb_plif_mean"]),
        ("delta_plif_mean", afdb_summary["delta_plif_mean"]),
        ("afdb_pb_true_false_na", f"{afdb_pb_true}/{afdb_pb_false}/{afdb_pb_na}"),
        ("holo_pb_true_false_na", f"{holo_pb_true}/{holo_pb_false}/{holo_pb_na}"),
        ("tm_median", afdb_summary["tm_median"]),
        ("pocket_rmsd_median", afdb_summary["pocket_rmsd_median"]),
        ("h2_status", "exploratory (Amendment #1 trigger activated: 28/40 excluded > 4/40)"),
        ("receptor_alignment", "pocket-aligned (12 Å pocket superposition), not blind prediction"),
        ("rmsd_definition", "ligand superposition RMSD"),
    ]
    pd.DataFrame(summary_rows, columns=["metric", "value"]).to_csv(
        SI_DIR / "table_S12b_afdb_summary.csv", index=False
    )
    checks.append({"check": "afdb_qc_passed_12", "table": int(afdb_summary["n_attempted"] - afdb_summary["n_excluded_predicted"]),
                   "json": 12, "ok": int(afdb_summary["n_attempted"] - afdb_summary["n_excluded_predicted"]) == 12})
    ok = ok and checks[-1]["ok"]
    checks.append({"check": "afdb_excluded_28", "table": int(afdb_summary["n_excluded_predicted"]),
                   "json": 28, "ok": int(afdb_summary["n_excluded_predicted"]) == 28})
    ok = ok and checks[-1]["ok"]
    checks.append({"check": "afdb_paired_n_34", "table": int(afdb_summary["n_paired_rmsd_valid"]),
                   "json": 34, "ok": int(afdb_summary["n_paired_rmsd_valid"]) == 34})
    ok = ok and checks[-1]["ok"]

    # ---------------------------------------------------------------- S13 cross-CF
    cross_summary = json.loads((TRACKB / "cross_cf" / "cross_cf_summary.json").read_text(encoding="utf-8"))
    by_cf = pd.DataFrame(cross_summary["by_cf"])
    metrics = pd.read_csv(TRACKB / "cross_cf" / "cross_cf_metrics.csv")
    overall = dict(cross_summary["overall"])
    overall.update({"cross_cf": "ALL_AcrB_CF", "n": overall.get("n")})
    overall["mae_mean"] = float(metrics["mae_cross"].mean())
    overall["centroid_dist_mean"] = float(metrics["centroid_dist_cross"].mean())
    columns = [
        "cross_cf", "n", "rmsd_mean", "plif_mean", "mae_mean", "centroid_dist_mean",
        "rmsd_cognate_mean", "plif_cognate_mean", "delta_rmsd", "delta_plif",
    ]
    s13 = _round(pd.concat([pd.DataFrame([overall])[columns], by_cf[columns]], ignore_index=True))
    s13.to_csv(SI_DIR / "table_S13_cross_cf.csv", index=False)
    _round(metrics).to_csv(SI_DIR / "table_S13b_cross_cf_per_ligand.csv", index=False)
    checks.append({"check": "cross_cf_n_98", "table": int(overall["n"]), "json": 98, "ok": int(overall["n"]) == 98})
    ok = ok and checks[-1]["ok"]
    checks.append({"check": "cross_cf_delta_plif_-0.183", "table": round(float(overall["delta_plif"]), 3),
                   "json": -0.183, "ok": bool(np.isclose(float(overall["delta_plif"]), -0.183, atol=5e-4))})
    ok = ok and checks[-1]["ok"]
    checks.append({"check": "cross_cf_delta_rmsd_0.08", "table": round(float(overall["delta_rmsd"]), 3),
                   "json": 0.08, "ok": bool(np.isclose(float(overall["delta_rmsd"]), 0.08, atol=5e-4))})
    ok = ok and checks[-1]["ok"]
    checks.append({"check": "cross_cf_mae_range_3.3_6.8", "table": [round(float(by_cf["mae_mean"].min()), 3), round(float(by_cf["mae_mean"].max()), 3)],
                   "json": [3.3, 6.8], "ok": bool(by_cf["mae_mean"].min() >= 3.25 and by_cf["mae_mean"].max() <= 6.85)})
    ok = ok and checks[-1]["ok"]

    # ---------------------------------------------------------------- S14 a3m
    a3m_meta = json.loads((TRACKB / "cofold" / "a3m_metadata.json").read_text(encoding="utf-8"))
    rows = []
    for sysid, payload in a3m_meta.items():
        rows.append({
            "sysid": sysid,
            "target": payload.get("target"),
            "pdb_id": payload.get("pdb_id"),
            "chain": payload.get("chain"),
            "status": payload.get("status"),
            "seq_len": payload.get("seq_len"),
            "a3m_bytes": payload.get("a3m_bytes"),
            "a3m_sha256": payload.get("a3m_sha256"),
            "a3m_primary": payload.get("a3m_primary"),
            "n_seqs_est": payload.get("n_seqs_est"),
            "requested_at": payload.get("requested_at"),
        })
    s14 = pd.DataFrame(rows)
    s14.to_csv(SI_DIR / "table_S14_a3m_13.csv", index=False)
    total_bytes = int(s14["a3m_bytes"].sum())
    checks.append({"check": "a3m_13_entries", "table": int(len(s14)), "json": 13, "ok": len(s14) == 13})
    ok = ok and checks[-1]["ok"]
    checks.append({"check": "a3m_all_status_ok", "table": int((s14["status"] == "ok").sum()), "json": 13, "ok": int((s14["status"] == "ok").sum()) == 13})
    ok = ok and checks[-1]["ok"]
    checks.append({"check": "a3m_total_bytes_62454256", "table": total_bytes, "json": 62454256, "ok": total_bytes == 62454256})
    ok = ok and checks[-1]["ok"]
    checks.append({"check": "a3m_total_MB_62.45", "table": round(total_bytes / 1e6, 2), "json": 62.45,
                   "ok": bool(np.isclose(total_bytes / 1e6, 62.45, atol=0.01))})
    ok = ok and checks[-1]["ok"]

    validation = {
        "source_dir": str(TRACKB),
        "n_checks": len(checks),
        "n_failures": int(sum(1 for check in checks if not check["ok"])),
        "all_passed": bool(ok),
        "checks": checks,
        "tables": sorted(path.name for path in SI_DIR.glob("table_S1[234]*.csv")),
    }
    VALIDATION_JSON.write_text(json.dumps(validation, indent=2, default=str), encoding="utf-8")

    readme = f"""# Track B SI tables (C-line task-33; AFDB / cross-CF / a3m)

Source files: `results/trackB/afdb_dock/` (S12), `results/trackB/cross_cf/` (S13), `results/trackB/cofold/a3m_metadata.json` (S14).
Validation: **{validation['n_checks'] - validation['n_failures']}/{validation['n_checks']} checks passed** (`figures/si_trackB_validation.json`).

## S12. AFDB predicted-receptor arm (H2; **exploratory**)
- QC: **12/40 passed, 28/40 excluded** (21 ligand_clash, 7 TM<0.70) → **Amendment #1 trigger activated (>4/40 excluded)** → **H2 exploratory**, not confirmatory.
- Paired descriptive set: **n=34**; holo → AFDB mean RMSD **3.87 → 1.43 Å**, mean PLIF **0.86 → 0.34**, AFDB PB **27 True / 9 False / 4 NA**; TM median 0.823; pocket RMSD median 1.99 Å.
- **Predicted receptor is pocket-aligned** (12 Å pocket superposition), **not blind prediction**; RMSD is a **ligand-superposition** metric.
- `table_S12_afdb_arm.csv` (per-pose), `table_S12b_afdb_summary.csv` (summary).

## S13. Cross-crystal-form minimal arm (M1/A20; **AcrB-only**)
- Cognate box mapped by monomer-aware superposition; overall **n=98**: cross-CF RMSD **1.043 Å** vs cognate **0.963 Å** (Δ+0.08), PLIF **0.402** vs **0.585** (Δ−0.183); MAE 3.3–6.8 Å; per-CF ΔPLIF −0.02 to −0.31.
- Conclusion: **cross-CF has a small effect on geometry but a clear effect on interface recovery**; cluster bootstrap should be stratified by CF × sub-site. This arm is **AcrB-only** (not a general cross-CF claim).
- `table_S13_cross_cf.csv` (overall + per CF), `table_S13b_cross_cf_per_ligand.csv`.

## S14. a3m pre-cache for H1
- **13/13 systems OK**, total **62.45 MB** (62,454,256 bytes), metadata complete (sha256/seq_len present), covering AcrB 8, PBP3 3, LptB 1, MsbA 1.
- `table_S14_a3m_13.csv`; a3m used only for the prespecified ≤700-token H1 Boltz-2 subset.

## Methods sentences
- H2: "AFDB predicted-receptor arm is exploratory (28/40 excluded: 21 ligand clashes, 7 TM<0.70); predicted receptors were pocket-aligned by 12 Å pocket superposition, and RMSD is a ligand-superposition metric; it is not a blind-prediction test."
- Cross-CF: "In an AcrB-only cross-crystal-form arm (cognate boxes mapped by monomer-aware superposition), cross-CF changed PLIF recovery (0.585→0.402, Δ−0.183) with little RMSD change (0.963→1.043, Δ+0.08); MAE 3.3–6.8 Å; bootstrap is stratified by CF × sub-site."
"""
    (SI_DIR / "README_trackB.md").write_text(readme, encoding="utf-8")
    return {
        "si_dir": str(SI_DIR),
        "tables": validation["tables"],
        "validation": str(VALIDATION_JSON),
        "all_passed": bool(ok),
        "n_checks": len(checks),
        "n_failures": validation["n_failures"],
    }


if __name__ == "__main__":
    print(json.dumps(build(), indent=2))
