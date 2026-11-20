"""Build v3 working SI draft tables S18-S24 (no manuscript_v2 edits)."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "manuscript/v3_working/si"
OUT.mkdir(parents=True, exist_ok=True)

# S18 ligand class review
lig = pd.read_csv(ROOT / "data/trackB/ligand_class_review.csv")
lig.to_csv(OUT / "table_S18_ligand_class_review.csv", index=False)

# S19 structure-source all-valid rows, annotated with QC subset
src = pd.read_csv(ROOT / "results/trackA/axes/structure_source_allvalid.csv")
src["qc_subset"] = src["excluded_predicted"].map(lambda v: "QC-passed" if str(v).strip().lower() == "false" else "QC-excluded")
src.to_csv(OUT / "table_S19_structure_source_qc.csv", index=False)

# S20 cross-species bootstrap corrected
cross = json.loads((ROOT / "results/trackA/axes/cross_species_bootstrap_corrected.json").read_text(encoding="utf-8"))
cross_rows = [
    {"item": "point_estimate", "auc_ec": cross["point"]["auc_ec"], "auc_pa": cross["point"]["auc_pa"], "delta_auc": cross["point"]["delta_auc_pa_minus_ec"], "ci_lo": None, "ci_hi": None, "method": "point estimates"},
    {"item": "reviewer_decoy_level", "auc_ec": None, "auc_pa": None, "delta_auc": cross["point"]["delta_auc_pa_minus_ec"], "ci_lo": cross.get("reviewer_authoritative", {}).get("decoy_level_delta_ci95", [None, None])[0], "ci_hi": cross.get("reviewer_authoritative", {}).get("decoy_level_delta_ci95", [None, None])[1], "method": "decoy-level bootstrap with replacement (reviewer authoritative)"},
    {"item": "reviewer_active_level", "auc_ec": None, "auc_pa": None, "delta_auc": cross["point"]["delta_auc_pa_minus_ec"], "ci_lo": cross.get("reviewer_authoritative", {}).get("active_level_delta_ci95", [None, None])[0], "ci_hi": cross.get("reviewer_authoritative", {}).get("active_level_delta_ci95", [None, None])[1], "method": "active-level bootstrap with replacement (reviewer authoritative)"},
    {"item": "old_ci_artifact", "auc_ec": None, "auc_pa": None, "delta_auc": cross["point"]["delta_auc_pa_minus_ec"], "ci_lo": cross.get("old_cross_species_csv_ci", [None, None])[0], "ci_hi": cross.get("old_cross_species_csv_ci", [None, None])[1], "method": "old isin()-dedup artifact; do not use"},
]
pd.DataFrame(cross_rows).to_csv(OUT / "table_S20_cross_species_bootstrap.csv", index=False)

# S21 ECFP4 test-only
ecfp = pd.read_csv(ROOT / "results/trackA/ecfp4_baseline/baseline_metrics_test_only.csv")
ecfp.to_csv(OUT / "table_S21_ecfp4_test_only.csv", index=False)

# S22 W12 complete-case summary
cc = json.loads((ROOT / "results/maes/w12/w12_complete_case.json").read_text(encoding="utf-8"))
rows22 = [
    {"analysis": "full", "n": 146, "n_eval_positive": 19, "auroc": 0.4956, "brier": 0.145, "note": "14 non-assessable rows treated as negative"},
    {"analysis": "complete_case", "n": cc["reviewer_authoritative"]["n"], "n_eval_positive": 19, "auroc": cc["reviewer_authoritative"]["auroc"], "brier": cc["reviewer_authoritative"]["brier"], "note": "pb_pass_rate present; reviewer-authoritative; local reproduction AUROC 0.4332 / Brier 0.1594"},
]
pd.DataFrame(rows22).to_csv(OUT / "table_S22_w12_complete_case.csv", index=False)

# S23 W12 power strata + CV folds
strata = pd.read_csv(ROOT / "results/maes/w12/w12_strata_table.csv")
strata.to_csv(OUT / "table_S23_w12_strata.csv", index=False)
folds = pd.read_csv(ROOT / "results/maes/w12/w12_cv_folds.csv")
folds.to_csv(OUT / "table_S23b_w12_cv_folds.csv", index=False)

# S24 anchor sensitivity
anc = json.loads((ROOT / "results/litpcba/anchor_sensitivity.json").read_text(encoding="utf-8"))
rows24 = [
    {"item": "AC2_primary", "value": "2/9 = 22.2%", "ci_or_note": "pre-registered n_val>=20; FAIL", "source": "anchor_AC.json"},
    {"item": "AC2_include_TP53", "value": "3/10 = 30.0%", "ci_or_note": "Wilson95 reviewer 6.3-54.7 (local two-sided Wilson 10.8-60.3); PASS at threshold only", "source": "anchor_sensitivity.json"},
    {"item": "AC3_population_SD", "value": "28.28%", "ci_or_note": "reported value", "source": "anchor_AC.json"},
    {"item": "AC3_sample_SD", "value": "34.64%", "ci_or_note": "MTORC1, 3 subsample seeds", "source": "anchor_sensitivity.json"},
    {"item": "AC1", "value": "15/15", "ci_or_note": "attempted and completed", "source": "anchor_AC.json"},
]
pd.DataFrame(rows24).to_csv(OUT / "table_S24_anchor_sensitivity.csv", index=False)

index = """# v3 SI table numbering — DRAFT (not final; v2 untouched)

> Working copy: `manuscript/v3_working/manuscript_v3_working.md`; these tables are drafts for the v3 revision, not yet cited in v2.
> v2 final SI remains S1–S17; v3 replaces/adjoints S18–S24 for reanalysis/disclosure.

| v3 table | content | source | status |
|---|---|---|---|
| S18 | Track B ligand-class review (40 poses; current vs suggested class, evidence, confidence) | `data/trackB/ligand_class_review.csv` | draft |
| S19 | Structure-source all-valid table (34 rows; QC-passed vs QC-excluded subset flag) | `results/trackA/axes/structure_source_allvalid.csv` | draft |
| S19-summary | All-valid vs QC-passed metrics (ΔRMSD/ΔPLIF, PB, n) | `results/trackA/axes/structure_source_qc_sensitivity.json` | draft (summary text/JSON) |
| S20 | Cross-species paired comparison, corrected bootstrap (point, decoy-level, active-level, old artifact) | `results/trackA/axes/cross_species_bootstrap_corrected.json` | draft |
| S21 | ECFP4 1-NN test-only baseline per target | `results/trackA/ecfp4_baseline/baseline_metrics_test_only.csv` | draft |
| S22 | W12 full-sample vs complete-case (n=146 vs n=132) | `results/maes/w12/w12_complete_case.json` | draft |
| S23 | W12 per-stratum estimates | `results/maes/w12/w12_strata_table.csv` | draft |
| S23b | W12 CV folds (n_train 116–118) | `results/maes/w12/w12_cv_folds.csv` | draft |
| S24 | LIT-PCBA anchor sensitivity (AC2 2/9 vs 3/10; AC3 population vs sample SD) | `results/litpcba/anchor_sensitivity.json` | draft |

**Caution:** All tables are pending the Lead/A-B decision and a v3 rewrite; do not cite these numbers in v2 or the submission package. Track B S18 shows that only 4/40 poses are true inhibitors; S19 shows 26/34 all-valid pairs are QC-excluded; S21 shows test-only median AUC 0.7745 (not 0.958); S22 shows complete-case AUROC 0.434 vs full 0.496.
"""
(OUT / "SI_v3_numbering_draft.md").write_text(index, encoding="utf-8")
print("wrote draft SI tables to", OUT)
for p in sorted(OUT.iterdir()):
    print(p.name, p.stat().st_size)
