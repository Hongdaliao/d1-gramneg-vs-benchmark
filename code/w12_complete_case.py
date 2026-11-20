import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
W12 = ROOT / "results/maes/w12"

# W12 complete-case reproduction
features = pd.read_csv(W12 / "manifest_features.csv")
pred = pd.read_csv(W12 / "w12_predictions_formal.csv")
merged = pred.merge(features[["pose_key", "method", "rmsd", "pb_pass_rate"]],
                    on=["pose_key", "method"], how="left", suffixes=("", "_manifest"))
merged["y"] = pd.to_numeric(merged["rmsd_manifest"], errors="coerce").lt(1).astype(int)
complete = merged[merged["pb_pass_rate"].notna()].copy()
local = {
    "n": int(len(complete)),
    "n_positive": int(complete["y"].sum()),
    "auroc": float(roc_auc_score(complete["y"], complete["MAES_primary_VM"])) if complete["y"].nunique() > 1 else None,
    "brier": float(brier_score_loss(complete["y"], complete["MAES_primary_VM"])),
}
complete_case = {
    "source": "w12_predictions_formal.csv joined to manifest_features.csv on pose_key+method; complete-case = pb_pass_rate not missing (14 rows excluded)",
    "reviewer_authoritative": {"n": 132, "auroc": 0.4341, "brier": 0.1591},
    "local_reproduction": local,
    "note": "The primary 146-row result treats 14 non-assessable rows (pb_pass_rate missing / no reference coordinates) as negative evaluation labels; the complete-case result is materially lower. Report both and state that the AUROC criterion is inconclusive underpowered, while the |rho| criterion is decisive.",
    "missing_rows": merged[merged["pb_pass_rate"].isna()][["pose_key", "method", "rmsd_manifest", "pb_pass_rate"]].to_dict("records"),
}
(W12 / "w12_complete_case.json").write_text(json.dumps(complete_case, indent=2), encoding="utf-8")

# W12 strata table and CV folds already generated; add a compact power JSON
strata = pd.read_csv(W12 / "w12_strata_table.csv")
folds = pd.read_csv(W12 / "w12_cv_folds.csv")
power = {
    "analysis": "primary V+M; full 146 rows",
    "n_observations": 146,
    "n_evaluation_positive": 19,
    "n_calibration_positive": 37,
    "n_clusters": 40,
    "auroc": 0.4956,
    "auroc_wild_ci95": [0.3182, 0.6722],
    "auroc_threshold": 0.65,
    "auroc_ci_includes_threshold": True,
    "complete_case": complete_case,
    "rho": -0.0745,
    "rho_ci95": [-0.3015, 0.1615],
    "rho_abs": 0.0745,
    "rho_abs_threshold": 0.40,
    "rho_ci_upper_below_threshold": True,
    "cv_folds": folds.to_dict("records"),
    "strata": strata.to_dict("records"),
    "underpowered_statement": "The W12 design has 19 evaluation positives; the AUROC 95% CI (0.318-0.672) includes 0.65, so the AUROC criterion is inconclusive rather than evidence of AUROC <0.65. The |rho| criterion is decisive (upper CI 0.161 <0.40). Per-stratum estimates reverse direction across method/receptor groups and are descriptive only.",
}
(W12 / "w12_power_disclosure.json").write_text(json.dumps(power, indent=2), encoding="utf-8")

print("wrote W12 complete-case and power disclosure")
print(json.dumps({"complete_case": local, "n_missing": len(complete_case["missing_rows"])}, indent=2))
