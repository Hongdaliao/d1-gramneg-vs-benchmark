"""Summarize final method-axis analysis: pooled paired deltas and S29-S31 tables."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(r"D:\harness\projects\d1_gramneg_vs")
OUT = ROOT / "results/trackA/method_axis/analysis"
TARGETS = ["LpxC_Ec", "LpxC_Pa", "DHFR_Ec", "Gyrase_Ec", "GyrA_Ec", "GyrB_Ec", "FabI_Ec", "TEM1_Ec"]
COMMON7 = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec", "GyrA_Ec", "TEM1_Ec"]
N_BOOT, SEED = 2000, 42
methods = ["diffdock_main", "boltz2", "drugclip", "unidock"]

metrics = pd.read_csv(OUT / "metrics_all_methods_final.csv")
paired = pd.read_csv(OUT / "paired_all_methods_final.csv")
cov = pd.read_csv(OUT / "coverage_all_methods_final.csv")
pooled = json.loads((OUT / "metrics_pooled_final.json").read_text(encoding="utf-8"))

# pooled paired deltas, target-level bootstrap
rows = []
for m in methods:
    d = paired[paired["method"] == m].copy()
    allowed = COMMON7 if m == "unidock" else TARGETS
    d = d[d["target"].isin(allowed)]
    arr = d["delta_vs_vina"].to_numpy(float)
    rng = np.random.default_rng(SEED)
    draws = [float(np.mean(arr[rng.integers(0, len(arr), len(arr))])) for _ in range(N_BOOT)]
    rows.append({"method": m, "scope": "Uni-Dock 7 targets (excl FabI)" if m == "unidock" else "main 8 targets",
                 "mean_delta_vs_vina": float(np.mean(arr)),
                 "ci_lo": float(np.percentile(draws, 2.5)), "ci_hi": float(np.percentile(draws, 97.5)),
                 "n_targets": len(arr)})
pd.DataFrame(rows).to_csv(OUT / "table_S31d_paired_pooled.csv", index=False)
(OUT / "paired_pooled_final.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")

# S29/S30 per-target + pooled summary
for m, name in (("drugclip", "S29_drugclip"), ("unidock", "S30_unidock")):
    d = metrics[metrics["method"] == m].copy()
    # LepB descriptive row from coverage (not in metrics)
    tgt = "LepB_Ec"
    if not ((d["target"] == tgt).any()):
        d = pd.concat([d, pd.DataFrame([{"target": tgt, "method": m, "n_ligands": 6, "n_analysis": 6,
                                          "n_active": 3, "n_decoy": 3, "auc": np.nan, "auc_ci_lo": np.nan,
                                          "auc_ci_hi": np.nan, "top5_recall": np.nan, "rho_vina": np.nan,
                                          "status": "descriptive only (n=6)"}])], ignore_index=True)
    d.to_csv(OUT / f"table_{name}_per_target.csv", index=False)
    p = pooled["methods"][m]["pooled"]
    pd.DataFrame([{"method": m, "scope": "Uni-Dock 7 targets (excl FabI)" if m == "unidock" else "main 8 targets",
                   "mean_auc": p["mean_auc"], "ci_lo": p["ci95"][0], "ci_hi": p["ci95"][1],
                   "n_targets": p["n_targets"], "n_ok": int(cov[cov["method"] == m]["ok"].sum()),
                   "n_na": int(cov[cov["method"] == m]["na"].sum())}]).to_csv(OUT / f"table_{name}_pooled.csv", index=False)

# remove zero Vina rows from paired table
paired[paired["method"] != "vina"].to_csv(OUT / "table_S31b_method_paired.csv", index=False)
cov.to_csv(OUT / "table_S31c_method_coverage.csv", index=False)
print("wrote final summaries")
print(pd.read_csv(OUT / "table_S31d_paired_pooled.csv").to_string(index=False))
