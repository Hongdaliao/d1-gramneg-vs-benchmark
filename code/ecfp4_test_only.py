import json
from pathlib import Path

import numpy as np
import pandas as pd

from maes.evaluation import auroc, enrichment_factor
from maes.tracka import load_score_table, _normalize_ligand_key

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/trackA/ecfp4_baseline"
TRACKA = ROOT / "results/trackA/full"
TARGETS = ["LpxC_Ec", "LpxC_Pa", "DHFR_Ec", "Gyrase_Ec", "GyrA_Ec", "GyrB_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]

rows = []
for target in TARGETS:
    b = pd.read_csv(BASE / f"{target}.csv")
    scores = load_score_table(TRACKA / target / f"{target}_scores.csv")
    clipped = scores["box_clipped"].map(lambda v: str(v).strip().lower() in {"true", "1", "1.0", "yes"})
    scores = scores[~clipped][["ligand_id_norm", "temporal_split", "label"]]
    b["ligand_id_norm"] = b["ligand_id"].map(_normalize_ligand_key)
    frame = b.merge(scores, on="ligand_id_norm", how="left", suffixes=("", "_score"))
    valid = frame[np.isfinite(frame["score"].astype(float))]
    test = valid[valid["temporal_split"].astype(str).str.lower() == "test"]
    dec = valid[valid["label"].astype(int) == 0]
    act = test[test["label"] == 1]
    if len(act) == 0 or len(dec) == 0:
        rows.append({"target": target, "n_test_act": len(act), "n_decoy": len(dec), "auc_test_only": float("nan")})
        continue
    y = np.r_[np.ones(len(act)), np.zeros(len(dec))]
    s = np.r_[act["score"].to_numpy(dtype=float), dec["score"].to_numpy(dtype=float)]
    rows.append({
        "target": target, "n_test_act": int(len(act)), "n_decoy": int(len(dec)),
        "auc_test_only": float(auroc(y, s)),
        "ef5_test_only": float(enrichment_factor(y, s, 0.05)),
    })
out = pd.DataFrame(rows)
out.to_csv(BASE / "baseline_metrics_test_only.csv", index=False)
summary = {
    "method": "ECFP4 1-NN score (max Tanimoto to temporal-train actives) evaluated only on temporal test actives vs all decoys; the published baseline_metrics.csv included train actives in evaluation (leakage).",
    "targets": rows,
    "median_auc_test_only": float(out["auc_test_only"].median()),
    "range_auc_test_only": [float(out["auc_test_only"].min()), float(out["auc_test_only"].max())],
    "n_targets_below_0.80": int((out["auc_test_only"] < 0.80).sum()),
}
(BASE / "baseline_metrics_test_only.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(out.to_string(index=False))
print(json.dumps({k: v for k, v in summary.items() if k != "targets"}, indent=2))
