import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
W12 = ROOT / "results/maes/w12"
report = json.loads((W12 / "w12_report_formal.json").read_text(encoding="utf-8"))

rows = []
for stratum_type, entries in report["strata"].items():
    for row in entries:
        rows.append({
            "analysis": "primary_V+M",
            "stratum_type": stratum_type,
            "stratum": row.get("stratum"),
            "n": row.get("n"),
            "auroc": row.get("auroc"),
            "rho_rmsd": row.get("rho_rmsd"),
            "note": row.get("note", ""),
        })
pd.DataFrame(rows).to_csv(W12 / "w12_strata_table.csv", index=False)

folds = []
for row in report["cv"]["folds"]:
    folds.append({
        "fold": row.get("fold"),
        "n_train": row.get("n_train"),
        "n_test": row.get("n_test"),
        "fit_ok": row.get("fit_ok"),
        "error": row.get("error", ""),
    })
pd.DataFrame(folds).to_csv(W12 / "w12_cv_folds.csv", index=False)

print(pd.DataFrame(rows).to_string(index=False))
print(pd.DataFrame(folds).to_string(index=False))
