import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "results/trackA/axes/structure_source.csv"
OUT = ROOT / "results/trackA/axes"

df = pd.read_csv(SRC)
for col in ("holo_rmsd", "afdb_rmsd", "holo_plif", "afdb_plif"):
    df[col] = pd.to_numeric(df[col], errors="coerce")
valid = df[df["afdb_rmsd"].notna() & df["holo_rmsd"].notna()].copy()
allvalid = valid.copy()
qcpassed = valid[valid["excluded_predicted"].astype(str).str.lower() == "false"].copy()
qcpassed.to_csv(OUT / "structure_source_qcpassed.csv", index=False)
allvalid.to_csv(OUT / "structure_source_allvalid.csv", index=False)


def boot_delta(frame, rng, n=2000):
    if len(frame) < 3:
        return [float("nan"), float("nan")]
    vals = []
    for _ in range(n):
        idx = rng.integers(0, len(frame), len(frame))
        sample = frame.iloc[idx]
        vals.append(float(sample["afdb_rmsd"].mean() - sample["holo_rmsd"].mean()))
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]


rng = np.random.default_rng(42)
summary = {
    "source": str(SRC),
    "n_rows_all_valid": int(len(allvalid)),
    "n_rows_qc_passed": int(len(qcpassed)),
    "n_excluded_predicted_true_among_valid": int((allvalid["excluded_predicted"].astype(str).str.lower() == "true").sum()),
    "n_rows_excluded_from_qcpass": int((allvalid["excluded_predicted"].astype(str).str.lower() == "true").sum()),
    "all_valid": {
        "holo_rmsd_mean": float(allvalid["holo_rmsd"].mean()),
        "afdb_rmsd_mean": float(allvalid["afdb_rmsd"].mean()),
        "delta_rmsd_mean": float(allvalid["afdb_rmsd"].mean() - allvalid["holo_rmsd"].mean()),
        "delta_rmsd_ci95": boot_delta(allvalid, rng),
        "holo_plif_mean": float(allvalid["holo_plif"].mean()),
        "afdb_plif_mean": float(allvalid["afdb_plif"].mean()),
        "delta_plif_mean": float(allvalid["afdb_plif"].mean() - allvalid["holo_plif"].mean()),
        "holo_pb_counts": allvalid["holo_pb"].astype(str).value_counts().to_dict(),
        "afdb_pb_counts": allvalid["afdb_pb"].astype(str).value_counts().to_dict(),
    },
    "qc_passed_only": {
        "holo_rmsd_mean": float(qcpassed["holo_rmsd"].mean()),
        "afdb_rmsd_mean": float(qcpassed["afdb_rmsd"].mean()),
        "delta_rmsd_mean": float(qcpassed["afdb_rmsd"].mean() - qcpassed["holo_rmsd"].mean()),
        "delta_rmsd_ci95": boot_delta(qcpassed, rng),
        "holo_plif_mean": float(qcpassed["holo_plif"].mean()),
        "afdb_plif_mean": float(qcpassed["afdb_plif"].mean()),
        "delta_plif_mean": float(qcpassed["afdb_plif"].mean() - qcpassed["holo_plif"].mean()),
        "holo_pb_counts": qcpassed["holo_pb"].astype(str).value_counts().to_dict(),
        "afdb_pb_counts": qcpassed["afdb_pb"].astype(str).value_counts().to_dict(),
    },
    "interpretation": "The manuscript's paired n=34 is the all-valid set, but 26 of those 34 are excluded_predicted=True. QC-passed-only n=8 is the appropriate confirmatory subset; direction in QC-passed subset should be reported separately. Do not present the all-valid geometry/PLIF effect as a clean predicted-receptor effect.",
}
(OUT / "structure_source_qc_sensitivity.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
