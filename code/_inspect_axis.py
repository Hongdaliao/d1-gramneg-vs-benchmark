import json
from pathlib import Path
import pandas as pd

out = Path(r"results/trackA/method_axis/analysis")
print("files", [p.name for p in sorted(out.iterdir())])
print("\nmetrics")
print(pd.read_csv(out / "metrics_per_target.csv").to_string(index=False))
print("\npaired")
print(pd.read_csv(out / "paired_bootstrap_per_target.csv").to_string(index=False))
cov = pd.read_csv(out / "coverage_na.csv")
print("\ncoverage", cov.groupby("target")["status"].value_counts().to_dict())
print("\nNA rows")
print(cov[cov["status"] != "ok"][["target", "complex_name", "ligand_id", "status", "attempts", "failed_at"]].to_string(index=False))
print("\nsinglechain")
print(pd.read_csv(out / "singlechain_status.csv").to_string(index=False))
print("\nboltz", json.load(open(out / "metrics_pooled.json"))["boltz2"])
