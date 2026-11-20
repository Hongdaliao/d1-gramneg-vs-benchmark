import pandas as pd
from pathlib import Path
out = Path(r"results/trackA/method_axis/analysis")
m = pd.read_csv(out / "metrics_per_target.csv")
print("methods", m["method"].value_counts().to_dict())
print("ECFP4")
print(m[m["method"] == "ECFP4"][["target", "auc", "auc_ci_lo", "auc_ci_hi", "top5_recall"]].to_string(index=False))
print("\nsingle comparison")
print(pd.read_csv(out / "singlechain_comparison.csv").to_string(index=False))
print("\nsingle coverage")
print(pd.read_csv(out / "singlechain_coverage.csv").to_string(index=False))
print("\nS28")
print(pd.read_csv(out / "table_S28_method_axis_singlechain.csv").to_string(index=False))
print("\nfiles", [p.name for p in sorted(out.iterdir())])
