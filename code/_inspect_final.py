import pandas as pd, json
from pathlib import Path
out = Path(r"results/trackA/method_axis/analysis")
print("pooled table")
print(pd.read_csv(out / "table_S31_method_pooled.csv").to_string(index=False))
print("\npaired")
print(pd.read_csv(out / "paired_all_methods_final.csv").to_string(index=False))
m = pd.read_csv(out / "metrics_all_methods_final.csv")
print("\nmetrics")
print(m[m["method"].isin(["drugclip", "unidock"])][["method", "target", "n_analysis", "auc", "auc_ci_lo", "auc_ci_hi", "top5_recall", "rho_vina"]].to_string(index=False))
print("\nfigure", Path(r"figures/v3/fig7_method_axis.png").stat().st_size)
