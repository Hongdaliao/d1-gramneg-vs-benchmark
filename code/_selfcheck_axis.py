import json
from pathlib import Path
import pandas as pd

ROOT = Path(r"D:\harness\projects\d1_gramneg_vs")
OUT = ROOT / "results/trackA/method_axis/analysis"
v3p = ROOT / "manuscript/v3/manuscript_v3.md"
draft = (ROOT / "manuscript/v3/method_axis_integration_draft.md").read_text(encoding="utf-8")
v3 = v3p.read_text(encoding="utf-8")
metrics_json = json.loads((OUT / "metrics_pooled.json").read_text(encoding="utf-8"))
pooled = metrics_json["pooled"]
metrics = pd.read_csv(OUT / "metrics_per_target.csv")
s25 = pd.read_csv(OUT / "table_S25_method_axis_per_target.csv")
s26 = pd.read_csv(OUT / "table_S26_method_axis_paired.csv")
s26b = pd.read_csv(OUT / "table_S26b_boltz2_paired.csv")
s27b = pd.read_csv(OUT / "table_S27b_boltz2_coverage.csv")
s28 = pd.read_csv(OUT / "table_S28_method_axis_singlechain.csv")
cov = pd.read_csv(OUT / "coverage_na.csv")
scna = pd.read_csv(OUT / "singlechain_coverage_na.csv")
figdata = pd.read_csv(OUT / "fig7_data.csv")
overlap = pd.read_csv(OUT / "ecfp4_leakage_overlap.csv")
amendment_sha = __import__("hashlib").sha256((ROOT / "notes/preregistration_amendment12_method_axis_draft.md").read_bytes()).hexdigest()
bz = pooled["boltz2"]

checks = {
    "main 806 expected": len(cov) == 806,
    "main 13 NA": int((cov["status"] != "ok").sum()) == 13,
    "no ECFP4 rows in metrics": int((metrics["method"] == "ECFP4").sum()) == 0,
    "no ECFP4 rows in S25": int((s25["method"] == "ECFP4").sum()) == 0,
    "S25 has DiffDock+Vina+Boltz": int((s25["method"] == "DiffDock").sum()) >= 8 and int((s25["method"] == "Vina").sum()) == 8 and int((s25["method"] == "Boltz-2").sum()) >= 8,
    "S26 paired 8": len(s26) == 8,
    "S26b Boltz paired 8": len(s26b) == 8,
    "S27b Boltz coverage 166": len(s27b) == 166,
    "S28 3 rows": len(s28) == 3,
    "singlechain 4 NA": len(scna) == 4,
    "singlechain 202 ok": int(s28["n_singlechain_ok"].sum()) == 202,
    "fig7 includes singlechain": bool(set(figdata["method"]) & {"DiffDock_singlechain"}),
    "fig7 includes Boltz": bool(set(figdata["method"]) & {"Boltz-2"}),
    "overlap table 9 rows": len(overlap) == 9,
    "amendment sha current": metrics_json["metadata"]["amendment_sha256"] == amendment_sha,
    "Boltz status ok": bz.get("status") == "ok",
    "Boltz 165 ok/1 NA": int(bz.get("n_ok")) == 165 and int(bz.get("n_na")) == 1,
    "Boltz pooled 0.718": abs(float(bz.get("pooled_auc_mean")) - 0.7183333333) < 1e-3,
    "Boltz CI 0.585-0.840": abs(bz["pooled_auc_ci95"][0] - 0.5903958333) < 1e-3 and abs(bz["pooled_auc_ci95"][1] - 0.84003125) < 1e-3,
    "DiffDock pooled 0.528": abs(pooled["diffdock_auc_mean"] - 0.527962) < 1e-4,
    "DiffDock CI 0.479-0.583": abs(pooled["diffdock_auc_ci_lo"] - 0.479109) < 1e-4 and abs(pooled["diffdock_auc_ci_hi"] - 0.582646) < 1e-4,
    "DiffDock delta -0.077": abs(pooled["delta_auc_diffdock_minus_vina_mean"] + 0.076811) < 1e-4,
    "DiffDock delta CI -0.233/+0.081": abs(pooled["delta_auc_ci_lo"] + 0.233125) < 1e-4 and abs(pooled["delta_auc_ci_hi"] - 0.080588) < 1e-4,
    "v3 bytes 57577": v3p.stat().st_size == 57577,
    "v3 sha 2c47ff": __import__("hashlib").sha256(v3p.read_bytes()).hexdigest().startswith("2c47ff"),
    "v3 has Methods 2.7": "### 2.7 Exploratory method-class comparison" in v3,
    "v3 has Results 3.10": "### 3.10 Exploratory method-class comparison" in v3,
    "v3 has Figure 7 caption": "**Figure 7. Exploratory method-class comparison" in v3,
    "v3 has S1-S28": "S1–S28" in v3 and "S26b" in v3,
    "v3 abstract <=250": len(v3.split("## Abstract")[1].split("## Keywords")[0].split()) <= 250,
    "draft no SENS_PENDING": "[SENS_PENDING]" not in draft,
    "draft no BOLTZ_PENDING": "[BOLTZ_PENDING]" not in draft,
    "draft ECFP4 withdrawn": "not recomputed" in draft and "0.927" not in draft and "0.698–1.000" not in draft,
    "draft Boltz 0.718": "0.718" in draft and "0.585–0.840" in draft,
    "draft singlechain GyrA": "0.488 to 0.667" in draft and "+0.178" in draft,
    "draft singlechain FabI": "0.469 to 0.444" in draft and "−0.025" in draft,
    "draft singlechain caveat": "does not change the Vina ranking" in draft,
    "draft limitation phrase": "single-chain sensitivity is reported in S28" in draft,
    "draft subset caveat": "not comparable to full-panel Table 1" in draft,
    "draft no leading claim": "does not support a leading-engine claim" in draft,
}
for k, v in checks.items():
    print(("OK " if v else "FAIL "), k)
print("ALL", all(checks.values()))
(OUT / "SELFCHECK.json").write_text(json.dumps({
    "all_ok": bool(all(checks.values())),
    "v3": {"bytes": v3p.stat().st_size, "mtime_epoch": v3p.stat().st_mtime,
           "sha256": __import__("hashlib").sha256(v3p.read_bytes()).hexdigest()},
    "checks": {k: bool(v) for k, v in checks.items()},
}, indent=2), encoding="utf-8")
print("v3 bytes", v3p.stat().st_size)
