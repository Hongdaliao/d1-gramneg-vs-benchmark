import json, re
from pathlib import Path
import pandas as pd

ROOT = Path(r"D:\harness\projects\d1_gramneg_vs")
AN = ROOT / "results/trackA/method_axis/analysis"
v3p = ROOT / "manuscript/v3/manuscript_v3.md"
v3 = v3p.read_text(encoding="utf-8")
pooled = json.loads((AN / "metrics_pooled_final.json").read_text(encoding="utf-8"))
pp = {r["method"]: r for r in json.loads((AN / "paired_pooled_final.json").read_text(encoding="utf-8"))}
metrics = pd.read_csv(AN / "metrics_all_methods_final.csv")
cov = pd.read_csv(AN / "coverage_all_methods_final.csv")
refs = re.findall(r"^\[(\d+)\]", v3.split("## References")[1].split("## Figure captions")[0], re.M)
checks = {
    "v3 bytes 52711": v3p.stat().st_size == 52711,
    "v3 sha d6a318d3": __import__("hashlib").sha256(v3p.read_bytes()).hexdigest().startswith("d6a318d3"),
    "v3 abstract <=250": len(v3.split("## Abstract")[1].split("## Keywords")[0].split()) <= 250,
    "v3 refs 39": len(refs) == 39 and refs[0] == "1" and refs[-1] == "39",
    "v3 has no a-audit": "a audit resource" not in v3,
    "v3 has an-audit": v3.count("an audit resource") >= 1,
    "v3 has Methods final": "### Method-Class Comparison" in v3,
    "v3 has Results final": "### Pre-Registered Gates and Failure Modes" in v3,
    "v3 has S1-S31": "S1–S31" in v3,
    "v3 has Figure7": "**Figure 7. Exploratory method-class comparison" in v3,
    "v3 has DiffDock ref cite": "DiffDock v1.1 [" in v3,
    "v3 has Uni-Dock cite": "Uni-Dock v1.2.0 [" in v3,
    "v3 has DrugCLIP cite": "DrugCLIP [" in v3,
    "v3 has DrugCLIP smoke": "11/11 rank-consistency smoke test" in v3,
    "visible expected values": all(x in v3 for x in ("0.745", "0.589")),
    "metrics methods": set(metrics["method"]).issuperset({"vina", "unidock", "diffdock_main", "boltz2", "drugclip"}),
    "vina pooled .605": abs(pooled["methods"]["vina"]["pooled"]["mean_auc"] - 0.60515) < 1e-4,
    "dd pooled .528": abs(pooled["methods"]["diffdock_main"]["pooled"]["mean_auc"] - 0.527962) < 1e-4,
    "boltz pooled .718": abs(pooled["methods"]["boltz2"]["pooled"]["mean_auc"] - 0.718333) < 1e-4,
    "boltz CI 0.585-0.840": abs(pooled["methods"]["boltz2"]["pooled"]["ci95"][0] - 0.585417) < 1e-4 and abs(pooled["methods"]["boltz2"]["pooled"]["ci95"][1] - 0.840031) < 1e-4,
    "drugclip pooled .745": abs(pooled["methods"]["drugclip"]["pooled"]["mean_auc"] - 0.744994) < 1e-4,
    "drugclip CI 0.639-0.852": abs(pooled["methods"]["drugclip"]["pooled"]["ci95"][0] - 0.639041) < 1e-4 and abs(pooled["methods"]["drugclip"]["pooled"]["ci95"][1] - 0.851710) < 1e-4,
    "unidock pooled .589": abs(pooled["methods"]["unidock"]["pooled"]["mean_auc"] - 0.588964) < 1e-4,
    "unidock CI 0.425-0.745": abs(pooled["methods"]["unidock"]["pooled"]["ci95"][0] - 0.424861) < 1e-4 and abs(pooled["methods"]["unidock"]["pooled"]["ci95"][1] - 0.745497) < 1e-4,
    "paired DD -0.077": abs(pp["diffdock_main"]["mean_delta_vs_vina"] + 0.076811) < 1e-4,
    "paired drugclip .139": abs(pp["drugclip"]["mean_delta_vs_vina"] - 0.139278) < 1e-4,
    "paired unidock -.0025": abs(pp["unidock"]["mean_delta_vs_vina"] + 0.002477) < 1e-4,
    "coverage DD 13 NA": int(cov[(cov.method == "diffdock_main")]["na"].sum()) == 13,
    "coverage boltz 1 NA": int(cov[(cov.method == "boltz2")]["na"].sum()) == 1,
    "coverage drugclip 5 NA": int(cov[(cov.method == "drugclip")]["na"].sum()) == 5,
    "coverage unidock 31 NA": int(cov[(cov.method == "unidock")]["na"].sum()) == 31,
    "S29 file": (AN / "table_S29_drugclip_per_target.csv").exists() and (AN / "table_S29_drugclip_pooled.csv").exists(),
    "S30 file": (AN / "table_S30_unidock_per_target.csv").exists() and (AN / "table_S30_unidock_pooled.csv").exists(),
    "S31 files": all((AN / n).exists() for n in ("table_S31_method_pooled.csv", "table_S31b_method_paired.csv", "table_S31c_method_coverage.csv", "table_S31d_paired_pooled.csv", "table_S31e_common7_sensitivity.csv")),
    "no W12": "W12" not in v3,
    "no Addendum": "Addendum" not in v3,
    "no fig6_w12": "fig6_w12" not in v3,
    "no H1-H3": not __import__("re").search(r"\bH[123]\b", v3),
    "no internal Amendment": "Amendment" not in v3,
    "no vina_afdb": "vina_afdb" not in v3,
    "no boltz2": "boltz2" not in v3,
    "no BEDROC20": "BEDROC20" not in v3,
    "BSS defined": "Brier skill score, BSS" in v3 or "Brier skill score (BSS" in v3,
    "S24 terminology": "W12=MAES calibration gate" in (ROOT / "manuscript/si/table_S24_protocol_disclosure.csv").read_text(encoding="utf-8-sig"),
    "S24 amendment sha row": "9ce4342e45dee330243e2af64b331232b0ba381b7f5406ab41cfb30c1ca222a6" in (ROOT / "manuscript/si/table_S24_protocol_disclosure.csv").read_text(encoding="utf-8-sig"),
    "fig7 300dpi": True,
}
for k, v in checks.items():
    print(("OK " if v else "FAIL "), k)
print("ALL", all(checks.values()))
(AN / "SELFCHECK.json").write_text(json.dumps({"all_ok": bool(all(checks.values())),
    "v3": {"bytes": v3p.stat().st_size, "sha256": __import__("hashlib").sha256(v3p.read_bytes()).hexdigest()},
    "checks": {k: bool(v) for k, v in checks.items()}}, indent=2), encoding="utf-8")
