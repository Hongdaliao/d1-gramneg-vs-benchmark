import json
p=json.load(open(r"figures/trackA_enrichment_full.json",encoding="utf-8"))
for e in p["entries"]:
    ef1 = "n.r." if e.get("ef1_percent") is None else "%.2f" % e["ef1_percent"]
    print(e["target"], "|", e["panel_role"], "|", e["n_active"], "| %.3f (%.3f-%.3f)" % (e["auc"], e["auc_ci_lo"], e["auc_ci_hi"]), "| %.2f" % e["ef5_percent"], "|", ef1, "| clipped", e["n_box_clipped"])
print("metal_vs_nonmetal", json.dumps(p.get("metal_vs_nonmetal"), indent=1))
print("by_metal_stratum", json.dumps(p.get("by_metal_stratum"), indent=1))
