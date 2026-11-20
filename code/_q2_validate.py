import json
from pathlib import Path

root = Path(".")
manuscript = root / "manuscript"
main = (manuscript / "main_skeleton.md").read_text(encoding="utf-8")
disc = (manuscript / "discussion_draft.md").read_text(encoding="utf-8")
trackb = (manuscript / "results_trackb_draft.md").read_text(encoding="utf-8")
methods = (manuscript / "methods_maes_draft.md").read_text(encoding="utf-8")
report = json.loads((root / "results/maes/w12/w12_report_formal.json").read_text(encoding="utf-8"))
pm = report["primary"]["metrics"]
pw = report["primary"]["wild_cluster_bootstrap"]
pr = report["primary"]["rho_cluster_bootstrap"]

banned = ["calibrated", "reusable", "method ranking", "external validation", "首个", "首次"]
all_text = "\n".join(p.read_text(encoding="utf-8") for p in sorted(manuscript.glob("*.md")))
banned_hits = {term: term.lower() in all_text.lower() for term in banned}
old_w12 = {
    "0.492": "0.492" in all_text,
    "|ρ| 0.026": "|ρ| 0.026" in all_text,
    "Brier 0.137": "Brier 0.137" in all_text,
    "ECE 0.131": "ECE 0.131" in all_text,
}

checks = {
    "q2_title": "Benchmarking structure-based virtual screening across Gram-negative antibacterial targets: target-specific enrichment, transfer axes and failure modes" in main,
    "q2_abstract_sentence": "We evaluate a membrane-aware pose-quality score (MAES) as a pre-registered exploratory arm" in main,
    "q2_contributions_trackB": "40 inhibitor cognate-redocking poses across 9 membrane targets" in main,
    "q2_contributions_dockings": "28,136 successfully docked ligands" in main,
    "q2_2_6_heading": "### 2.6 MAES as an exploratory pose-quality score: calibration gate failed (Fig 6)" in main,
    "q2_2_6_final": "AUROC 0.496" in main and "|ρ(MAES,RMSD)| 0.075" in main and "19 eval positives" in main,
    "q2_fig6_negative": "Calibration-gate negative result for MAES" in main,
    "q2_ledger_rejected": "REJECTED — W12 formal NO-GO" in main and "probability/ranking" in main,
    "q2_conclusion": "Across a two-regime Gram-negative panel, enrichment is target-specific" in disc and "AC2 2/9 = 22.2%" in disc and "MAES calibration (AUROC 0.496, |ρ| 0.075)" in disc,
    "w12_final_metrics": abs(pm["auroc_primary"] - 0.4956485702445089) < 1e-9 and abs(pm["rho_maes_rmsd"] - (-0.0745122157034245)) < 1e-9,
    "w12_final_ci": abs(pw["ci_lo"] - 0.3182415270647716) < 1e-9 and abs(pr["ci_hi"] - 0.16148198458621854) < 1e-9,
    "w12_manifest_sha": report.get("manifest_sha256") == "e4ae5c23461474f773a016ec415e9da8e5dd6ff10bbef55eda73c726cb4fa3ca",
    "w12_no_go": report.get("verdict") == "NO_GO" and report.get("n_unique_pose_keys") == 146 and report.get("n_duplicates") == 0,
    "banned_terms_absent": not any(banned_hits.values()),
    "old_w12_numbers_absent": not any(old_w12.values()),
    "trackb_w12_final": "AUROC **0.496**" in trackb and "|ρ(MAES,RMSD)| **0.075**" in trackb,
    "methods_w12_final": "AUROC **0.496**" in methods and "|ρ(MAES,RMSD)| **0.075**" in methods,
    "discussion_w12_final": "AUROC **0.496**" in disc and "|ρ(MAES,RMSD)| **0.075**" in disc,
}
checks = {key: bool(value) for key, value in checks.items()}
out = {
    "all_passed": bool(all(checks.values())),
    "verdict": "NO_GO",
    "checks": checks,
    "banned_hits": banned_hits,
    "old_w12_hits": old_w12,
    "report": str(root / "results/maes/w12/w12_report_formal.json"),
    "figure": str(root / "figures/fig6_w12_calibration_formal.png"),
    "note": "Q2 reframing applied; final W12 manifest sha e4ae5c23...; numbers 0.496/0.075/0.145/0.124.",
}
(root / "figures/w12_validation.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
(root / "figures/task36_formal_validation.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print(json.dumps(out, indent=2))
