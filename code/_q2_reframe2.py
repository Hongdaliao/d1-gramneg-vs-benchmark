from pathlib import Path

p = Path(r"manuscript/main_skeleton.md")
lines = p.read_text(encoding="utf-8").splitlines()

new_contrib = "4. Contributions (framed as integration, not primacy): (i) a two-regime panel for Gram-negative antibacterial targets (**9-target Track A, 1,748 actives / 28,136 successfully docked ligands**; **Track B: 40 inhibitor cognate-redocking poses across 9 membrane targets**); (ii) per-target quantification of three transfer axes on that panel (temporal, cross-species, structure-source; R6), showing **axis- and target-dependent** transfer rather than a single transferability score; (iii) three reproducible benchmark/tool failure modes — metal coordination, boron atom typing, DUD-E-style decoy bias — separated from target biology; (iv) a negative-result-driven reporting protocol (pre-registered gates, per-target cluster CIs, 2D baselines, no post-hoc threshold/method/data changes), with MAES retained only as an exploratory/descriptive pose-quality analysis."
new_26 = "- MAES is reported as an exploratory/descriptive pose-quality analysis only: the pre-registered W12 calibration gate fails on 146 pose observations (AUROC 0.496, wild-CI 0.318–0.672; |ρ(MAES,RMSD)| 0.075, cluster-CI −0.302–0.161; Brier 0.145, ECE 0.124; 19 eval positives), so we make no probability, ranking or MAES-as-contribution claim, and present the calibration curves as a negative-result diagnostic (full V+M+S exploratory only)."
new_fig = "| Fig 6 | **Calibration-gate negative result for MAES (reliability/Brier/ECE; AUROC 0.496; |ρ| 0.075; 146 poses)** — exploratory/descriptive only, not a probability score or ranking figure | WP4 | **NO_GO; retained as negative result** |"
new_ledger = '| "MAES provides a probability/ranking score" | AUROC/Brier/ECE + ablation + cross-target validation | **REJECTED — W12 formal NO-GO (AUROC 0.496, |ρ| 0.075); descriptive only, Q2 fallback** |'

for i, line in enumerate(lines):
    if line.startswith("4. Contributions (framed as integration"):
        lines[i] = new_contrib
        print("contrib replaced")
    if line.startswith("- Definition; calibration curves") or line.startswith("- **W12 calibration gate (formal"):
        lines[i] = new_26 if line.startswith("- Definition;") else ""
        print("2.6 line handled", i + 1)
    if line.startswith("| Fig 6 | MAES calibration + IRT"):
        lines[i] = new_fig
        print("fig row replaced")
    if line.startswith('| "MAES is calibrated and reusable"'):
        lines[i] = new_ledger
        print("ledger replaced")

# remove any empty line left in 2.6 from old second bullet
cleaned = [line for line in lines if line != ""] if False else lines
p.write_text("\n".join(cleaned) + "\n", encoding="utf-8")
print("written")
