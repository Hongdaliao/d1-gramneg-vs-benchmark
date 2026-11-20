from pathlib import Path
import re

p = Path(r"code/maes/build_jcim_style_draft.py")
t = p.read_text(encoding="utf-8")

new_abs = """Gram-negative antibacterial discovery faces a dual data problem: many targets have few potent public actives, whereas membrane machines that are structurally tractable seldom provide inhibitor-bound co-crystal poses. Structure-based virtual screening (SBVS) is therefore used under exactly the conditions in which it is least validated. Deep-learning docking and co-folding methods promise better accuracy, but add training-data dependence, GPU cost, receptor-representation sensitivity and benchmark-leakage risk, and balanced head-to-head comparisons remain scarce. We present a two-regime benchmark of nine Gram-negative targets. Track A is a potency-driven panel (1,748 curated actives; 28,136 docked ligands) evaluated by enrichment. Track B is a 40-pose co-crystal benchmark with expert re-annotation (4 inhibitors, 24 antibiotic substrates, 3 covalent β-lactams, 8 fragments and 1 unresolved); a subset-balanced axis compares Vina, Uni-Dock, DiffDock, Boltz-2 and DrugCLIP on 806 ligands, with a 166-ligand Boltz-2 subset and a 206-complex single-chain arm. Track A enrichment is target-specific (AUC 0.300–0.852), and pooled subset AUC was 0.745 for DrugCLIP, 0.718 for Boltz-2, 0.605 for Vina, 0.589 for Uni-Dock (seven targets) and 0.528 for DiffDock; no engine dominated every target. The test-only ECFP4 baseline reached median 0.775, and 7/8 targets matched or exceeded Vina, showing synthetic-decoy inflation. Pre-registered gates failed or were inconclusive: LIT-PCBA did not robustly meet its criteria, and MAES was NO_GO on rank agreement with underpowered AUROC. These results favor per-target, subset-balanced reporting with coverage, decoy ratios, receptor QC and pre-registered gates, not a single ranking. Data and code will be released."""
t = re.sub(r'abstract = """.*?"""', 'abstract = """' + new_abs + '"""', t, count=1, flags=re.S)

# Data Sets figure callouts
t = t.replace("leaving 1,710 scored actives (1,631 primary) in the analysis table.",
              "leaving 1,710 scored actives (1,631 primary) in the analysis table. Figure 1 shows the data landscape, and Figure 2 summarises the two-regime panel composition.")
# Method citation + Figure 7
t = t.replace("Uni-Dock v1.2.0 under Windows Subsystem for Linux", "Uni-Dock v1.2.0 [19] under Windows Subsystem for Linux")
t = t.replace("Pooled AUC on the balanced subset was 0.745", "Figure 7 summarises the method-class comparison. Pooled AUC on the balanced subset was 0.745")
# Gates figure/citation additions
t = t.replace("The MAES calibration gate used 146 pose observations", "**Figure 6** displays the calibration-gate result. The MAES calibration gate used 146 pose observations")
t = t.replace("Brier score 0.145 and expected calibration error 0.124.", "Brier score 0.145 [25] and expected calibration error 0.124 [26].")
t = t.replace("Three benchmark/tool failure modes recur. Metal coordination", "**Figure 5** shows the failure-mode comparison. Three benchmark/tool failure modes recur. Metal coordination")
t = t.replace("Metal coordination is target-specific rather than a stable penalty:", "Metal coordination is target-specific rather than a stable penalty [34-37]:")
# Discussion Figure4 callout and cross-species expansion
t = t.replace("The results do not support a single reliability ranking.", "**Figure 4** summarises the exploratory axes. The results do not support a single reliability ranking.")
# Methods citations and headings
t = t.replace("per-target p-values are not used for inference.", "per-target p-values are not used for inference, and Benjamini–Hochberg values [22] are not used; reporting follows DOME recommendations [29].")
t = t.replace("**Methods**\n\n", "")
t = re.sub(r'^\*\*(.+?)\.\*\* ', r'### \1\n\n', t, flags=re.M)
# Table3 marker corrected
t = t.replace('marker3 = "\\n\\n**Method-Class Comparison**"', 'marker3 = "\\n\\n### Method-Class Comparison"')
p.write_text(t, encoding="utf-8")
print("patched builder")
