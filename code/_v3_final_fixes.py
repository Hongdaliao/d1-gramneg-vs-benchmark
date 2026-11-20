from pathlib import Path

p = Path("manuscript/v3/manuscript_v3.md")
t = p.read_text(encoding="utf-8")

repls = [
    ("exhaustiveness 4 for EF1-enabled targets (LpxC_Ec, LpxC_Pa and DHFR_Ec)",
     "exhaustiveness 4 for EF1-enabled targets (LpxC_Ec, LpxC_Pa, DHFR_Ec and Gyrase_Ec)"),
    ("Uni-Dock [19] was evaluated for throughput but was not deployable on the available Windows host;",
     "Uni-Dock [19] was considered but could not be deployed on the available Windows host;"),
    ("AcrB contributes 34 of 55 non-polymer poses but only 34 inhibitor-class poses.",
     "AcrB contributes 55 non-polymer poses in total, of which 34 are inhibitor-class."),
    ("| Complete-case V+M | 132 | 19 | 0.434 (not bootstrapped; reviewer-authoritative) | — | 0.159 | — |",
     "| Complete-case V+M | 132 | 19 | 0.434 (complete-case sensitivity; no CI reported) | — | 0.159 | — |"),
    ('[34] Zhang X, Wong S E, Lightstone F C. Toward fully automated high performance computing drug discovery: a massively parallel virtual screening pipeline for docking and molecular mechanics/generalized Born surface area rescoring to improve enrichment. J Chem Inf Model, 2014, 54(1): 324–337. doi:10.1021/ci4005145',
     '[34] Ganesh H, Mada S, Aryal S, et al. MBind: a web-based platform for metalloprotein and nonmetalloprotein docking with optional ML docking integration. Molecules, 2026, 31(15): 2703. doi:10.3390/molecules31152703'),
    ("Source: `figures/trackA_enrichment_full.json`, `data/target_panel.csv`, `data/trackB/ligand_class_review.csv`; figure file `figures/v3/fig2_benchmark_design.png` / `.pdf`.",
     "Source: `figures/trackA_enrichment_full.json`, `results/trackA/full/target_metadata.csv`, `data/trackB/ligand_class_review.csv`; figure file `figures/v3/fig2_two_regime_panel.png` / `.pdf`."),
    ("Source: `figures/fig3_trackA_enrichment_full.png`.", "Source: `figures/v3/fig3_trackA_enrichment.png` / `.pdf`."),
    ("**Figure 5. Ligand-class composition and test-only 2D comparison.** Track B poses are classified as 4 inhibitors, 24 antibiotic substrates, 3 covalent β-lactams, 8 fragment-like and 1 unresolved; the right panel compares the test-only ECFP4 1-NN baseline with Vina (median AUC 0.7745; 7/8 comparable targets 2D >= Vina). Source: `data/trackB/ligand_class_review.csv`, `results/trackA/ecfp4_baseline/baseline_metrics_test_only.csv`, `figures/trackA_enrichment_full.json`; figure file `figures/v3/fig5_ligand_class_2d.png` / `.pdf`.",
     "**Figure 5. Three failure modes and benchmark-validity comparison.** Metal-stratified Vina AUC, Track B ligand-class composition (4 inhibitors, 24 antibiotic substrates, 3 covalent β-lactams, 8 fragment-like, 1 unresolved) and the test-only ECFP4 1-NN versus Vina comparison (median AUC 0.7745; 7/8 comparable targets 2D >= Vina). Source: `figures/trackA_enrichment_full.json`, `data/trackB/ligand_class_review.csv`, `results/trackA/ecfp4_baseline/baseline_metrics_test_only.csv`; figure file `figures/v3/fig5_failure_modes.png` / `.pdf`."),
    ("Source: `figures/v3/fig6_w12_calibration_formal.png` / `.pdf`.",
     "Source: `figures/v3/fig6_w12_calibration.png` / `.pdf`."),
    ("**Figure 5** shows the Track B ligand-class composition and the test-only 2D comparison; the three failure modes are described below.",
     "**Figure 5** shows the metal-stratified AUC, the Track B ligand-class composition and the test-only 2D comparison; the three failure modes are described below."),
]
for old, new in repls:
    if old not in t:
        print("MISSING:", old[:120])
    else:
        t = t.replace(old, new)
p.write_text(t, encoding="utf-8")
print("v3 final figure/reference/method fixes applied")
