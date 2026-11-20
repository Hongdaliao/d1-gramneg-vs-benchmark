from pathlib import Path

p = Path("manuscript/manuscript_v2.md")
t = p.read_text(encoding="utf-8")

# M1 bridge + m3 auxiliary footnote after Table 1
needle = "| LepB_Ec | auxiliary | 6 | 0.292 (0.067–0.575) | 3.67 | n.r. |"
assert needle in t
t = t.replace(needle, needle + "\n\nOf 1,748 curated actives (1,648 primary), 38 (Gyrase_Ec 17 + LepB 21) fell outside the docking box and were excluded, leaving 1,710 scored actives (1,631 primary) in Table 1. Auxiliary values (GyrA_Ec and LepB_Ec) are descriptive and are not included in any aggregate or trend.")

# Add references [37],[38] after [36]
needle = "[36] Balieiro A M, Fernandes G M B, Girdhari L, et al. Why computational predictions fall short in subclass B1 metallo-beta-lactamase inhibitor discovery: progress, pitfalls, and translational opportunities. ACS Infect Dis, 2026, 12(11): 3573–3583. doi:10.1021/acsinfecdis.6c00608"
assert needle in t
t = t.replace(needle, needle + """

[37] Lu W, Zhang J, Huang W, et al. DynamicBind: predicting ligand-specific protein–ligand complex structure with a deep equivariant generative model. Nat Commun, 2024, 15(1): 1071. doi:10.1038/s41467-024-45461-2

[38] Cai H, Shen C, Jian T, et al. CarsiDock: a deep learning paradigm for accurate protein–ligand docking and screening based on large-scale pre-training. Chem Sci, 2024, 15(4): 1449–1471. doi:10.1039/d3sc05552c""")

# B3 captions: Fig2 content, sources; Fig4 source; Fig5 failure modes
old2 = "**Figure 2. Benchmark design, decoy construction and reporting flow.** Schematic of target curation, decoy generation, docking, enrichment metrics, cluster bootstrap and the three transfer axes; the design separates potency-driven and structure-driven evaluations while keeping a common reporting protocol."
new2 = "**Figure 2. Two-regime data panel composition.** Track A targets (7 primary + 2 auxiliary) are shown with curated active counts, docking-box size and observed decoy ratio; the right panel summarises Track A (1,748 curated actives; 28,136 successfully docked ligands) and Track B (40 inhibitor poses from five of nine curated membrane-target entries). Source: `figures/trackA_enrichment_full.json`, `data/target_panel.csv`, `data/trackA/box_clipped_ligands.csv`; figure file `figures/fig2_benchmark_design.png` / `.pdf`."
assert old2 in t
t = t.replace(old2, new2)

old4 = "**Figure 4. Three-axis transferability.** Temporal (full holdout versus novel scaffold), cross-species (LpxC_Ec versus LpxC_Pa on 951 shared scored decoys) and structure-source (holo versus AFDB, paired n = 34) summaries; the three axes show axis- and target-dependent behaviour rather than a single transferability score."
new4 = "**Figure 4. Three-axis transferability.** Temporal (full holdout versus novel scaffold; Delta AUC with cluster bootstrap CI), cross-species (LpxC_Ec versus LpxC_Pa on 951 shared scored decoys), and structure-source (holo versus AFDB, paired n = 34) summaries; the three axes show axis- and target-dependent behaviour rather than a single transferability score. Source: `results/trackA/axes/temporal.csv`, `cross_species.csv`, `structure_source.csv`; figure file `figures/fig4_three_axis_transfer.png` / `.pdf`."
assert old4 in t
t = t.replace(old4, new4)

old5 = "**Figure 5. Track B cognate redocking and structure-source comparison.** Pose RMSD, PLIF recovery and PoseBusters outcomes for the 40 inhibitor poses, with the exploratory AFDB arm and the cross-crystal-form sensitivity shown separately."
new5 = "**Figure 5. Three benchmark/tool failure modes.** Metal-stratified AUC (four metal primary versus three non-metal primary targets; GyrA and LepB excluded from aggregates), boron rescue sensitivity (with B→C surrogate versus exclude-surrogate for TEM1_Ec and LepB_Ec), and the ECFP4 1-NN versus Vina AUC benchmark-validity comparison. Source: `figures/trackA_enrichment_full.json`, `manuscript/si/table_S5_boron_surrogate.csv`, `results/trackA/ecfp4_baseline/baseline_metrics.csv`; figure file `figures/fig5_trackB_pose_quality.png` / `.pdf`."
assert old5 in t
t = t.replace(old5, new5)

p.write_text(t, encoding="utf-8")
print("M1/references/captions appended")
