from pathlib import Path

p = Path(r"manuscript/v3/manuscript_v3_jcim_style.md")
t = p.read_text(encoding="utf-8")

repls = [
# blockers
("Ligands were prepared at pH 7.4 with RDKit and converted to PDBQT with Meeko using the Vina force field.",
 "Ligands were prepared from canonical SMILES with RDKit (ETKDG embedding, MMFF optimisation, hydrogens added) in neutral form and converted to PDBQT with Meeko using the Vina force field; Meeko failures were recovered with MGLTools prepare_ligand4.py."),
("Holo receptors were prepared with Meeko, preserving cofactors and metals where present.",
 "Holo receptors were prepared with MGLTools 1.5.7 prepare_receptor4.py (with -A bonds_hydrogens and -U nphs_lps), preserving cofactors and metals where present."),
("LpxC Zn stripping changed AUC only from 0.222 to 0.220.",
 "LpxC Zn stripping changed AUC only from 0.220 (Zn-retained) to 0.222 (Zn-stripped)."),
# m8
("Targets with fewer than 50 eligible actives were retained but flagged descriptive in the method axis. Nine targets were retained",
 "LepB, with 3 eligible actives, is descriptive only. Nine targets were retained"),
# M4
("Cluster-aware confidence intervals are wider than ligand-level intervals for scaffold-concentrated targets, and per-target p-values are not used for inference because matched decoys are not independent.",
 "For GyrA, where one scaffold cluster accounts for about 80.8% of actives, the cluster interval was close to the ligand-level interval; resampling was not clustered by chemical series, which is a limitation, and per-target p-values are not used for inference because matched decoys are not independent."),
# M1
("the combination rule was fixed before the calibration gate. Membrane consistency",
 "the combination rule was fixed before the calibration gate. The calibration label includes PoseBusters while the validity sub-score is derived from PoseBusters; this circularity and the nesting of calibration RMSD<2 Å with evaluation RMSD<1 Å are disclosed. S was neutral because phi came from reference structures. Membrane consistency"),
# m3
("The nested design reuses the first ten DiffDock actives per target and their matched decoys, so Boltz-2 is not an independent 806-ligand screen.",
 "The nested design reuses the first ten DiffDock actives per target and their matched decoys, so Boltz-2 is not an independent 806-ligand screen. A 16 GB GPU limited co-folding to 13 of the 40 co-crystal systems (26 observations = 13×2)."),
# m5
("Four of eight receptors had mode-1 redock RMSD above 3 Å (GyrA 13.54 Å, TEM1 11.26 Å, LepB 6.59 Å and GyrB 3.11 Å), and this limitation is carried into the metal and structure-source analyses.",
 "Four of eight receptors had mode-1 redock RMSD above 3 Å (GyrA 13.54 Å, TEM1 11.26 Å, LepB 6.59 Å and GyrB 3.11 Å); best-of-nine redocking recovered 2/8 receptors, and this limitation is carried into the metal and structure-source analyses."),
# m4
("The complete-case sensitivity, stratified AUROC reversals (Vina 0.572 versus Vina-AFDB 0.454 and holo 0.572 versus predicted-receptor 0.454; PBP3 0.107 with rho +0.811) and the 14 non-assessable rows are reported in Table 2 and SI.",
 "The complete-case sensitivity and the 14 non-assessable rows are reported in Table 2 and SI; the stratified AUROC reversals (Vina 0.572 versus Vina-AFDB 0.454 and holo 0.572 versus predicted-receptor 0.454; PBP3 0.107 with rho +0.811) are in SI (S22b)."),
# M5
("The questions were framed before the method-class results were inspected, and the answers are reported as descriptive contrasts rather than a ranking. First,",
 "The four framing questions were defined before the post-hoc method-class arm was added; the arm is exploratory, not a pre-registered test. First,"),
("Each question was answered descriptively; none was used to select a winning method after the fact.",
 "The answers are reported as descriptive contrasts rather than a ranking."),
# M2
("No wet-lab validation was performed, and the co-crystal panel is a cognate redocking set rather than a prospective screening test.",
 "No wet-lab validation was performed, and the co-crystal panel is a cognate redocking set rather than a prospective screening test. The co-folding-subset and membrane-MD analyses were not reported (definitions/status in S24)."),
# M6 + m6
("the all-valid n = 34 version contains 26 QC-excluded pairs and is not used as the headline; using it would change the RMSD/PLIF conclusion.",
 "the all-valid n = 34 version contains 26 QC-excluded pairs and is not used as the headline; using it would change the headline magnitudes (ΔRMSD −2.436 vs −3.493 Å; ΔPLIF −0.519 vs −0.289) but not their direction."),
# m7 duplicate interface
("Only 3 of 40 poses are interface-bound, and receptors are monomeric; Track B therefore tests",
 "Receptors are monomeric; Track B therefore tests"),
("Only 3/40 poses are interface-bound, and AcrB contributes 34 of the 40 main-panel poses, so apparent",
 "AcrB contributes 34 of the 40 main-panel poses, so apparent"),
# M3
("and to retain negative results rather than selecting a single score. No post-hoc threshold",
 "and to retain negative results rather than selecting a single score. An alternative by-stratum aggregation (metal n=5 mean 0.588 vs non-metal n=4 mean 0.498) reverses direction and is therefore not used. No post-hoc threshold"),
# m1/m2/m6 Results additions
("The apparent 2D competitiveness is a property of synthetic decoys filtered to ECFP4 maximum Tanimoto below 0.4 versus all actives, not evidence that docking is unnecessary. Absolute enrichment values are therefore not portable to real screening libraries, and the leaky full-set baseline is retained only as a reproducibility record in SI.",
 "The apparent 2D competitiveness is a property of synthetic decoys filtered to ECFP4 maximum Tanimoto below 0.4 versus all actives, not evidence that docking is unnecessary. Absolute enrichment values are therefore not portable to real screening libraries, and the leaky full-set baseline (median AUC 0.958) is retained only as a reproducibility record in SI. The temporal axis gave median ΔAUC +0.019 (only 2/9 CIs excluded zero); cross-species LpxC_Ec AUC was 0.281 (95% CI 0.265–0.297) and LpxC_Pa 0.810 (0.798–0.821), with paired ΔAUC +0.5288 (corrected decoy-level CI 0.503–0.554). Figure 4 summarises these axes."),
# M7 cross-species duplicate + old isin artefact
("The cross-species comparison of LpxC_Ec and LpxC_Pa on 951 shared decoys gives ΔAUC +0.5288 (corrected decoy-level 95% CI 0.503–0.554), a paired comparison of two active sets rather than a receptor-transfer experiment. The cross-species paired comparison (ΔAUC +0.5288, corrected CI 0.503–0.554) is not a transfer test; the two LpxC homologues differ strongly, and the interval describes the paired active-set contrast.",
 "The cross-species comparison of LpxC_Ec and LpxC_Pa on 951 shared decoys gives ΔAUC +0.5288 (corrected decoy-level 95% CI 0.503–0.554), a paired comparison of two active sets rather than a receptor-transfer experiment. An earlier interval of 0.5135–0.5433 was an isin() artifact and is withdrawn."),
# M7 CPU duplicate
("CPU screening remains inexpensive, but GPU arms add scheduling and reproducibility costs, and pretrained DrugCLIP and Boltz-2 require model-weight and environment disclosure. ", ""),
# m7 title
("Manuscript version: v3 JCIM-style rewrite draft (2026-11-20)",
 "Manuscript version: v3 JCIM-style (2026-11-20)"),
]
for old, new in repls:
    if old not in t:
        print("MISSING:", old[:120])
    else:
        t = t.replace(old, new, 1)

# Replace the duplicated single-chain/no-engine paragraph tail.
old_tail = t[t.index("No engine dominates every target."):t.index("\n\n", t.index("No engine dominates every target."))]
new_tail = ("Per-target ordering varied. DiffDock was higher than Vina for GyrB_Ec (0.682 versus 0.420) and TEM1_Ec (0.598 versus 0.523), but Vina was higher for LpxC_Pa (0.857 versus 0.437) and GyrA_Ec (0.869 versus 0.488). DrugCLIP and Boltz-2 were higher on several targets but lower on GyrA_Ec, and Uni-Dock was nearly tied with Vina overall yet strongly target-dependent. The comparison is therefore target-specific, with no leading method claim. DrugCLIP uses pretrained weights (CC BY-NC 4.0), so train-test overlap cannot be excluded; Uni-Dock positive affinities are treated as not assessable, FabI is excluded from its pooled estimate, and 22 TEM1 boron actives use the pre-specified `_B2C` surrogate. None of these differences supports a leading claim: confidence intervals for the paired differences include zero for four of five methods. The single-chain arm changed DiffDock AUC from 0.488 to 0.667 on GyrA (Δ 0.178; 95% CI 0.021–0.362) and from 0.469 to 0.444 on FabI (Δ −0.025; 95% CI −0.188 to 0.146); it changes only the DiffDock receptor representation, not the Vina ranking. No single engine dominated across targets.")
t = t.replace(old_tail, new_tail, 1)

p.write_text(t, encoding="utf-8")
print("patched", p.stat().st_size, len(t.split()), "words")
