# Benchmarking structure-based virtual screening across Gram-negative antibacterial targets: target-specific enrichment, co-crystal ligand pose limitations and failure modes

**Running title:** Gram-negative SBVS benchmark and failure modes

**Authors:** Yutao Zhao^1, Hongda Liao^2,*

**Affiliations:** ^1 Southern University of Science and Technology, Shenzhen, Guangdong, China; ^2 Shenzhen Bay Laboratory, Shenzhen, Guangdong, China

**Corresponding author:** Hongda Liao, Shenzhen Bay Laboratory, Shenzhen, Guangdong, China; E-mail: hongdaliao@outlook.com


## Abstract

Gram-negative antibacterial discovery faces a dual data problem: many targets have few potent public actives, whereas membrane machines that are structurally tractable seldom provide inhibitor-bound co-crystal poses. Structure-based virtual screening (SBVS) is therefore used under exactly the conditions in which it is least validated. Deep-learning docking and co-folding methods promise better accuracy, but add training-data dependence, GPU cost, receptor-representation sensitivity and benchmark-leakage risk, and balanced head-to-head comparisons remain scarce. We present a two-regime benchmark of nine Gram-negative targets. Track A is a potency-driven panel (1,748 curated actives; 28,136 docked ligands) evaluated by enrichment. Track B is a 40-pose co-crystal benchmark with expert re-annotation (4 inhibitors, 24 antibiotic substrates, 3 covalent β-lactams, 8 fragments and 1 unresolved); a subset-balanced axis compares Vina, Uni-Dock, DiffDock, Boltz-2 and DrugCLIP on 806 ligands, with a 166-ligand Boltz-2 subset and a 206-complex single-chain arm. Track A enrichment is target-specific (AUC 0.300–0.852), and pooled subset AUC was 0.745 for DrugCLIP, 0.718 for Boltz-2, 0.605 for Vina, 0.589 for Uni-Dock (seven targets) and 0.528 for DiffDock; no engine dominated every target. The test-only ECFP4 baseline reached median 0.775, and 7/8 targets matched or exceeded Vina, showing synthetic-decoy inflation. Pre-registered gates failed or were inconclusive: LIT-PCBA did not robustly meet its criteria, and MAES was NO_GO on rank agreement with underpowered AUROC. These results favor per-target, subset-balanced reporting with coverage, decoy ratios, receptor QC and pre-registered gates, not a single ranking. Data and code will be released.

## Keywords

Gram-negative bacteria; structure-based virtual screening; benchmark; co-crystal ligand poses; enrichment; decoy bias; ligand class; failure modes

## Introduction

Antimicrobial resistance is a leading global health threat, and the WHO bacterial priority pathogen list identifies Gram-negative pathogens as an urgent target for new classes [1,2]. The lipopolysaccharide transporter LptB2FGC has recently been targeted by two new antibiotic classes [3,4], renewing interest in Gram-negative targets that are structurally tractable but chemically and biologically difficult. SBVS is the default computational triage for such targets, yet its benchmarks assume abundant potency data and soluble proteins [5-7]. Gram-negative discovery violates both assumptions: public actives are sparse for many targets, and the most attractive membrane machines are often represented by substrate- or fragment-bound structures rather than inhibitor-bound complexes.

Existing benchmarking practice has advanced on several fronts but remains fragmented. DUD-E, DEKOIS and LIT-PCBA established matched-decoy and target-based evaluation [5-9], while recent co-folding and AI-docking studies have expanded the structural side [10-20]. AlphaFold-enabled docking has asked whether predicted receptors can support hit identification [21,22], and pretrained screening models have moved toward ultra-large chemical spaces [23,24]. These studies rarely combine potency-rich and structure-rich regimes under one reporting protocol, and they seldom quantify how receptor preparation, decoy construction, target identity and method class interact. Single-engine benchmarks, absent pre-registration, unobserved decoy ratios and unreported receptor quality can therefore produce rankings that do not transfer.

Toward this goal, we present a benchmark of two complementary Gram-negative regimes with regime-appropriate endpoints. Track A is a nine-target enrichment panel, and Track B is a co-crystal-ligand pose benchmark whose ligand classes are re-annotated before interpretation. We add a subset-balanced method axis spanning Vina, Uni-Dock, DiffDock, Boltz-2 and DrugCLIP, with ECFP4 retained only as a full-panel 2D reference. We pre-register acceptance criteria for an external-consistency anchor and a MAES calibration gate, report negative as well as positive results, and ask what can and cannot be concluded about docking reliability across Gram-negative targets.

## Results

### Data Sets

Track A contains nine Gram-negative antibacterial targets: LpxC from E. coli and P. aeruginosa, DHFR, Gyrase, GyrA, GyrB, FabI, TEM-1 and LepB. Seven targets are primary; GyrA and LepB are auxiliary because of low active counts or low-confidence enrichment. Public activity data were curated from ChEMBL 37 with a pChEMBL threshold of 6, activity-type harmonisation, salt stripping and duplicate removal, yielding 1,748 actives (1,648 primary) and 28,136 successfully docked ligands. Receptors were prepared from holo PDB structures with co-crystal-ligand centres; box sizes were 26 Å for LpxC_Ec, Gyrase, GyrB and LepB and 22 Å for the remaining targets. Gyrase_Ec used the GyrB_Ec receptor, so the metal primary set involves three unique receptors rather than four. Thirty-eight actives (17 primary and 21 auxiliary) fell outside the docking box and were excluded, leaving 1,710 scored actives (1,631 primary) in the analysis table. The temporal split and cross-species comparison use the same panel: 951 shared scored decoys support a paired LpxC_Ec versus LpxC_Pa comparison, and the temporal axis compares a full time-split holdout with a novel-scaffold subset requiring ECFP4 maximum Tanimoto below 0.4 to training. These axes are descriptive and are reported in Figure 4 and Table 3. Figure 1 shows the data landscape, and Figure 2 summarises the two-regime panel composition.

Track B contains 40 co-crystal ligand poses from five of the nine curated membrane-target entries: AcrB 34, PBP3 3, LptB 1, MexB 1 and MsbA 1. BamA and LptD entries contributed no inhibitor-class pose that passed curation, and lipid or substrate analogues were used only for sensitivity analyses (68 clean poses in total). The co-crystal panel is therefore intentionally heterogeneous rather than a matched inhibitor set. Ligand classes were reviewed against RCSB chemical-component names: 4 true pump inhibitors (5QE, 5QF, 5QG and MBX), 24 antibiotic substrates (MIY 10, FUA 7, DM2, DXT, ERY, LFX, RBT, CLM and NOV), 3 covalent β-lactams (AZR, 0W0 and PFV), 8 fragment-like ligands (1K8×2, 1KE, A1AN8, A1AOE, A1AOF, A1AON and XE9) and 1 unresolved/low-confidence ligand (95X; annotated as a substrate in Table S18, not counted as an inhibitor). AcrB has 55 non-polymer poses in total; 34 carry the automated inhibitor label, but only four poses are true pump inhibitors after expert re-annotation. Track B is therefore a co-crystal-ligand pose benchmark with a documented ligand-class re-annotation, not an inhibitor-only benchmark.

The method axis uses a balanced eligible subset from the same panel: an active is eligible when it has a successful Vina pose, was not box-clipped and has at least one matched successful decoy. With seed 42, DiffDock uses 806 rows (403 active plus 403 decoy; eight targets 50+50 and LepB 3+3), Boltz-2 uses a nested 166 rows (eight targets 10+10 and LepB 3+3), and DrugCLIP and Uni-Dock use the same 806 ligands. A single-chain sensitivity arm adds 206 complexes for GyrA, FabI and LepB. The main eight-target pooled analysis covers Vina, DiffDock, Boltz-2 and DrugCLIP; Uni-Dock is pooled over seven targets because FabI has selective positive-affinity missingness, and a common seven-target set is reported for all methods. All subset endpoints are balanced and are not comparable to the full-panel Table 1.

### Track A Benchmark: Per-Target Enrichment

Per-target Vina enrichment spans a wide range (Table 1 and Figure 3). Table 1 reports all nine targets, including the two auxiliary targets. The seven primary targets run from AUC 0.300 (LpxC_Ec) to 0.852 (LpxC_Pa); the metal primary set (n = 4: GyrB, Gyrase, LpxC_Ec and LpxC_Pa) has mean AUC 0.521, median 0.466 and mean EF5% 1.780, whereas the non-metal set (n = 3: DHFR, FabI and TEM1) has mean AUC 0.567, median 0.606 and mean EF5% 1.028. EF1% is reported only when at least 150 actives are available. GyrA (AUC 0.856) and LepB (AUC 0.292) are auxiliary and are excluded from both aggregates and from any trend test. Across all nine targets the AUC range is 0.292–0.856 and EF5% ranges from 0.00 (GyrB) to 7.47 (GyrA); these values are descriptive because the auxiliary targets contain only 6–73 actives. The metal aggregate involves three receptors because Gyrase and GyrB share a receptor, so it is descriptive rather than a receptor-replication test.

**Table 1. Track A per-target enrichment (Vina, primary decoy set).** AUC is shown with active-cluster bootstrap 95% CI; EF1% is reported only when n_active >= 150; n.r. = not reported. GyrA and LepB are auxiliary.
| Target | Role | n_active | AUC (95% CI) | EF5% | EF1% |
|---|---|---:|---|---:|---:|
| LpxC_Ec | primary | 222 | 0.300 (0.268–0.334) | 0.09 | 0.00 |
| LpxC_Pa | primary | 418 | 0.852 (0.836–0.869) | 5.64 | 4.57 |
| DHFR_Ec | primary | 383 | 0.666 (0.642–0.689) | 1.83 | 2.08 |
| Gyrase_Ec | primary | 389 | 0.578 (0.548–0.608) | 1.39 | 0.77 |
| GyrB_Ec | primary | 75 | 0.354 (0.298–0.407) | 0.00 | n.r. |
| FabI_Ec | primary | 53 | 0.606 (0.519–0.688) | 0.37 | n.r. |
| TEM1_Ec | primary | 91 | 0.430 (0.370–0.488) | 0.88 | n.r. |
| GyrA_Ec | auxiliary | 73 | 0.856 (0.811–0.902) | 7.47 | n.r. |
| LepB_Ec | auxiliary | 6 | 0.292 (0.067–0.575) | 3.67 | n.r. |

The LpxC paralogues illustrate the target specificity most clearly. LpxC_Ec is inverted (AUC 0.300, 95% CI 0.268–0.334), whereas LpxC_Pa is strongly enriched (AUC 0.852, 95% CI 0.836–0.869); the difference survives property matching and nearest-decoy analysis and is not explained by exclusion of the largest P. aeruginosa scaffolds. GyrA is also enriched (0.856, 95% CI 0.811–0.902), while GyrB (0.354, 95% CI 0.298–0.407) and LepB (0.292, 95% CI 0.067–0.575) are not. These contrasts argue against a single Gram-negative enrichment factor and against interpreting a global metal effect.

Sensitivity analyses bound the main conclusions. A three-tier LpxC_Ec box experiment gave AUC values of 0.224, 0.289 and 0.333 at 22, 26 and 28 Å; ranking was stable (Spearman rho 0.949 for 26 versus 28 Å), but the absolute box effect (0.109) exceeds the metal/non-metal mean difference (0.046), so the 26 Å box is retained as primary without attributing a causal metal effect. Exhaustiveness 4 versus 8 and decoy ratios of 10, 15 and 20:1 produced small absolute changes (|ΔAUC| ≤ 0.0116; Tables S4b and S8). For GyrA, where one scaffold cluster accounts for about 80.8% of actives, the cluster interval was close to the ligand-level interval; resampling was not clustered by chemical series, which is a limitation, and per-target p-values are not used for inference because matched decoys are not independent. Vinardo rescoring changed absolute AUC (LpxC_Ec 0.300 to 0.554, LpxC_Pa 0.852 to 0.918 and Gyrase 0.578 to 0.518) while target rankings remained moderately concordant (Spearman rho 0.65–0.84), so the single Vina engine is a sensitivity choice rather than a claim of scoring optimality.

The test-only ECFP4 1-nearest-neighbour baseline provides the benchmark-validity check. On temporal test actives alone, it reaches median AUC 0.7745 (range 0.4992–1.0000); five of nine targets are below 0.80, yet seven of eight comparable targets match or exceed Vina. The apparent 2D competitiveness is a property of synthetic decoys filtered to ECFP4 maximum Tanimoto below 0.4 versus all actives, not evidence that docking is unnecessary. Absolute enrichment values are therefore not portable to real screening libraries, and the leaky full-set baseline (median AUC 0.958) is retained only as a reproducibility record in SI. The temporal axis gave median ΔAUC +0.019 (only 2/9 CIs excluded zero); cross-species LpxC_Ec AUC was 0.281 (95% CI 0.265–0.297) and LpxC_Pa 0.810 (0.798–0.821), with paired ΔAUC +0.5288 (corrected decoy-level CI 0.503–0.554). Figure 4 summarises these axes.

### Track B Benchmark: Co-Crystal Ligand Poses

Cognate redocking recovered 9 of 36 poses below 2 Å, with median ligand RMSD 4.44 Å and median PLIF recovery 0.896. Strict PoseBusters physics and chemistry checks gave 32 passed, 4 failed and 4 not assessable because no reference coordinates were available. Receptors are monomeric; Track B therefore tests pose and interface recovery under cognate conditions, not cross-docking or membrane-embedded assembly. Only 3 of the 40 poses are interface-bound.

The ligand-class composition is itself a benchmark-validity result. Of the 40 poses, 4 are true pump inhibitors, 24 are antibiotic substrates, 3 are covalent β-lactams, 8 are fragment-like and 1 is unresolved/low-confidence. AcrB contributes 55 non-polymer poses in total; 34 carry the automated inhibitor label, but expert re-annotation identifies only four true inhibitors. AcrB contributes 34 of the 40 main-panel poses, so apparent ligand-class effects are also target-composition effects; the class medians are descriptive and not tested. Class-stratified descriptive medians (Table S18b; small subgroups, no inference) were: true inhibitors n = 4, median RMSD 0.95 Å and PLIF 1.00; antibiotic substrates n = 24, 3.84 Å and 0.79; covalent β-lactams n = 3, 5.67 Å and 0.86; fragments n = 8 (4 with valid RMSD), 4.66 Å and 1.00; unresolved/substrate-like n = 1, 5.43 Å and 1.00.

The predicted-receptor arm passed QC for 12 of 40 systems and excluded 28 of 40 (21 ligand clashes and 7 TM-score below 0.70). Restricting to the eight QC-passed pairs, holo-to-AFDB mean RMSD changed from 5.24 to 1.74 Å (Δ −3.493 Å, 95% CI −4.711 to −2.064) and PLIF recovery dropped from 0.818 to 0.528 (Δ −0.289, 95% CI −0.418 to −0.180); PoseBusters passed 7/1 for holo and 7/1 for AFDB. The all-valid n = 34 version (mean RMSD 3.87→1.43 Å; ΔPLIF −0.519) contains 26 QC-excluded pairs and is SI-only. The predicted-receptor hypothesis was negative/underpowered: the original acceptable-pose-rate endpoint difference was 0, the endpoint was replaced after the pre-specified QC trigger, and predicted receptors do not support a confirmatory transfer claim.

**Table 3. Exploratory axis summary.**
| Axis | Key metric | Interpretation |
|---|---|---|
| Temporal (novel scaffold) | median ΔAUC +0.019; Gyrase_Ec −0.147; DHFR_Ec +0.084; 2/9 CIs exclude 0 | descriptive/mixed |
| Cross-species (LpxC Ec vs Pa) | ΔAUC +0.5288; corrected decoy-level CI 0.503–0.554 | paired comparison, not transfer |
| Structure source (QC-passed n=8) | ΔRMSD −3.493 (CI −4.711 to −2.064); ΔPLIF −0.289 (CI −0.418 to −0.180) | exploratory; predicted-receptor hypothesis negative |

### Method-Class Comparison

The method-class comparison is post-hoc and subset-balanced. It compares Vina and Vinardo as a second scoring function, Uni-Dock v1.2.0 [25] under Windows Subsystem for Linux with the same receptor, box and parameters as Vina, DiffDock v1.1 with ESM2-650M, Boltz-2 with ligand_iptm as the primary score, and DrugCLIP with the official 6_folds model. ECFP4 is not a subset-scored method: its test-only reference was not recomputed because 68–94% of subset actives overlap the temporal-train reference for seven of eight comparison targets (TEM1 38%, LepB 100% at n=3; S25), so the 2D reference remains the full-panel result in §3.3. Endpoints are per-target AUC, top-5 active recall and Spearman rho; LepB is descriptive only. Per-target exceptions were substantial: DiffDock exceeded Vina on GyrB and TEM1 but fell below Vina on LpxC_Pa and GyrA, whereas Boltz-2 and DrugCLIP were lower on GyrA than on several other targets. The main pooled set covers Vina, DiffDock, Boltz-2 and DrugCLIP across eight targets; Uni-Dock is pooled over seven targets excluding FabI, and a common seven-target set is reported for all methods.

Figure 7 summarises the method-class comparison. Pooled AUC on the balanced subset was 0.745 for DrugCLIP (95% CI 0.639–0.852), 0.718 for Boltz-2 (0.585–0.840), 0.605 for Vina, 0.589 for Uni-Dock on seven targets (0.425–0.745) and 0.528 for DiffDock (0.479–0.583). Pairwise differences versus Vina were −0.077 for DiffDock (95% CI −0.233 to 0.081), 0.089 for Boltz-2 (−0.120 to 0.303), 0.139 for DrugCLIP (−0.001 to 0.291) and −0.002 for Uni-Dock (−0.013 to 0.008). Coverage was complete for Vina, 13/806 not assessable for DiffDock, 1/166 for Boltz-2, 5/806 for DrugCLIP (all GyrA decoys with RDKit 3D generation failure) and 31/806 for Uni-Dock (30 positive-affinity and one no-model result; FabI 86/14). The 7-target common sensitivity gave means of 0.605 for Vina, 0.536 for DiffDock, 0.735 for Boltz-2, 0.768 for DrugCLIP and 0.589 for Uni-Dock.

Per-target ordering varied. DiffDock was higher than Vina for GyrB_Ec (0.682 versus 0.420) and TEM1_Ec (0.598 versus 0.523), but Vina was higher for LpxC_Pa (0.857 versus 0.437) and GyrA_Ec (0.869 versus 0.488). DrugCLIP and Boltz-2 were higher on several targets but lower on GyrA_Ec, and Uni-Dock was nearly tied with Vina overall yet strongly target-dependent. The comparison is therefore target-specific, with no leading method claim. DrugCLIP uses pretrained weights (CC BY-NC 4.0), so train-test overlap cannot be excluded; Uni-Dock positive affinities are treated as not assessable, FabI is excluded from its pooled estimate, and 22 TEM1 boron actives use the pre-specified `_B2C` surrogate. None of these differences supports a leading claim: confidence intervals for the paired differences include zero for four of five methods. The single-chain arm changed DiffDock AUC from 0.488 to 0.667 on GyrA (Δ 0.178; 95% CI 0.021–0.362) and from 0.469 to 0.444 on FabI (Δ −0.025; 95% CI −0.188 to 0.146); it changes only the DiffDock receptor representation, not the Vina ranking. No single engine dominated across targets.

### Pre-Registered Gates and Failure Modes

The LIT-PCBA anchor was run on all 15 sets. The anchor's acceptance criteria AC1–AC7 are listed in Table 4. AC1 passed (15/15), AC2 was 2/9 = 22.2% of the nine evaluable sets (GBA and KAT2A), below the 30% threshold; including TP53 gives 3/10 = 30.0% (Wilson 95% CI 6.3–54.7). AC3 failed under both SD definitions (maximum CV 28.28%, reported as 28.3% under population SD, and 34.64% under sample SD, MTORC1), while AC4, AC6 and AC7 passed; AC5 was not applicable because the anchor was CPU-only. The 15-set descriptive view gave EF5% ≥ 2 for 5/15 (33.3%), median AUC 0.547 and median EF5% 1.19; EF1% was reported for ALDH1 (1.21) and VDR (1.28) only, and the anchor remains an appendix consistency check, not a validation claim.

**Table 4. LIT-PCBA anchor ACs (appendix).**
| AC | Criterion | Primary result | Sensitivity/notes |
|---|---|---|---|
| AC1 | 15/15 attempted; >=10/15 completed | 15/15 attempted and completed | PASS |
| AC2 | >=30% of evaluable sets with EF5% >= 2 | 2/9 = 22.2% | including TP53: 3/10 = 30.0% (Wilson 6.3–54.7); not robust |
| AC3 | EF5% CV across three inactive subsamples <=10% | 28.3% population SD | 34.64% sample SD; FAIL either way |
| AC4 | all 15 rows retained | 15/15 rows | PASS |
| AC5 | AI arm if run | not applicable (CPU-only) | n/a |
| AC6 | freeze precedes results | freeze before results | PASS |
| AC7 | one-command reproduction | demo reproducible | PASS |

**Figure 6** displays the calibration-gate result. The MAES calibration gate used 146 pose observations across four arms: Vina (holo receptors) 40, Vinardo 40, Vina with predicted receptors 40 and Boltz-2 26. The primary validity-plus-membrane score gave AUROC 0.496 (wild-cluster CI 0.318–0.672), signed rho −0.0745 (|rho| 0.075; cluster CI −0.302 to 0.161), Brier score 0.145 [26] and expected calibration error 0.124 [27]. The rank criterion is decisive (upper CI 0.161 < 0.40), whereas the AUROC criterion is inconclusive and underpowered (19 positives). Excluding 14 non-assessable rows leaves n = 132, AUROC 0.434 and Brier score 0.159 (Brier skill score, BSS, −0.281 full-sample and −0.290 complete-case; negative values indicate worse than the no-skill baseline). The complete-case sensitivity and the 14 non-assessable rows are reported in Table 2 and SI; the stratified AUROC reversals (Vina 0.572 versus Vina-AFDB 0.454 and holo 0.572 versus predicted-receptor 0.454; PBP3 0.107 with rho +0.811) are in SI (S22b). MAES is retained only as a descriptive pose-quality score, not as a calibrated probability.

**Table 2. MAES calibration gate (frozen plan; primary V+M).** Manifest sha256 e4ae5c23461474f773a016ec415e9da8e5dd6ff10bbef55eda73c726cb4fa3ca. CIs are wild cluster bootstrap (AUROC) and cluster bootstrap (rho), with 40 clusters.
| Analysis | n | Evaluation positives | AUROC (95% CI) | rho (absolute; 95% CI) | Brier | ECE |
|---|---:|---:|---|---|---:|---:|
| Primary V+M | 146 | 19 | 0.496 (0.318–0.672) | 0.075 (−0.302 to 0.161) | 0.145 | 0.124 |
| Complete-case V+M | 132 | 19 | 0.434 (complete-case sensitivity; no CI reported) | — | 0.159 | — |
| Exploratory V+M+S | 146 | 19 | 0.495 (0.318–0.671) | 0.072 (−0.301 to 0.166) | 0.145 | 0.124 |

**Figure 5** shows the failure-mode comparison. These three modes are reported together because they change conclusions in different ways. Three benchmark/tool failure modes recur. Metal coordination is target-specific rather than a stable penalty [28-31]: the metal primary set overlaps the non-metal set, Gyrase and GyrB share a receptor, receptor-preparation quality is uneven (4/8 top-ranked mode-1 redock RMSD >3 Å) and box-size effects exceed the metal difference; LpxC Zn stripping changed AUC only from 0.220 (Zn-retained) to 0.222 (Zn-stripped). Second, Vina 1.2.5 lacks boron atom typing, so boron-containing actives required a pre-specified B-to-C isosteric surrogate. Third, synthetic decoys filtered to ECFP4 maximum Tanimoto below 0.4 make the test-only 2D baseline competitive with Vina on seven of eight targets. The benchmark response is to report per-target estimates, observed decoy ratios, coverage, receptor QC and pre-registered gates, and to retain negative results rather than selecting a single score. An alternative by-stratum aggregation (metal n=5 mean 0.588 vs non-metal n=4 mean 0.498) reverses direction and is therefore not used. No post-hoc threshold, method or observation was changed after these results were observed.

## Discussion

This work addresses a gap between how SBVS is used in Gram-negative antibacterial discovery and how it is benchmarked. It separates benchmark-construction artefacts from target biology because the same panel supports both a full enrichment table and a balanced paired subset. Potency-rich panels and co-crystal-ligand panels answer different questions, and neither alone describes whether a docking score can be trusted for a new Gram-negative target. By placing both regimes under one reporting protocol and adding a balanced method axis, the study separates target-specific enrichment from benchmark-construction artefacts and from method-class differences.

The four framing questions were defined before the post-hoc method-class arm was added; the arm is exploratory, not a pre-registered test. First, does a potency-driven panel give stable per-target enrichment, or do a few targets dominate the aggregate? Second, does a co-crystal-ligand benchmark contain the inhibitor classes its name implies? Third, do modern docking and co-folding methods change the reliability ordering on a balanced subset? Fourth, do pre-registered gates and failure-mode audits expose limitations that a single headline number would hide? The answers are reported as descriptive contrasts rather than a ranking.

Several caveats follow from the panel and data choices. The unbalanced target sizes also mean that pooled means weight each target equally, not each active. The nine-target panel is broad for Gram-negative SBVS but still small, with scored actives per target ranging from 6 to 418 after box-clipping. The balanced method subset deliberately sacrifices panel-wide coverage for matched active-decoy pairs, and its 1:1 composition is not comparable to the full-panel enrichment values; pooled subset values such as Vina 0.605 are interpreted only as method-class contrasts. The LIT-PCBA anchor was CPU-only and used a capped inactive ratio, with ALDH1 at 0.74:1; AC2 and AC3 are sensitive to reasonable analysis choices. No wet-lab validation was performed, and the co-crystal panel is a cognate redocking set rather than a prospective screening test. The co-folding-subset and membrane-MD analyses were not reported (definitions/status in S24).

**Figure 4** summarises the exploratory axes. The results do not support a single reliability ranking. Target-specific enrichment spans AUC 0.300–0.852, and the LpxC paralogues move in opposite directions. Temporal variation is mixed: only two of nine targets have confidence intervals excluding zero for the novel-scaffold comparison, so the temporal axis is descriptive, not a transfer test. The cross-species comparison of LpxC_Ec and LpxC_Pa on 951 shared decoys gives ΔAUC +0.5288 (corrected decoy-level 95% CI 0.503–0.554), a paired comparison of two active sets rather than a receptor-transfer experiment. An earlier interval of 0.5135–0.5433 was an isin() artifact and is withdrawn. The practical implication is that a target-specific pilot and a test-only 2D baseline remain necessary before a docking ranking is trusted.

Target-specific bias and confounding are too large to collapse into one correction. Metal coordination differs between LpxC homologues despite similar zinc sites, and the metal aggregate is confounded by Gyrase/GyrB receptor sharing, receptor-preparation quality and box size. The observed decoy ratio was 9.82–19.98 rather than the planned 20:1, and decoys were filtered to ECFP4 maximum Tanimoto below 0.4 versus actives; both choices affect absolute enrichment. Cluster-aware intervals and per-target reporting reduce overconfidence but do not remove these biases, and the metal/non-metal summary is therefore reported as descriptive only.

Structure-source dependence is another axis that resists a simple conclusion. The QC-passed holo-versus-AFDB comparison (n = 8) shows a large RMSD improvement alongside reduced interface recovery, while the all-valid n = 34 version contains 26 QC-excluded pairs and is not used as the headline; using it would change the headline magnitudes (ΔRMSD −2.436 vs −3.493 Å; ΔPLIF −0.519 vs −0.289) but not their direction. The predicted-receptor hypothesis is negative/underpowered and its endpoint was replaced after the pre-specified QC trigger. In the co-crystal panel, ligand classes are dominated by antibiotic substrates and fragments, and AcrB contributes 55 non-polymer poses but only 34 automated inhibitor labels, of which only four are true inhibitors. A method that appears strong on the pose set may therefore be recovering substrate-like or fragment-like geometry rather than inhibitor binding. The single-chain DiffDock sensitivity changes receptor representation only and leaves the Vina ordering unchanged.

Practical utility is bounded by cost and data access. The practical cost is asymmetric: CPU-only Vina remains inexpensive, whereas DiffDock and Boltz-2 require GPU time (and Boltz-2's nested subset is deliberately small), and GPU arms require scheduling, model-weight provenance and deterministic environment disclosure. DrugCLIP adds a pretrained model whose training overlap with public actives cannot be excluded; Uni-Dock adds GPU docking with positive-affinity failures and a boron surrogate requirement. These costs do not disqualify the methods, but they change the appropriate claim: a method can be useful for triage without being a portable enrichment predictor. All curated data, decoy manifests, docking scores, analysis scripts and figure-generation code will be released upon publication, and the negative gates remain part of the resource.

We suggest four checks before a Gram-negative docking ranking is trusted. Run a per-target pilot rather than relying on a single aggregate; include a test-only 2D baseline and report the observed decoy ratio; verify receptor preparation, box size and metal handling for the chemistry at hand; and report per-target estimates with cluster-aware intervals plus pre-registered decision gates. These checks are themselves limited: they require enough actives for per-target estimation, they do not replace prospective validation, and they cannot remove pretrained-model overlap. They are intended as an audit protocol, not a guarantee of prospective success; per-target pilots and confidence intervals require enough actives, and matched-decoy controls can themselves introduce bias.

Future work should extend the panel to additional Gram-negative targets, prospective screening libraries and membrane-embedded receptor states; test whether ligand-class-aware co-crystal sets improve pose-benchmark validity; and evaluate calibration and ranking under distribution shift rather than only on matched subsets. The released panel, decoy manifests, scores and analysis code are intended as an audit resource for such work and for groups comparing docking or co-folding tools before committing to a synthesis campaign. A prospective test across several target classes with pre-registered enrichment and pose endpoints would be the appropriate next step. The benchmark therefore functions as a decision aid, not as a replacement for target-specific experimental validation, and its value will grow as more prospective target panels and negative results are deposited.

## Methods

### Data curation and panel

Track A targets were selected from ChEMBL 37 using a pChEMBL threshold of 6, followed by activity-type harmonisation, salt stripping and duplicate removal. LepB, with 3 eligible actives, is descriptive only. Nine targets were retained (seven primary and two auxiliary); scored actives per target range from 6 to 418 after exclusion. GyrA and LepB were designated auxiliary because of low active counts or low-confidence enrichment. Public actives totalled 1,748 (1,648 primary), and 28,136 ligands were successfully docked. Receptors were prepared from holo PDB entries with co-crystal-ligand centres; box sizes were 26 Å for LpxC_Ec, Gyrase, GyrB and LepB and 22 Å for the remaining targets. Auxiliary targets were retained for descriptive contrast but not for metal/non-metal aggregates, and a single box centre was used per receptor.

### Ligand preparation and matched decoys

Ligands were prepared from canonical SMILES with RDKit (ETKDG embedding, MMFF optimisation, hydrogens added) in neutral form and converted to PDBQT with Meeko using the Vina force field; Meeko failures were recovered with MGLTools prepare_ligand4.py. Decoys followed the DUD-E-style property-matched protocol with three independent sets per target; the primary analysis used set s1, and sets s2 and s3 were used for stability checks. The protocol targeted up to 20 matched decoys per active, but the observed inactive:active ratio was 9.82–19.98 (median 15.15), and the observed ratio is reported per target. Candidate decoys were filtered to ECFP4 maximum Tanimoto below 0.4 versus all target actives before property matching. Ligands with no successful Vina pose or with box clipping were excluded from the method-axis subset, and all decoy sets were constructed before method scores were inspected.

### Receptor preparation and QC

Holo receptors were prepared with MGLTools 1.5.7 prepare_receptor4.py (with -A bonds_hydrogens and -U nphs_lps), preserving cofactors and metals where present. A receptor manifest recorded grid centres and box sizes; redocking QC used the top-ranked mode-1 pose and a best-of-nine check. Four of eight receptors had mode-1 redock RMSD above 3 Å (GyrA 13.54 Å, TEM1 11.26 Å, LepB 6.59 Å and GyrB 3.11 Å); best-of-nine redocking recovered 2/8 receptors, and this limitation is carried into the metal and structure-source analyses. Predicted receptors used AlphaFold DB models aligned to the holo pocket; QC excluded 28 of 40 systems (21 ligand clashes and 7 TM-score below 0.70), leaving eight QC-passed pairs for the headline comparison. The same receptor manifest was used across Vina, Uni-Dock and DiffDock main to avoid engine-specific receptor preparation.

### Vina and enrichment metrics

AutoDock Vina 1.2.5 [32] was run with a single centre and box per receptor, seed 42, exhaustiveness 4 for EF1-enabled targets (LpxC_Ec, LpxC_Pa, DHFR_Ec and Gyrase_Ec) and exhaustiveness 2 for non-EF1 targets; sensitivity arms used exhaustiveness 8. The primary endpoint was EF5% = (hits in top k / n_active) / (k / N), with k = max(1, round(0.05 N)); EF1% was reported when n_active ≥ 150, and ROC-AUC and BEDROC (α = 20) were secondary metrics. Uncertainty used a 2,000-resample active-level cluster bootstrap, with decoys grouped by matched active; per-target p-values are not used for inference, and Benjamini–Hochberg values [33] are not used; reporting follows DOME recommendations [34].

### Uni-Dock

Uni-Dock v1.2.0 was run under Windows Subsystem for Linux with the same receptor PDBQT, box and parameters as Vina (`--scoring vina --exhaustiveness 4 --seed 42 --num_modes 3`). Positive affinities were treated as not assessable, one FabI ligand produced no model, and 22 TEM1 boron-containing actives used the pre-specified `_B2C` surrogate; FabI is excluded from the pooled Uni-Dock estimate because of selective positive-affinity missingness. The common seven-target sensitivity set is reported for all methods. The same 806 ligands were run without result filtering, and all failures are reported in S30.

### DiffDock

DiffDock v1.1 [35] used ESM2-650M, two samples per complex and rank-1 selection by confidence, with the same 806-ligand subset and seed 42. The main arm used the multichain holo receptor; a single-chain sensitivity arm used the box-centre-nearest CA chain for GyrA, FabI and LepB (206 complexes). DiffDock input retained HETATM cofactors, but the DiffDock/ProDy parser models protein CA only; cofactors are therefore not scored. Thirteen complexes failed and are treated as not assessable. No pose post-processing or rescoring was applied; confidence rank 1 was used exactly as produced.

### Boltz-2

Boltz-2 [11] used ligand_iptm as the primary score and confidence_score as a secondary score, with empty MSA input (single-sequence mode; no MSA server) and no metal or cofactor ion entity. The nested subset contains 166 ligands (eight targets 10+10 and LepB 3+3), run with default parameters and seed 42; 165 succeeded and one GyrA decoy failed after two attempts. Top-5 recall and EF were not reported for this 20-ligand-per-target design. The nested design reuses the first ten DiffDock actives per target and their matched decoys, so Boltz-2 is not an independent 806-ligand screen. A 16 GB GPU limited co-folding to 13 of the 40 co-crystal systems (26 observations = 13×2).

### DrugCLIP

DrugCLIP [36] used the official 6_folds model with native-ligand 6 Å pocket encoding, which excludes metals and cofactors; the same 806 ligands were screened with seed 42. The local deployment passed an 11/11 rank-consistency smoke test against the official reference before screening; five GyrA decoys failed RDKit 3D conformer generation and are reported as not assessable. The model and outputs are distributed under CC BY-NC 4.0, and train-test overlap with public actives cannot be excluded. Scores were computed without target-specific fine-tuning, and the five failed decoys remain in coverage rather than being replaced.

### Pose-quality metrics

Track B used symmetry-aware ligand RMSD, PLIF recovery, PoseBusters physics/chemistry checks and interface-contact counts [16]. MAES combines a validity term from PoseBusters with a membrane-consistency term from the OPM reference plane [37] and a state-consistency term based on TM-helix angle, pocket volume and gating distance; the combination rule was fixed before the calibration gate. The calibration label includes PoseBusters while the validity sub-score is derived from PoseBusters; this circularity and the nesting of calibration RMSD<2 Å with evaluation RMSD<1 Å are disclosed. S was neutral because phi came from reference structures. Membrane consistency was assessed from reference-relative geometry without molecular-dynamics simulations [38,39]. PoseBusters was used only for physics/chemistry validity; RMSD-like checks were excluded to avoid circularity.

### Statistics and bootstrap

All confidence intervals are percentile bootstrap unless stated otherwise. Main Track A cluster intervals use 2,000 resamples; the MAES calibration-gate wild-cluster AUROC uses 4,999 Rademacher draws and its rho uses 2,000 cluster bootstrap resamples; corrected cross-species intervals use 10,000 draws. Pooled method-axis AUCs are unweighted target-level means with target-level bootstrap (2,000 draws); per-target intervals resample ligands within target, and paired comparisons resample matched subset pairs rather than comparing independent confidence intervals. No multiple-comparison correction or significance claim is made; per-target estimates are the primary descriptive unit. We report minimum detectable effect sizes (MDES) rather than post-hoc power: Track A detects approximately 0.93–0.99 × between-target SD; Track B effective n_eff is 20–35 and detects only ≥15–25 percentage-point acceptable-rate differences, or MAES differences of about 0.10–0.16.

### Pre-registered gates and amendments

Two decision gates were fixed before their analyses: the LIT-PCBA external-consistency anchor (AC1–AC7; Table 4) and the MAES calibration gate (AUROC ≥ 0.65 and |rho| ≥ 0.40). The anchor used a 20:1 inactive:active target, a cap of 1,000 inactives per set, three inactive subsample seeds and all 15 sets; ALDH1 reached only 0.74:1 because of the cap. The MAES gate used 146 pose observations and five-fold grouped cross-validation by system. The method-axis addenda and amendment hashes are listed in S24; no threshold, method or sample was changed after observing the gate results. The amendments and addenda document the pre-specified analysis plan and all deviations.

### Data and software release

The released data, decoy manifests, docking scores, analysis scripts, figure-generation code and intermediate tables are described in the back-matter Data and Software Availability statement; SI tables S1–S31 accompany the submission.

## Author contributions

Yutao Zhao: Conceptualization, Methodology, Software, Investigation, Data curation, Visualization, Writing – original draft. Hongda Liao: Conceptualization, Supervision, Resources, Writing – review & editing.

## Acknowledgments

Funding: none. During the preparation of this work the authors used DeepSeek V4.1 for language editing and for reference checking, used AI assistance for data-processing script preparation, and prepared the figures with Python matplotlib. These tools were used to improve the clarity of the text, to verify that the cited references are reported accurately, and to generate the plots. The authors reviewed and verified all content and take full responsibility for the integrity of the manuscript, and no content was generated by AI without author verification.

## Competing interests

The authors declare no competing financial interest.

## Data and Software Availability

**Datasets.** All curated data, decoy manifests, docking scores and intermediate tables will be released upon publication. The MAES calibration-gate manifest and report are identified by sha256 e4ae5c23461474f773a016ec415e9da8e5dd6ff10bbef55eda73c726cb4fa3ca. Intermediate tables (including per-target enrichment, sensitivity analyses, anchor ACs, MAES calibration-gate metrics, ligand-class review, QC sensitivity and corrected cross-species intervals) are provided as SI tables S1–S31 (including sub-tables S4a/S4b, S13b, S16b, S18b, S19b, S22b, S22c, S26b, S27b, S31b, S31c, S31d and S31e).

**Software.** Analysis scripts, figure-generation code, environment manifests and the raw score/coverage files used for every pooled estimate will be released upon publication.

## Supporting Information

31 supporting tables (S1–S31) covering per-target enrichment, sensitivity analyses, method-class metrics, coverage/NA and calibration details (XLSX: SI_tables_S1-S31.xlsx); the workbook includes sub-tables S4a/S4b, S13b, S16b, S18b, S19b, S22b/S22c, S26b, S27b and S31b–S31e.

## References

[1] Antimicrobial Resistance Collaborators. Global burden of bacterial antimicrobial resistance in 2019: a systematic analysis. Lancet, 2022, 399: 629–655. doi:10.1016/S0140-6736(21)02724-0

[2] World Health Organization. WHO Bacterial Priority Pathogens List, 2024. Geneva: WHO, 2024.

[3] Zampaloni C, Mattei P, Bleicher K, et al. A novel antibiotic class targeting the lipopolysaccharide transporter. Nature, 2024, 625: 566–571. doi:10.1038/s41586-023-06873-0

[4] Pahil K S, Gilman M S A, Baidin V, et al. A new antibiotic traps lipopolysaccharide in its intermembrane transporter. Nature, 2024, 625: 572–577. doi:10.1038/s41586-023-06799-7

[5] Mysinger M M, Carchia M, Irwin J J, Shoichet B K. Directory of Useful Decoys, Enhanced (DUD-E): better ligands and decoys for better benchmarking. J Med Chem, 2012, 55(14): 6582–6594. doi:10.1021/jm300687e

[6] Vogel S M, Bauer M R, Boeckler F M. DEKOIS: demanding evaluation kits for objective in silico screening — a versatile tool for benchmarking docking programs and scoring functions. J Chem Inf Model, 2011, 51(10): 2650–2665. doi:10.1021/ci2001549

[7] Tran-Nguyen V K, Jacquemard C, Rognan D. LIT-PCBA: an unbiased data set for machine learning and virtual screening. J Chem Inf Model, 2020, 60(9): 4263–4273. doi:10.1021/acs.jcim.0c00155

[8] Sieg J, Flachsenberg F, Rarey M. In need of bias control: evaluating chemical data for machine learning in structure-based virtual screening. J Chem Inf Model, 2019, 59(3): 947–961. doi:10.1021/acs.jcim.8b00712

[9] Bauer M R, Ibrahim T M, Vogel S M, Boeckler F M. Evaluation and optimization of virtual screening workflows with DEKOIS 2.0 — a public library of challenging docking benchmark sets. J Chem Inf Model, 2013, 53(6): 1447–1462. doi:10.1021/ci400115b

[10] Abramson J, Adler J, Dunger J, et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. Nature, 2024, 630: 493–500. doi:10.1038/s41586-024-07487-w

[11] Passaro S, Corso G, Wohlwend J, et al. Boltz-2: towards accurate and efficient binding affinity prediction. bioRxiv, 2025. doi:10.1101/2025.06.14.659707

[12] Obendorf L, Doering N P, Knaus P, Wolber G. Benchmarking protein–ligand cofolding models: correct ligand placement is decoupled from accurate protein conformation. J Chem Inf Model, 2026, 66(21): 14188–14199. doi:10.1021/acs.jcim.6c02169

[13] Dobson L, Tusnády G E, Tompa P. Regularly updated benchmark sets for statistically correct evaluations of AlphaFold applications. Brief Bioinform, 2025, 26(2): bbaf104. doi:10.1093/bib/bbaf104

[14] Truong C M, To V T, Janel N, et al. Enhancing virtual screening of cystathionine beta-synthase inhibitors: benchmarking target-specific machine-learning scoring functions against state of the art AI docking and co-folding approaches. J Cheminform, 2026, 18(1): 141. doi:10.1186/s13321-026-01262-x

[15] Ahmed F, Soellner M B, Brooks C L III. Real-world assessment of machine-learned docking using bioassay-derived benchmarks. J Chem Inf Model, 2026, 66(15): 8752–8759. doi:10.1021/acs.jcim.5c03020

[16] Buttenschoen M, Morris G M, Deane C M. PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences. Chem Sci, 2024, 15(9): 3130–3139. doi:10.1039/d3sc04185a

[17] Yang L Y, Gupta S, Mullininx L N, McShan A C. The LiPP benchmark set for modeling lipid–protein complexes: comparison of co-folding and docking methods. J Chem Inf Model, 2026, 66(12): 7231–7251. doi:10.1021/acs.jcim.6c01457

[18] Oheda H, Inoue M, Ekimoto T, Yamane T, Ikeguchi M. Co-folding of membrane proteins and lipid molecules improves membrane-protein structure prediction accuracy. openRxiv preprint, posted 2026-03-25. doi:10.64898/2026.03.24.713847

[19] Lu W, Zhang J, Huang W, et al. DynamicBind: predicting ligand-specific protein–ligand complex structure with a deep equivariant generative model. Nat Commun, 2024, 15(1): 1071. doi:10.1038/s41467-024-45461-2

[20] Cai H, Shen C, Jian T, et al. CarsiDock: a deep learning paradigm for accurate protein–ligand docking and screening based on large-scale pre-training. Chem Sci, 2024, 15(4): 1449–1471. doi:10.1039/d3sc05552c

[21] Gu S, Yang Y, Zhao Y, et al. Evaluation of AlphaFold2 structures for hit identification across multiple scenarios. J Chem Inf Model, 2024, 64(9): 3630–3639. doi:10.1021/acs.jcim.3c01976

[22] Wong F, Krishnan A, Zheng E J, et al. Benchmarking AlphaFold-enabled molecular docking predictions for antibiotic discovery. Mol Syst Biol, 2022, 18(9): e11081. doi:10.15252/msb.202211081

[23] Pan Z, Yang X, Ge W, et al. A task-adaptive multimodal pretrained framework for antibiotic virtual screening with joint activity and cytotoxicity prediction. J Chem Inf Model, 2026, 66(18): 11793–11806. doi:10.1021/acs.jcim.6c02659

[24] García-Ortegón M, Seal S, Geddes E, et al. Transfer learning enables discovery of sub-micromolar antibacterials for ESKAPE pathogens from ultra-large chemical spaces. Chem Sci, 2025, 16(45): 21518–21533. doi:10.1039/d5sc03055b

[25] Yu Y, Cai C, Wang J, Bo Z, Zhu Z, Zheng H. Uni-Dock: GPU-accelerated docking enables ultralarge virtual screening. J Chem Theory Comput, 2023, 19(11): 3336–3345. doi:10.1021/acs.jctc.2c01145

[26] Brier G W. Verification of forecasts expressed in terms of probability. Mon Weather Rev, 1950, 78: 1–3.

[27] Guo C, Pleiss G, Sun Y, Weinberger K Q. On calibration of modern neural networks. Proceedings of the 34th International Conference on Machine Learning, PMLR 70: 1321–1330, 2017.

[28] Santos-Martins D, Forli S, Ramos M J, Olson A J. AutoDock4Zn: an improved AutoDock force field for small-molecule docking to zinc metalloproteins. J Chem Inf Model, 2014, 54(8): 2371–2379. doi:10.1021/ci500209e

[29] Ganesh H, Mada S, Aryal S, et al. MBind: a web-based platform for metalloprotein and nonmetalloprotein docking with optional ML docking integration. Molecules, 2026, 31(15): 2703. doi:10.3390/molecules31152703

[30] Zhang H, Zhang X, Su Q, et al. MetalloDock: decoding metalloprotein–ligand interactions via physics-aware deep learning for metalloprotein drug discovery. J Am Chem Soc, 2026, 148(3): 3086–3101. doi:10.1021/jacs.5c15876

[31] Balieiro A M, Fernandes G M B, Girdhari L, et al. Why computational predictions fall short in subclass B1 metallo-beta-lactamase inhibitor discovery: progress, pitfalls, and translational opportunities. ACS Infect Dis, 2026, 12(11): 3573–3583. doi:10.1021/acsinfecdis.6c00608

[32] Eberhardt J, Santos-Martins D, Tillack A F, Forli S. AutoDock Vina 1.2.0: new docking methods, expanded force field, and Python bindings. J Chem Inf Model, 2021, 61: 3891–3898. doi:10.1021/acs.jcim.1c00203

[33] Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc B, 1995, 57: 289–300.

[34] Walsh I, Fishman D, Garcia-Gasulla D, et al. DOME: recommendations for supervised machine learning validation in biology. Nat Methods, 2021, 18(10): 1122–1127. doi:10.1038/s41592-021-01205-4

[35] Corso G, Stärk H, Jing B, et al. DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking. International Conference on Learning Representations (ICLR), 2023. arXiv:2210.01776

[36] Jia Y, Gao B, Tan J, et al. Deep contrastive learning enables genome-wide virtual screening. Science, 2026, 391(6781): eads9530. doi:10.1126/science.ads9530

[37] Lomize M A, Lomize A L, Pogozheva I D, Mosberg H I. OPM: orientations of proteins in membranes database. Bioinformatics, 2006, 22: 623–625. doi:10.1093/bioinformatics/btk023

[38] White S H, Wimley W C. Membrane protein folding and stability: physical principles. Annu Rev Biophys Biomol Struct, 1999, 28: 319–365. doi:10.1146/annurev.biophys.28.1.319

[39] Abraham M J, Murtola T, Schulz R, et al. GROMACS: high performance molecular simulations through multi-level parallelism from laptops to supercomputers. SoftwareX, 2015, 1–2: 19–25. doi:10.1016/j.softx.2015.06.001

## Figure captions

**Figure 1. Data landscape of the two-regime Gram-negative panel.** Potency-rich soluble enzymes (Track A) and structure-rich membrane machines (Track B) are shown with public active counts, structure coverage and the two evaluation regimes; the figure motivates the joint panel rather than treating either regime as a substitute for the other.

**Figure 2. Two-regime panel composition and ligand-class breakdown.** Track A targets (7 primary + 2 auxiliary) are shown with curated active counts, docking-box size and observed decoy ratio; the right panel summarises Track A (1,748 curated actives; 28,136 successfully docked ligands) and Track B's 40 co-crystal ligand poses (4 true inhibitors, 24 antibiotic substrates, 3 covalent β-lactams, 8 fragment-like, 1 unresolved/low-confidence ligand (95X; annotated as substrate in Table S18, not counted as an inhibitor); AcrB 55 non-polymer poses; 34 automated inhibitor labels (4 true inhibitors after re-annotation)).

**Figure 3. Track A per-target enrichment.** Heatmap and per-target AUC/EF5 values for the nine targets; auxiliary targets are shown separately and are not included in the metal/non-metal aggregates.

**Figure 4. Exploratory axes.** Temporal (full holdout versus novel scaffold; ΔAUC with cluster bootstrap CI), cross-species (LpxC_Ec versus LpxC_Pa on 951 shared scored decoys; corrected decoy-level CI 0.5033–0.5543) and structure-source (QC-passed holo versus AFDB, n=8; ΔRMSD −3.493, ΔPLIF −0.289) summaries. None is a strict transfer test.

**Figure 5. Three failure modes and benchmark-validity comparison.** Metal-stratified Vina AUC, Track B ligand-class composition (4 inhibitors, 24 antibiotic substrates, 3 covalent β-lactams, 8 fragment-like, 1 unresolved/low-confidence ligand (95X; annotated as substrate in Table S18, not counted as an inhibitor)) and the test-only ECFP4 1-NN versus Vina comparison (median AUC 0.7745; 7/8 comparable targets 2D >= Vina).

**Figure 6. MAES calibration gate and complete-case sensitivity.** Reliability and formal metrics for the pre-registered gate; full sample n=146 (AUROC 0.496; 95% CI 0.318–0.672; |rho| 0.075; Brier 0.145; ECE 0.124) and complete-case n=132 (AUROC 0.434; Brier 0.159; BSS −0.281 full-sample (146 rows, missing RMSD as negative; complete-case −0.290)). The |rho| criterion is decisive; AUROC is inconclusive/underpowered. MAES is retained only as an exploratory/descriptive pose-quality analysis.

**Figure 7. Exploratory method-class comparison (subset-balanced).** (A) Per-target AUC for Vina, Uni-Dock, DiffDock, Boltz-2 and DrugCLIP; (B) pooled AUC with target-level bootstrap CI (main 8 targets; Uni-Dock 7 targets); (C) paired per-target ΔAUC versus Vina; (D) DiffDock single-chain sensitivity for GyrA and FabI (LepB descriptive; 206 attempted/202 ok/4 NA); (E) common seven-target sensitivity means; (F) method coverage (ok rows). Boltz-2 and DrugCLIP are exploratory/pretrained and no leading or complementary-superiority claim is made. All endpoints are subset-balanced and not comparable to full-panel Table 1.