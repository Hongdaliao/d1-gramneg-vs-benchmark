from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
V2 = ROOT / "manuscript/manuscript_v2.md"
OUT = ROOT / "manuscript/v3/manuscript_v3.md"
text = V2.read_text(encoding="utf-8")
refs_start = text.index("## References")
refs_end = text.index("## Figure captions")
refs = text[refs_start:refs_end].rstrip()

new = r'''# Benchmarking structure-based virtual screening across Gram-negative antibacterial targets: target-specific enrichment, co-crystal ligand pose limitations and failure modes

Manuscript version: v3-draft (2026-11-19)

**Running title:** Gram-negative SBVS benchmark and failure modes

**Authors and affiliations:** withheld for review.

**Corresponding author:** correspondence will be handled through the journal submission system.

## Abstract

Structure-based virtual screening (SBVS) is a workhorse of antibacterial lead discovery, but its reliability across Gram-negative targets, and across the two very different data regimes those targets occupy, has not been organised into a common quantitative framework. We analyse two complementary regimes with regime-appropriate endpoints. Track A is a potency-driven nine-target panel (1,748 curated actives; 1,648 primary; 28,136 successfully docked ligands) evaluated by enrichment metrics. Track B is a structure-rich co-crystal ligand pose benchmark: its 40 poses come from five of nine curated membrane-target entries (AcrB 34, PBP3 3, LptB 1, MexB 1, MsbA 1), but only four are true pump inhibitors; 24 are antibiotic substrates, three are covalent β-lactams, eight are fragment-like and one is unresolved. We report per-target enrichment, three exploratory axes, and three benchmark/tool failure modes: metal coordination (target-specific, confounded by receptor and box effects), boron atom typing, and synthetic-decoy bias. The ECFP4 1-NN baseline evaluated on temporal test actives alone reaches median AUC 0.775 (range 0.499–1.000; 5/9 targets <0.80), and 7/8 comparable targets match or exceed Vina, so absolute enrichment is not portable. Cross-species comparison is a paired comparison (ΔAUC +0.5288, corrected decoy-level CI 0.503–0.554), not a transfer test. The structure-source axis, restricted to eight QC-passed pairs, shows ΔRMSD −3.493 Å and ΔPLIF −0.289. The pre-registered MAES calibration gate is NO_GO: |rho| CI upper 0.161 is below 0.40, while AUROC is inconclusive/underpowered (19 positives; complete-case n=132, AUROC 0.434). All data and code will be released upon publication.

## Keywords

Gram-negative bacteria; structure-based virtual screening; benchmark; co-crystal ligand poses; enrichment; decoy bias; ligand class; failure modes

## 1. Introduction

Antimicrobial resistance is a leading global health threat [1], and the WHO bacterial priority pathogen list identifies Gram-negative pathogens as an urgent target for new classes [2]. The lipopolysaccharide transporter LptB2FGC has recently been targeted by two new antibiotic classes [3,4], renewing interest in Gram-negative targets that are structurally tractable but chemically and biologically difficult.

SBVS is the default computational triage for such targets, yet its benchmarks assume both abundant potency data and soluble proteins [5-7]. DUD-E, DEKOIS and LIT-PCBA have been foundational [5-7,30], but they were not designed to represent membrane machines with sparse public activity data, nor to test how performance varies across time, species or receptor source. Recent co-folding and AI-docking benchmarks have expanded the structural side [8-16], including membrane-protein–lipid systems [15,16], but they generally evaluate pose generation or ligand placement rather than jointly organising potency-rich and structure-rich regimes under a common reporting protocol.

Two literatures meet here. On the ligand side, antibiotic virtual screening has moved toward large chemical spaces and multimodal pretraining [26,27]; these approaches are valuable but do not provide a structure-based, per-target assessment of docking reliability. On the structure side, AlphaFold-enabled docking and predicted-structure comparisons have asked whether predicted receptors can support hit identification [31,32], while membrane pose benchmarks such as LiPP and CoMPLip have addressed membrane-protein–ligand placement [15,16]. What remains missing is a common empirical framework that (i) places potency-rich soluble enzymes and potency-poor membrane machines on the same evaluation grid, (ii) measures transfer along explicitly defined axes, and (iii) reports negative as well as positive results with pre-registered gates.

We address this gap with two complementary panels rather than a single quantitative framework. Track A is a potency-driven benchmark of nine targets (seven primary and two auxiliary), comprising 1,748 public actives (1,648 primary) and 28,136 successfully docked ligands; its endpoint is enrichment. Track B is a co-crystal ligand pose benchmark of 40 poses across five of nine curated membrane-target entries; its endpoint is pose geometry and interface recovery, not inhibitor enrichment. The pose set is deliberately reported with its ligand-class composition because it is dominated by antibiotic substrates and fragment-like ligands rather than true pump inhibitors. On this panel we describe three exploratory axes: temporal (time-split holdout versus a novel-scaffold subset), cross-species (a paired comparison of LpxC from E. coli versus P. aeruginosa on 951 shared scored decoys), and structure source (holo versus predicted AFDB receptors, restricted to eight QC-passed pairs). We also report three benchmark/tool failure modes — metal coordination, boron atom typing and DUD-E-style decoy bias — and two pre-registered negative gates: the LIT-PCBA external-consistency anchor and the MAES calibration gate.

Our contributions are fourfold: (i) a two-regime benchmark with regime-appropriate endpoints, not a shared quantitative grid; (ii) exploratory, per-target descriptions of temporal variation, a cross-species paired comparison and a QC-passed structure-source comparison; (iii) three reproducible benchmark/tool failure modes and a test-only 2D baseline that shows synthetic-decoy bias rather than a 2D-versus-docking conclusion; and (iv) a negative-result-driven reporting protocol (pre-registered gates, per-target cluster CIs, observed decoy ratios, no post-hoc threshold or method changes), with MAES retained only as an exploratory/descriptive pose-quality analysis.

## 2. Methods

### 2.1 Panel and data regimes

Track A contains nine Gram-negative antibacterial targets: LpxC (E. coli and P. aeruginosa), DHFR, Gyrase, GyrA, GyrB, FabI, TEM-1 and LepB. Seven targets are primary; GyrA and LepB are auxiliary because of low active counts or low-confidence enrichment. Public activity data were retrieved from ChEMBL 37 and curated with a pChEMBL >= 6 threshold, activity-type harmonisation, salt stripping and duplicate removal. The final panel contains 1,748 curated actives (1,648 primary) and 28,136 successfully docked ligands. Of these, 38 (17 primary: Gyrase_Ec; 21 auxiliary: LepB) fell outside the docking box and were excluded, leaving 1,710 scored actives (1,631 primary) in the analysis table. Receptors were prepared from holo PDB structures with co-crystal ligand centres; box sizes were 26 Å for LpxC_Ec, Gyrase, GyrB and LepB, and 22 Å for the remaining targets. Gyrase_Ec used the GyrB_Ec receptor; the metal primary set therefore involves three unique receptors, not four.

Track B contains 40 co-crystal ligand poses contributed by five of the nine curated membrane-target entries (AcrB 34; PBP3 3; LptB 1; MexB 1; MsbA 1). The panel additionally includes BamA and LptD entries (E. coli and P. aeruginosa) for which no inhibitor-class co-crystal pose passed curation, plus lipid/substrate analogues used only for sensitivity analyses (68 clean poses in total). Ligand classes were reviewed against RCSB chemical-component names: **4 true pump inhibitors** (5QE/5QF/5QG/MBX), **24 antibiotic substrates** (MIY 10, FUA 7, DM2, DXT, ERY, LFX, RBT, CLM, NOV), **3 covalent β-lactams** (AZR/0W0/PFV), **8 fragment-like ligands** (1K8×2, 1KE, A1AN8, A1AOE, A1AOF, A1AON, XE9) and **1 unresolved** (95X). AcrB has 55 non-polymer poses in total, of which 34 are inhibitor-class; Track B is therefore a co-crystal ligand pose benchmark with class stratification, not an inhibitor-only benchmark.

### 2.2 Docking, decoys and enrichment metrics

Docking used AutoDock Vina 1.2.5 [18] with exhaustiveness 4 for EF1-enabled targets (LpxC_Ec, LpxC_Pa and DHFR_Ec) and exhaustiveness 2 for non-EF1 targets; the 1,000-ligand sensitivity arms used exhaustiveness 8 and the LIT-PCBA anchor used exhaustiveness 4. Seed 42 and a single centre/box per receptor were used. Decoys followed the DUD-E-style property-matched protocol [5], with three independent decoy sets per target; the primary analysis used set s1, and sets s2/s3 were used for stability checks. The protocol targeted up to 20 matched decoys per active, but the observed inactive:active ratio was 9.82–19.98 (median 15.15), and the observed ratio is reported per target. Candidate decoys were filtered to ECFP4 max Tanimoto < 0.4 versus all target actives before property matching; this construction is a benchmark-validity feature, not an independent 2D claim.

The primary enrichment metric is EF5% = (hits in top k / n_active) / (k / N), with k = max(1, round(0.05 N)); EF1% is reported only when n_active >= 150. ROC-AUC and BEDROC20 are secondary metrics. Uncertainty is reported with a 2,000-resample active-level cluster bootstrap (decoys grouped by matched active); in GyrA the largest scaffold cluster accounts for about 80.8% of actives and the cluster CI is close to the ligand-level CI. Bootstrap resamples were not clustered by chemical series, which is a limitation. The 2D ECFP4 1-nearest-neighbour baseline was evaluated both on the full ligand set (leaky, because training actives remained in the evaluation set) and on temporal test actives only; the test-only numbers are reported in the main text and the leaky version is retained only in SI as a reproducibility record.

### 2.3 Exploratory axes

The temporal axis uses a time-split design following prior time-split docking evaluations [37,38], with earlier benchmark context [10,31]; test actives were defined by deposition year and a strict novel-scaffold subset required ECFP4 maximum Tanimoto < 0.4 to training and no shared Bemis–Murcko scaffold. Positives were test actives and negatives were the same target decoys; full temporal holdout and novel-scaffold AUC/EF5 were compared with 2,000-resample cluster bootstrap.

The cross-species analysis compares LpxC_Ec and LpxC_Pa on the 951 shared scored decoys as a **paired comparison of two targets/active sets**, not a cross-receptor transfer test. Point estimates and a corrected decoy-level bootstrap (decoys resampled with replacement) are reported; the previous confidence interval computed with `isin()` filtering is retained only as an artefact note in SI. AcrB versus MexB is descriptive only.

The structure-source analysis compares holo and predicted AFDB receptors for **QC-passed pairs only (n=8)**; the all-valid set (n=34) is retained in SI for transparency but contains 26 pairs with `excluded_predicted=True` and is not used as the headline. Predicted receptors are pocket-aligned by 12 Å pocket superposition; RMSD is a ligand-superposition metric. This axis is exploratory because predicted-receptor QC excluded 28/40 systems (21 ligand clashes and 7 TM-score <0.70). H2 is negative/underpowered: the original acceptable-pose-rate endpoint (MAES >= 0.5) difference was 0; the endpoint was replaced by paired RMSD/PLIF after the Amendment #1 trigger, which is disclosed as a protocol deviation.

### 2.4 Track B pose quality and co-crystal ligand composition

Pose quality was evaluated with symmetry-aware ligand RMSD, PLIF/contact recovery and PoseBusters physics/chemistry checks [14]. Only physics/chemistry checks were used; RMSD-like checks were excluded to avoid circularity. Membrane consistency was assessed relative to the OPM reference plane [17], with neutral, flagged values when no reliable plane was available. No molecular-dynamics subset is included in the primary analysis; membrane environment is assessed from reference-relative geometry and physics checks [22,23].

Cognate redocking is explicitly distinguished from cross-docking: each ligand is docked into its own holo receptor, and receptors are monomeric. Interface-bound poses are only 3/40 (source: data/trackB/pose_chain_contacts.csv). The ligand-class composition (4 inhibitors, 24 antibiotic substrates, 3 covalent β-lactams, 8 fragments, 1 unresolved) is itself a benchmark-validity result: a structurally rich membrane-target set is dominated by substrates and fragment-like ligands, so it cannot be described as an inhibitor recovery benchmark.

Receptor-preparation quality is a confounder: `receptor_manifest.csv` records 4/8 receptors with top-1 redock >3 Å (GyrA 11.26 Å, TEM1 11.26 Å, LepB 6.59 Å, GyrB 3.11 Å). The structure-source and metal analyses are interpreted with this limitation.

### 2.5 MAES and the pre-registered W12 calibration gate

MAES combines three pose-level sub-scores: validity V (PoseBusters pass rate plus an optional strain penalty), membrane consistency M (reference-relative membrane-plane geometry) and state consistency S (TM-helix angle, pocket volume and gating distance). The combination rule and sub-score weights were fixed a priori. A candidate probability score was pre-registered, but its calibration gate failed; MAES is reported as a descriptive pose-quality score only.

The W12 gate used 146 pose observations across four method/arm labels: vina 40, vinardo 40, vina_afdb 40 and boltz2 26 (13 systems × 2 diffusion samples; 13 of the 40 inhibitor systems were within the Boltz-2 budget). All phi values came from reference structures, so S was neutral; the primary score is V+M and full V+M+S is exploratory. Labels were strictly separated but the calibration label includes PoseBusters while V is derived from PoseBusters; this circularity and the nesting of calibration RMSD<2 Å with evaluation RMSD<1 Å are disclosed. Evaluation used five-fold grouped cross-validation by system (n_train 116–118; n_test 28–30), wild cluster bootstrap for AUROC (4,999 draws; 40 clusters) and cluster bootstrap for rho (2,000 draws). The pre-registered Go criteria were AUROC >= 0.65 and |rho(MAES, RMSD)| >= 0.40. No post-hoc threshold, method or sample change was made.

### 2.6 LIT-PCBA external-consistency check

The LIT-PCBA anchor [7] contains 15 target sets. Amendment A2 used a 20:1 inactive-to-active ratio, a cap of 1,000 inactives per set, three inactive subsample seeds, Vina exhaustiveness 4, seed 42 and all 15 sets attempted. AC2 was evaluated on the nine sets with n_val_act >= 20; AC3 used the coefficient of variation across three inactive subsamples. EF1% was reported only for ALDH1 and VDR; ADRB2 was descriptive only. The anchor was CPU-only (AC5 = not applicable) and is an appendix external-consistency check, not a validation claim. Because of the 1,000-inactive cap, ALDH1 had an inactive:active ratio of 0.74:1 instead of 20:1; this is disclosed.

### 2.7 Statistics and reproducibility

All confidence intervals use percentile bootstrap unless stated otherwise; main Track A cluster CIs use 2,000 resamples, W12 wild-cluster AUROC uses 4,999 Rademacher draws and W12 rho uses 2,000 cluster bootstrap resamples. Calibration was assessed with Brier score [24] and expected calibration error [25]; reporting follows DOME recommendations [28]. The per-target p-values in the enrichment JSON are labelled invalid under matched-decoy clustering; Benjamini–Hochberg values are not used for inference. We report MDES/descriptive power rather than post-hoc power: Track A detects approximately 0.93–0.99 × between-target SD; Track B effective n_eff is 20–35 and detects only >=15–25 percentage-point acceptable-rate differences, or MAES differences of about 0.10–0.16. Code, manifest hashes and intermediate tables are recorded with the release; the W12 manifest is frozen at sha256 e4ae5c23461474f773a016ec415e9da8e5dd6ff10bbef55eda73c726cb4fa3ca.

## 3. Results

### 3.1 Two complementary panels

**Figure 1** shows the data landscape and **Figure 2** summarises panel composition. Track A contains 1,748 curated actives (1,648 primary) and 28,136 successfully docked ligands; 1,710 actives (1,631 primary) remain after box-clipping. Track B contains 40 co-crystal ligand poses from five of the nine curated membrane-target entries (AcrB 34; PBP3 3; LptB 1; MexB 1; MsbA 1). The category composition is 4 true inhibitors, 24 antibiotic substrates, 3 covalent β-lactams, 8 fragment-like ligands and 1 unresolved ligand. The two panels answer different questions: enrichment for Track A and pose/interface recovery for Track B; no shared quantitative endpoint is claimed.

### 3.2 Target-specific enrichment in Track A

**Table 1 and Figure 3** summarise per-target Vina enrichment. The seven primary targets span AUC 0.300 (LpxC_Ec) to 0.852 (LpxC_Pa), with no simple metal/non-metal dichotomy. The metal primary set (n = 4: GyrB, Gyrase, LpxC_Ec, LpxC_Pa) has mean AUC 0.521, median 0.466 and mean EF5% 1.780; the non-metal set (n = 3: DHFR, FabI, TEM1) has mean AUC 0.567, median 0.606 and mean EF5% 1.028. GyrA and LepB are auxiliary and were not included in either aggregate. Because Gyrase_Ec and GyrB_Ec share a receptor, the metal aggregate involves three receptors; this is descriptive and not a receptor-replication test. **Table 1** reports 1,710 scored actives (1,631 primary); 38 actives (17 primary: Gyrase_Ec; 21 auxiliary: LepB) were excluded because they fell outside the docking box. The LpxC homologues illustrate target specificity: LpxC_Ec (Zn) has AUC 0.300, whereas LpxC_Pa (Zn) has AUC 0.852. This difference persisted in property-stratified and nearest-decoy matched analyses and was not explained by chemotype concentration after excluding the largest Pa scaffolds; receptor-preparation quality remains a confounder, not an excluded explanation.

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

Of 1,748 curated actives (1,648 primary), 38 (17 primary: Gyrase_Ec; 21 auxiliary: LepB) fell outside the docking box and were excluded, leaving 1,710 scored actives (1,631 primary) in Table 1. Auxiliary values (GyrA_Ec and LepB_Ec) are descriptive and are not included in any aggregate or trend.

### 3.3 Benchmark validity: test-only 2D baseline and decoy bias

The original ECFP4 1-NN analysis retained training actives in the evaluation set and reached median AUC 0.958. After restricting evaluation to temporal test actives, the test-only baseline reaches **median AUC 0.7745 (range 0.4992–1.0000)**; 5/9 targets are below 0.80 (TEM1 0.499, GyrB 0.572, Gyrase 0.653, GyrA 0.705, FabI 0.775). On 7/8 comparable targets, the test-only 2D baseline matches or exceeds Vina. This is not evidence that docking is useless; it shows that synthetic decoys filtered to ECFP4 max Tanimoto < 0.4 versus actives make the 2D task artificially easy, so absolute enrichment values are not portable to real screening. The leaky full-set numbers are retained only in SI as a reproducibility record, and the original “0.82–1.00 / near-ceiling” statement is withdrawn.

### 3.4 Sensitivities and robustness

The 1,000-ligand LpxC_Ec box-size three-tier check gave AUC 0.224, 0.289 and 0.333 for 22, 26 and 28 Å; ranking was stable (Spearman rho = 0.949 for 26 versus 28 Å), but the absolute box effect (0.109) is larger than the primary metal/non-metal mean difference (0.046). The 26 Å box is retained as primary because it was selected by a ligand-size audit; 22 Å, which truncates the largest ligands, is an informative sensitivity. Exhaustiveness 4 versus 8 and decoy-ratio 10/15/20:1 checks showed small absolute changes (|Delta AUC| <= 0.0116; Tables S4b and S8). Applying Vinardo as a second scoring function changed absolute enrichment (for example, LpxC_Ec 0.300 to 0.554; LpxC_Pa 0.852 to 0.918; Gyrase 0.578 to 0.518) while rankings remained moderately concordant (Spearman rho 0.65–0.84). A single scoring function is not the ceiling of physics-based scoring, and these analyses remain in SI rather than the primary claim. The Zinc-stripped LpxC_Ec control gave AUC 0.220 versus 0.222 for Zn-retained, i.e. no meaningful change; this argues against a simple Zn-coordination explanation but does not identify the cause.

### 3.5 Track B: co-crystal ligand poses and structure-source comparison

Track B cognate redocking recovered 9/36 poses below 2 Å (median RMSD 4.44 Å), with median PLIF recovery 0.896; strict PoseBusters physics/chemistry checks gave 32 passed, 4 failed and 4 not assessable (no reference coordinates). The ligand-class composition is 4 true inhibitors, 24 antibiotic substrates, 3 covalent β-lactams, 8 fragment-like ligands and 1 unresolved; AcrB contributes 34 of 55 non-polymer poses but only 34 inhibitor-class poses. The AFDB predicted-receptor arm passed QC for 12/40 systems and excluded 28/40. For **QC-passed pairs only (n=8)**, holo-to-AFDB mean RMSD changed from 5.24 to 1.74 Å (Δ **−3.493 Å**, 95% CI −4.711 to −2.064) and mean PLIF recovery dropped from 0.818 to 0.528 (Δ **−0.289**, 95% CI −0.418 to −0.180); PoseBusters passed 7/1 for holo and 7/1 for AFDB. The all-valid n=34 set (ΔRMSD −2.436 Å; ΔPLIF −0.519) contains 26 `excluded_predicted=True` pairs and is retained only in SI. H2 is negative/underpowered: the original acceptable-pose-rate endpoint difference was 0, and the endpoint was replaced after Amendment #1; predicted receptors do not support a confirmatory transfer claim.

### 3.6 W12 calibration gate: NO_GO

**Figure 6** displays the calibration-gate result. The pre-registered gate used 146 pose observations (vina 40, vinardo 40, vina_afdb 40, boltz2 26; 19 pose-level evaluation positives; 40 systems). The primary V+M score gave AUROC 0.496 (wild cluster CI 0.318–0.672), signed rho −0.0745 (|rho| 0.075; cluster CI −0.302 to 0.161), Brier 0.145 and ECE 0.124. The **|rho| criterion is decisive** (upper CI 0.161 < 0.40). The **AUROC criterion is inconclusive**: 0.496 is below 0.65, but the CI includes 0.65 and only 19 positives give low power. Complete-case analysis excluding 14 non-assessable rows (n=132) gives AUROC 0.434 and Brier 0.159 (BSS −0.28), versus 0.496/0.145 in the full sample; the full-sample result treats missing rows as negative. Stratified AUROC estimates reverse direction (vina 0.572, vina_afdb 0.454; holo 0.572, predicted_afdb 0.454; PBP3 0.107 with rho +0.811). MAES is therefore retained as an exploratory/descriptive pose-quality analysis, not as a calibrated probability or a demonstrated failure of the AUROC criterion.

**Table 2. W12 calibration gate (frozen plan; primary V+M).** Manifest sha256 e4ae5c23461474f773a016ec415e9da8e5dd6ff10bbef55eda73c726cb4fa3ca. CIs are wild cluster bootstrap (AUROC) and cluster bootstrap (rho), with 40 clusters.

| Analysis | n | Evaluation positives | AUROC (95% CI) | |rho| (95% CI) | Brier | ECE |
|---|---:|---:|---|---|---:|---:|
| Primary V+M | 146 | 19 | 0.496 (0.318–0.672) | 0.075 (−0.302 to 0.161) | 0.145 | 0.124 |
| Complete-case V+M | 132 | 19 | 0.434 (not bootstrapped; reviewer-authoritative) | — | 0.159 | — |
| Exploratory V+M+S | 146 | 19 | 0.495 (0.318–0.671) | 0.072 (−0.301 to 0.166) | 0.145 | 0.124 |

### 3.7 Exploratory axes

**Figure 4** summarises the exploratory axes. **Temporal.** Across the nine targets, the novel-scaffold subset did not uniformly degrade temporal enrichment: median ΔAUC was +0.019, with Gyrase_Ec showing a clear decrease (−0.147, 95% CI −0.212 to −0.089) and DHFR_Ec showing an increase (+0.084, 95% CI +0.031 to +0.143); other confidence intervals included zero. Only two of nine targets had CIs excluding zero, so the temporal axis is descriptive/mixed and is not a transfer test.

**Cross-species.** On 951 shared scored decoys, LpxC_Ec gave AUC 0.281 (95% CI 0.265–0.297) and LpxC_Pa gave 0.810 (0.798–0.821); the paired ΔAUC was +0.5288. The corrected decoy-level 95% CI is **0.5033–0.5543** (active-level 0.4961–0.5601); the previously reported 0.5135–0.5433 interval was an `isin()` de-duplication artefact and is withdrawn. This is a paired comparison of two targets/active sets, not a cross-receptor transfer test; it shows that the two homologues differ strongly in discriminability.

**Structure source.** Using QC-passed pairs only (n=8), holo-to-AFDB mean RMSD changed from 5.24 to 1.74 Å (Δ−3.493, 95% CI −4.711 to −2.064) and PLIF from 0.818 to 0.528 (Δ−0.289, 95% CI −0.418 to −0.180). The all-valid n=34 version (ΔRMSD −2.436; ΔPLIF −0.519) is SI-only and contains 26 QC-excluded pairs. H2 is negative/underpowered; the original acceptable-rate endpoint difference was 0. The axis is exploratory and does not support a geometry-improves/interface-degrades headline.

**Consistency.** The three axes are descriptive and axis-dependent; none is a strict transfer test. The data do not support a single transferability score. **Table 3** summarises the corrected descriptive values.

**Table 3. Exploratory axis summary.**

| Axis | Key metric | Interpretation |
|---|---|---|
| Temporal (novel scaffold) | median ΔAUC +0.019; Gyrase_Ec −0.147; DHFR_Ec +0.084; 2/9 CIs exclude 0 | descriptive/mixed |
| Cross-species (LpxC Ec vs Pa) | ΔAUC +0.5288; corrected decoy-level CI 0.503–0.554 | paired comparison, not transfer |
| Structure source (QC-passed n=8) | ΔRMSD −3.493 (CI −4.711 to −2.064); ΔPLIF −0.289 (CI −0.418 to −0.180) | exploratory; H2 negative |

### 3.8 LIT-PCBA external-consistency check

The anchor was run on all 15 LIT-PCBA sets (15/15 attempted and completed). AC1 passed; AC2 was 2/9 = 22.2% of the pre-registered evaluable sets (GBA and KAT2A only) and therefore below the 30% threshold, but this result is fragile: including TP53 (16 actives) gives 3/10 = 30.0%, exactly at the threshold (Wilson 95% CI 6.3–54.7%). AC3 failed using either definition, with maximum CV 28.3% under population SD and **34.64%** under sample SD (MTORC1); AC4, AC6 and AC7 passed; AC5 was not applicable because the anchor was CPU-only. The 15-set descriptive view gave EF5% >= 2 for 5/15 = 33.3%, median AUC 0.547 and median EF5% 1.19; EF1% was reported only for ALDH1 (1.21) and VDR (1.28). The primary AC2 judgement remains the pre-registered nine evaluable sets (2/9 = 22.2%), but the manuscript states that the 30% threshold was not robustly met rather than claiming a stable FAIL. The ALDH1 cap gave an inactive:active ratio of 0.74:1. The anchor is an appendix external-consistency check, not a validation claim.

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

### 3.9 Failure modes and benchmark validity

**Figure 5** shows the three failure modes and the test-only 2D comparison. First, metal coordination remains a known limitation [33-36], but the panel shows target-specific behaviour (LpxC_Ec 0.300 versus LpxC_Pa 0.852) and is confounded by receptor identity (Gyrase=GyrB), receptor-preparation quality (4/8 top-1 redock >3 Å) and box size (0.109 effect versus 0.046 metal/non-metal mean difference). The primary aggregate (n=4 metal mean 0.521, median 0.466, EF5 1.780; n=3 non-metal mean 0.567, median 0.606, EF5 1.028) is descriptive; the alternative by-stratum aggregate (metal n=5 mean 0.588; non-metal n=4 mean 0.498) reverses direction and is disclosed only as a caveat. Metal is therefore not claimed as an independent effect. Second, boron atom typing is missing in Vina 1.2.5; 25 boron-containing actives were handled with a B-to-C isosteric surrogate and an exclude-surrogate sensitivity is reported. Third, synthetic decoys filtered to ECFP4 max Tanimoto < 0.4 versus actives make the test-only 2D baseline competitive with Vina on 7/8 comparable targets; this is a benchmark-construction finding, not evidence that docking is unnecessary. Absolute EF/AUC values are not portable.

## 4. Discussion

Our principal finding is that enrichment is target-specific and that the benchmark's quantitative claims must be separated from its construction artefacts. The seven primary Track A targets span AUC 0.300–0.852; the metal primary set (n=4; mean 0.521, median 0.466, EF5 1.780) overlaps the non-metal set (n=3; mean 0.567, median 0.606, EF5 1.028). The LpxC paralogues illustrate the point directly: LpxC_Ec is inverted whereas LpxC_Pa is strongly enriched. This difference survives property matching and chemotype exclusion, but receptor-preparation quality remains a confounder, so it should not be assigned a single mechanism.

The two panels use different endpoints. Track A measures enrichment; Track B measures co-crystal ligand pose and interface recovery. Track B is not an inhibitor benchmark: among 40 poses, only four are true pump inhibitors, 24 are antibiotic substrates, three are covalent β-lactams, eight are fragment-like and one is unresolved; AcrB contributes 34 of 55 non-polymer poses but only 34 inhibitor-class poses. Reporting this composition is itself a benchmark-validity result: structure-rich membrane targets can be poor in inhibitor-class co-crystal poses.

The exploratory axes are descriptive. Temporal variation is mixed (2/9 CIs exclude zero). The cross-species analysis is a paired comparison of two homologues (ΔAUC +0.5288; corrected CI 0.503–0.554), not a transfer test. The structure-source comparison is restricted to eight QC-passed pairs (ΔRMSD −3.493; ΔPLIF −0.289); the all-valid n=34 version contains 26 QC-excluded pairs and is SI-only. H2 is negative/underpowered and its endpoint was replaced after Amendment #1. None of these axes supports a single transferability score.

Benchmark validity and tool limitations carry the main methodological message. The test-only ECFP4 baseline reaches median AUC 0.7745 (range 0.499–1.000; 5/9 <0.80) and matches or exceeds Vina on 7/8 comparable targets after synthetic-decoy filtering; the original full-set near-ceiling result was inflated by training actives in evaluation. Metal coordination [33-36], boron atom typing and receptor/box confounders are reported alongside it. The correct conclusion is that absolute enrichment is not portable and that 2D competitiveness under synthetic decoys is a warning about benchmark construction, not evidence for a new method.

The negative pre-registered gates are part of the result. The LIT-PCBA anchor passes coverage, freeze and reproducibility checks but does not robustly reach AC2 (2/9; 3/10 with TP53) and fails AC3 under population and sample SD (28.3%/34.64%). The W12 gate is NO_GO because the |rho| criterion is decisive (upper CI 0.161 < 0.40); the AUROC criterion is inconclusive/underpowered (19 positives; complete-case n=132, AUROC 0.434; BSS −0.28). No thresholds, methods or observations were changed after observing these results; MAES is retained only as a descriptive pose-quality analysis.

For practice, the results argue for per-target reporting with cluster-aware intervals, observed decoy ratios, test-only 2D baselines, and explicit metal/receptor/box coverage. A single docking score should not be treated as a portable predictor of antibacterial enrichment; target-specific validation and failure-mode audits remain necessary.

## 5. Limitations

This study has no wet-lab validation. Track A is limited to nine targets and 1,748 curated actives (1,710 scored after box-clipping), with scored actives per target ranging from 6 to 418 after exclusion. Track B is cognate redocking into monomeric receptors; it is a co-crystal ligand pose benchmark, not an inhibitor benchmark: only 4/40 poses are true inhibitors, 24 are antibiotic substrates, 3 are covalent β-lactams, 8 are fragment-like and 1 is unresolved; only 3/40 are interface-bound. The structure-source axis is restricted to 8 QC-passed pairs; the all-valid n=34 set includes 26 QC-excluded pairs. H2 is negative/underpowered and its endpoint was replaced after Amendment #1; H1 and H3 were not reported. The LIT-PCBA anchor was CPU-only under a reduced inactive ratio, with ALDH1 at 0.74:1, and AC2/AC3 are not robust to reasonable sensitivity choices. The W12 MAES gate is NO-GO on |rho| and inconclusive on AUROC; MAES is descriptive only. Metal analyses are confounded by Gyrase=GyrB receptor sharing, receptor-preparation quality (4/8 top-1 redock >3 Å), box-size effects (0.109 versus a 0.046 metal/non-metal difference) and auxiliary-target inclusion; the by-stratum aggregate reverses direction and is disclosed only as a caveat. Decoy ratios were 9.82–19.98 rather than the planned 20:1; decoys were filtered to ECFP4 max Tanimoto <0.4; non-EF1 targets used exhaustiveness 2; bootstrap counts vary by analysis (2,000 main; 4,999 W12 wild; 2,000 W12 rho); the per-target p-values are invalid for inference; Benjamini–Hochberg values are not used; and per-target MDES/descriptive power is reported instead of post-hoc power. The 16 GB GPU limits co-folding to 13 of the 40 inhibitor systems.

## 6. Conclusion

Across two complementary Gram-negative panels with regime-appropriate endpoints, enrichment is target-specific, the co-crystal ligand pose set is dominated by antibiotic substrates and fragment-like ligands rather than true inhibitors (4/24/3/8/1), and the exploratory temporal, cross-species and structure-source axes are descriptive rather than transfer tests. The test-only 2D baseline (median AUC 0.7745; 7/8 comparable targets matching or exceeding Vina) shows that synthetic-decoy construction inflates absolute enrichment; metal coordination, boron atom typing and decoy bias are reported as benchmark/tool failure modes. Both pre-registered gates are negative or inconclusive: the LIT-PCBA anchor does not robustly reach AC2 and fails AC3 under both SD definitions, and the MAES W12 gate fails the |rho| criterion while AUROC is underpowered/inconclusive (complete-case AUROC 0.434; BSS −0.28). MAES remains descriptive only. These results motivate per-target, uncertainty-aware reporting rather than a single ranking or probability score, and they do not support claims of a shared evaluation grid or a single transferability metric.

## Data and code availability

All curated data, decoy manifests, docking scores, analysis scripts and figure-generation code will be released upon publication. The W12 calibration manifest and report are identified by sha256 e4ae5c23461474f773a016ec415e9da8e5dd6ff10bbef55eda73c726cb4fa3ca. Intermediate tables (including per-target enrichment, sensitivity analyses, anchor ACs, W12 metrics, ligand-class review, QC sensitivity and corrected cross-species intervals) are provided as SI tables S1–S24, including sub-tables S4a/S4b, S13b, S16b, S19b, S22b and S22c.

''' + refs + r'''

## Figure captions

**Figure 1. Data landscape of the two-regime Gram-negative panel.** Potency-rich soluble enzymes (Track A) and structure-rich membrane machines (Track B) are shown with public active counts, structure coverage and the two evaluation regimes; the figure motivates the joint panel rather than treating either regime as a substitute for the other. Source: `figures/fig1_data_landscape.png`.

**Figure 2. Two-regime panel composition and ligand-class breakdown.** Track A targets (7 primary + 2 auxiliary) are shown with curated active counts, docking-box size and observed decoy ratio; the right panel summarises Track A (1,748 curated actives; 28,136 successfully docked ligands) and Track B's 40 co-crystal ligand poses (4 true inhibitors, 24 antibiotic substrates, 3 covalent β-lactams, 8 fragment-like, 1 unresolved; AcrB 34/55 non-polymer poses). Source: `figures/trackA_enrichment_full.json`, `data/target_panel.csv`, `data/trackB/ligand_class_review.csv`; figure file `figures/v3/fig2_benchmark_design.png` / `.pdf`.

**Figure 3. Track A per-target enrichment.** Heatmap and per-target AUC/EF5 values for the nine targets; auxiliary targets are shown separately and are not included in the metal/non-metal aggregates. Source: `figures/fig3_trackA_enrichment_full.png`.

**Figure 4. Exploratory axes.** Temporal (full holdout versus novel scaffold; ΔAUC with cluster bootstrap CI), cross-species (LpxC_Ec versus LpxC_Pa on 951 shared scored decoys; corrected decoy-level CI 0.5033–0.5543) and structure-source (QC-passed holo versus AFDB, n=8; ΔRMSD −3.493, ΔPLIF −0.289) summaries. None is a strict transfer test. Source: `results/trackA/axes/temporal.csv`, `cross_species_bootstrap_corrected.json`, `structure_source_qcpassed.csv`; figure file `figures/v3/fig4_three_axis_transfer.png` / `.pdf`.

**Figure 5. Ligand-class composition and test-only 2D comparison.** Track B poses are classified as 4 inhibitors, 24 antibiotic substrates, 3 covalent β-lactams, 8 fragment-like and 1 unresolved; the right panel compares the test-only ECFP4 1-NN baseline with Vina (median AUC 0.7745; 7/8 comparable targets 2D >= Vina). Source: `data/trackB/ligand_class_review.csv`, `results/trackA/ecfp4_baseline/baseline_metrics_test_only.csv`, `figures/trackA_enrichment_full.json`; figure file `figures/v3/fig5_ligand_class_2d.png` / `.pdf`.

**Figure 6. W12 calibration gate and complete-case sensitivity.** Reliability and formal metrics for the pre-registered gate; full sample n=146 (AUROC 0.496; 95% CI 0.318–0.672; |rho| 0.075; Brier 0.145; ECE 0.124) and complete-case n=132 (AUROC 0.434; Brier 0.159; BSS −0.28). The |rho| criterion is decisive; AUROC is inconclusive/underpowered. MAES is retained only as an exploratory/descriptive pose-quality analysis. Source: `figures/v3/fig6_w12_calibration_formal.png` / `.pdf`.
'''
OUT.write_text(new, encoding="utf-8")
print("wrote", OUT, "chars", len(new))
