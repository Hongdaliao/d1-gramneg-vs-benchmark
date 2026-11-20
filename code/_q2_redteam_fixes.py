from pathlib import Path

p = Path("manuscript/manuscript_v2.md")
t = p.read_text(encoding="utf-8")

# B1 Track B target list
old = "Track B contains 40 inhibitor-class cognate-redocking poses across nine membrane targets (AcrB, MexB, MsbA, LptB, LptD, BamA, PBP3, LpxC and MraY), plus lipid/substrate analogues for sensitivity analyses (68 clean poses in total). Each co-crystal ligand was docked into its own holo receptor (cognate redocking), not into a different receptor structure."
new = "Track B contains 40 inhibitor-class cognate-redocking poses contributed by five of the nine curated membrane-target entries (AcrB, 34; PBP3, 3; LptB, 1; MexB, 1; MsbA, 1). The panel additionally includes BamA and LptD entries (E. coli and P. aeruginosa) for which no inhibitor-class co-crystal pose passed curation, plus lipid/substrate analogues used only for sensitivity analyses (68 clean poses in total). Each co-crystal ligand was docked into its own holo receptor (cognate redocking), not into a different receptor structure."
assert old in t
t = t.replace(old, new)

old = "Track B is a structure-driven cognate-redocking benchmark of 40 inhibitor poses across nine membrane targets, with 68 clean poses available for sensitivity analyses."
new = "Track B is a structure-driven cognate-redocking benchmark of 40 inhibitor poses contributed by five of the nine curated membrane-target entries (AcrB, PBP3, LptB, MexB and MsbA), with 68 clean poses available for sensitivity analyses."
assert old in t
t = t.replace(old, new)

old = "Track B contributes 40 inhibitor cognate-redocking poses across nine membrane targets and 68 clean poses for sensitivity analyses."
new = "Track B contributes 40 inhibitor cognate-redocking poses from five of the nine curated membrane-target entries (AcrB, 34; PBP3, 3; LptB, 1; MexB, 1; MsbA, 1) and 68 clean poses for sensitivity analyses."
assert old in t
t = t.replace(old, new)

# B2 Table 2 exploratory CI
old = "| Exploratory V+M+S | 146 | 19 | 0.495 (0.308–0.665) | 0.072 (−0.226 to 0.285) | 0.145 | 0.124 |"
new = "| Exploratory V+M+S | 146 | 19 | 0.495 (0.318–0.671) | 0.072 (−0.301 to 0.166) | 0.145 | 0.124 |"
assert old in t
t = t.replace(old, new)
old = "**Table 2. W12 calibration gate (frozen plan; primary V+M).** Manifest sha256 e4ae5c23461474f773a016ec415e9da8e5dd6ff10bbef55eda73c726cb4fa3ca."
new = "**Table 2. W12 calibration gate (frozen plan; primary V+M).** Manifest sha256 e4ae5c23461474f773a016ec415e9da8e5dd6ff10bbef55eda73c726cb4fa3ca. CIs are wild cluster bootstrap (AUROC) and cluster bootstrap (rho), with 40 clusters."
assert old in t
t = t.replace(old, new)

# M2 Abstract curated vs scored range
old = "a potency-driven regime of soluble enzymes with 27–418 public actives (LpxC, DHFR, GyrA/B, FabI, LepB, TEM-1)"
new = "a potency-driven regime of soluble enzymes with 27–418 curated actives (6–418 after box-clipping; LpxC, DHFR, GyrA/B, FabI, LepB, TEM-1)"
assert old in t
t = t.replace(old, new)

# M2 limitation range
old = "The Track A panel is limited to nine targets and 1,748 actives, and the number of public actives varies from 6 to 418."
new = "The Track A panel is limited to nine targets and 1,748 curated actives; the number of scored actives per target is 6–418 after box-clipping."
assert old in t
t = t.replace(old, new)

# M3 Vinardo coverage
old = "For LpxC_Ec and DHFR, a subset of ligands was also rescored with Vinardo on the same Vina poses; this is a second scoring function, not an independent docking engine."
new = "For LpxC_Ec, LpxC_Pa, Gyrase_Ec and DHFR_Ec, ligands were rescored with Vinardo on the same Vina poses (score_only; second scoring function, not an independent engine)."
assert old in t
t = t.replace(old, new)

# M4 membrane MD claim -> not included
old = "Membrane consistency was assessed relative to the OPM reference plane [17], with neutral, flagged values when no reliable plane was available. A 100 ns GROMACS membrane-MD subset [22,23] was used only as a secondary membrane-environment sensitivity and is not part of the primary claim."
new = "Membrane consistency was assessed relative to the OPM reference plane [17], with neutral, flagged values when no reliable plane was available. No molecular-dynamics subset is included in the primary analysis; membrane environment is assessed from reference-relative geometry and physics checks [22,23]."
assert old in t
t = t.replace(old, new)

# M5 time-split citations
old = "The temporal axis uses the time-split protocol used in prior docking benchmarks [8,10,11,31], with test actives defined by deposition year"
new = "The temporal axis uses a time-split design following prior time-split docking evaluations [37,38], with earlier benchmark context [10,31], and test actives defined by deposition year"
assert old in t
t = t.replace(old, new)

# M6 PoseBusters wording
old = "strict PoseBusters physics/chemistry checks passed 32/4/4."
new = "strict PoseBusters physics/chemistry checks gave 32 passed, 4 failed and 4 not assessable (no reference coordinates)."
assert old in t
t = t.replace(old, new)

# m1 source for 3/40
old = "Interface-bound poses are only 3/40, so the benchmark measures cognate pose recovery rather than generalisation to new pockets."
new = "Interface-bound poses are only 3/40 (source: data/trackB/pose_chain_contacts.csv), so the benchmark measures cognate pose recovery rather than generalisation to new pockets."
assert old in t
t = t.replace(old, new)

# m2 55 poses wording
old = "with up to 55 curated ligand-bound poses per target (AcrB, from multiple PDB co-crystals)"
new = "with up to 55 curated inhibitor poses in the most pose-rich target (AcrB; multiple PDB co-crystals)"
assert old in t
t = t.replace(old, new)

# m4 decoy ratio SI refs
old = "Exhaustiveness 4 versus 8 and decoy-ratio 10/15/20:1 checks showed small absolute changes (|Delta AUC| <= 0.0116)."
new = "Exhaustiveness 4 versus 8 and decoy-ratio 10/15/20:1 checks showed small absolute changes (|Delta AUC| <= 0.0116; Tables S4b and S8)."
assert old in t
t = t.replace(old, new)

# m5 data availability sub-tables
old = "Intermediate tables (including per-target enrichment, sensitivity analyses, anchor ACs and W12 metrics) are provided as SI tables S1–S17."
new = "Intermediate tables (including per-target enrichment, sensitivity analyses, anchor ACs and W12 metrics) are provided as SI tables S1–S17 (including sub-tables S4a/S4b, S13b and S16b)."
assert old in t
t = t.replace(old, new)

p.write_text(t, encoding="utf-8")
print("v2 red-team text fixes applied")
