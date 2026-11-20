from pathlib import Path

p = Path("manuscript/manuscript_v2.md")
t = p.read_text(encoding="utf-8")

old_abs = "Here we jointly analyse Gram-negative antibacterial targets across two complementary data regimes: a potency-driven regime of soluble enzymes with 27–418 public actives (LpxC, DHFR, GyrA/B, FabI, LepB, TEM-1) and a structure-rich regime of membrane machines with up to 55 ligand-bound structures per target but few public actives (MsbA, AcrB, BamA, LptD, LptB, MexB, PBP3). Under one protocol, we assess docking (Vina, Uni-Dock), contrastive screening (DrugCLIP) and co-folding (Boltz-2 within our hardware budget), and we quantify transferability along three axes on the same panel: temporal (time-split), cross-species (Escherichia coli versus Pseudomonas aeruginosa homologues) and structure source (holo, apo and predicted receptors)."
new_abs = "Here we jointly analyse Gram-negative antibacterial targets across two complementary data regimes: a potency-driven regime of soluble enzymes with 27–418 public actives (LpxC, DHFR, GyrA/B, FabI, LepB, TEM-1) and a structure-rich regime of membrane machines with up to 55 curated ligand-bound poses per target (AcrB, from multiple PDB co-crystals) but few public actives (MsbA, AcrB, BamA, LptD, LptB, MexB, PBP3). Under one protocol, we assess Vina docking with Vinardo as a second scoring function and an ECFP4 1-NN baseline, and use Boltz-2 co-folding only for 13 systems in the W12 calibration arm; Uni-Dock and DrugCLIP were not run because of environment/compute constraints. We quantify transferability along three axes on the same panel: temporal (time-split), cross-species (Escherichia coli versus Pseudomonas aeruginosa homologues) and structure source (holo versus predicted AFDB receptors; apo receptors were prepared but not evaluated)."
assert old_abs in t
t = t.replace(old_abs, new_abs)

old_exh = "Docking used AutoDock Vina 1.2.5 [18] with exhaustiveness 8 for the main panel (4 for sensitivity arms), seed 42 and a single centre/box per receptor."
new_exh = "Docking used AutoDock Vina 1.2.5 [18] with exhaustiveness 4 for EF1-enabled targets and exhaustiveness 2 for non-EF1 targets in the main panel, exhaustiveness 8 for the 1,000-ligand sensitivity arms (LpxC_Ec and DHFR), seed 42 and a single centre/box per receptor; the LIT-PCBA anchor used exhaustiveness 4."
assert old_exh in t
t = t.replace(old_exh, new_exh)

# Methods: state apo not evaluated explicitly in the structure-source paragraph
old_struct = "The structure-source axis compares holo and predicted AFDB receptor arms for the 34 paired Track B poses with valid paired RMSD; predicted receptors are pocket-aligned by 12 Å pocket superposition, and RMSD is a ligand-superposition metric."
new_struct = "The structure-source axis compares holo and predicted AFDB receptor arms for the 34 paired Track B poses with valid paired RMSD; predicted receptors are pocket-aligned by 12 Å pocket superposition, and RMSD is a ligand-superposition metric. Apo receptors were prepared but were not evaluated, so the axis is holo versus predicted only."
assert old_struct in t
t = t.replace(old_struct, new_struct)

p.write_text(t, encoding="utf-8")
print("v2 claim fixes applied")
