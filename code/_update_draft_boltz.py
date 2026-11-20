from pathlib import Path

p = Path(r"manuscript/v3/method_axis_integration_draft.md")
t = p.read_text(encoding="utf-8")
repls = [
    ("> Status: **draft for integration after Boltz-2 is complete**; do not paste into v3 until numbers are final.",
     "> Status: **Boltz-2 complete; ready for v3 integration after Lead diff review**; do not paste into v3 until the diff is checked."),
    ("The existing ECFP4 1-NN baseline (test-only version), used as a 2D reference, not as a new method.",
     "The existing ECFP4 1-NN baseline was not recomputed on this subset because of temporal-train overlap (see Addendum 12a); the 2D reference remains the full-panel test-only result in §3.3."),
    ("(iv) the existing ECFP4 1-NN baseline is **not** recomputed on this subset, because 68–94% of subset actives overlap its temporal-train reference for the eight comparison targets (GyrB 68% to LpxC 94%; TEM1 38%, LepB 100% at n=3), so a subset-specific ECFP4 AUC would be a self-similarity artefact.",
     "(iv) the existing ECFP4 1-NN baseline is **not** recomputed on this subset, because 68–94% of subset actives overlap its temporal-train reference for seven of eight comparison targets (GyrB 68% to LpxC 94%; TEM1 38%, LepB 100% at n=3), so a subset-specific ECFP4 AUC would be a self-similarity artefact."),
    ("ECFP4:** not recomputed on this balanced subset because 68–94% of subset actives overlap its temporal-train reference for the eight comparison targets (TEM1 38% and LepB 100% at n=3 are exceptions; full per-target overlap in the S25 footnote)",
     "ECFP4:** not recomputed on this balanced subset because 68–94% of subset actives overlap its temporal-train reference for seven of eight comparison targets (TEM1 38% and LepB 100% at n=3 are exceptions; full per-target overlap in the S25 footnote)"),
    ("**Boltz-2 metrics:** `[BOLTZ_PENDING]`.",
     "**Boltz-2:** on the 166-row nested subset (**165 ok / 1 NA**; GyrA_Ec decoy failed after two attempts), pooled AUC was **0.718 (target-level 95% CI 0.590–0.840)**, per-target AUC 0.367–0.960, with a paired ΔAUC versus Vina of **+0.089 (95% CI −0.122 to +0.311)**; observed signed rho ranged −0.702 to +0.411 with no consistent direction, and median ligand_iptm 0.437–0.923. Top-5 recall/EF were not reported for n=20; Boltz-2 is exploratory and target-specific (GyrA_Ec 0.367 vs Vina 0.878) and does not support a leading-engine claim; **no single engine dominates across targets**."),
    ("**Figure 7 caption (draft)**", "**Figure 7 caption (final)**"),
    ("(E) Boltz-2 `[BOLTZ_PENDING]`.", "(E) Boltz-2 `ligand_iptm` per-target AUC (166-row nested subset; 165 ok/1 NA; pooled AUC 0.718, 95% CI 0.590–0.840)."),
    ("Source: `results/trackA/method_axis/analysis/fig7_data.csv` and `figures/fig7_method_axis_draft.png` / `.pdf` (draft).",
     "Source: `results/trackA/method_axis/analysis/fig7_data.csv`; figure file `figures/v3/fig7_method_axis.png` / `.pdf`."),
    ("## Pending list\n\n- `[BOLTZ_PENDING]` all Boltz-2 metrics: per-target AUC, pooled AUC/CI, rho, confidence distribution; no complete Boltz-2 output/score table yet (current aggregate is partial).\n- Fig 7 panel E and final 300 dpi figure not yet produced (panel D now uses single-chain data; S28 updated).\n- Replace dry-run DiffDock numbers only after Boltz final data; keep S25–S28 consistent with v3 references.\n- The method-axis draft does not modify any signed v3 file.",
     "## Pending list\n\n- No data pending: DiffDock main, single-chain and Boltz-2 are complete; Boltz n=165 ok/1 NA recorded in S27b.\n- Final 300 dpi Fig 7 and the v3 integration diff/apply are the remaining steps; S25–S28 must stay consistent with the v3 text.\n- The method-axis draft does not modify any signed v3 file."),
]
for old, new in repls:
    if old not in t:
        print("MISSING:", old[:160])
    else:
        t = t.replace(old, new)
p.write_text(t, encoding="utf-8")
print("draft updated", p.stat().st_size)
