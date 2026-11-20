import json
from pathlib import Path
import pandas as pd

OUT = Path(r"results/trackA/method_axis/analysis")
pooled = json.loads((OUT / "metrics_pooled.json").read_text(encoding="utf-8"))
p = pooled["pooled"]
bz = p["boltz2"]
metrics = pd.read_csv(OUT / "metrics_per_target.csv")
cov = pd.read_csv(OUT / "coverage_na.csv")
bna = pd.read_csv(OUT / "boltz2_coverage_na.csv")
single = pd.read_csv(OUT / "singlechain_comparison.csv")
scna = pd.read_csv(OUT / "singlechain_coverage_na.csv")
overlap = pd.read_csv(OUT / "ecfp4_leakage_overlap.csv")

brows = [r for r in bz["metrics"] if r["target"] != "LepB_Ec" and r["auc"] == r["auc"]]
brho = [r["rho_confidence_vina"] for r in brows if r["rho_confidence_vina"] == r["rho_confidence_vina"]]
conf = bz.get("confidence_distribution", [])
cmed = [c["median_ligand_iptm"] for c in conf]

lines = [
    "# Amendment #12 method-axis analysis — final status (before E recheck)",
    "",
    "> DiffDock main + single-chain + Boltz-2 complete; ECFP4 subset reference withdrawn under Addendum 12a (leakage). No GPU used by this analysis.",
    "",
    "## Provenance",
    f"- Amendment (including 12a) sha256: `{pooled['metadata']['amendment_sha256']}`",
    f"- DiffDock subset sha256: `{pooled['metadata']['subset_diffdock_sha256']}`",
    f"- Boltz subset sha256: `{pooled['metadata']['subset_boltz2_sha256']}`",
    f"- DiffDock main: 806 expected / {pooled['metadata']['n_diffdock_ok']} ok / {pooled['metadata']['n_diffdock_na']} NA",
    f"- Boltz-2: {bz['n_ok']} ok / {bz['n_na']} NA of {bz['n_observed']} rows",
    "",
    "## Main DiffDock/Vina (8 targets; LepB descriptive)",
    f"- DiffDock pooled AUC: **{p['diffdock_auc_mean']:.3f}** (target-level 95% CI {p['diffdock_auc_ci_lo']:.3f}–{p['diffdock_auc_ci_hi']:.3f})",
    f"- Vina pooled AUC (same paired ligands): **{p['vina_auc_mean']:.3f}**",
    f"- Paired ΔAUC DiffDock − Vina: **{p['delta_auc_diffdock_minus_vina_mean']:.3f}** (95% CI {p['delta_auc_ci_lo']:.3f} to {p['delta_auc_ci_hi']:.3f})",
    "",
    "## Boltz-2 (nested 166; ligand_iptm)",
    f"- Pooled AUC: **{bz['pooled_auc_mean']:.3f}** (target-level 95% CI {bz['pooled_auc_ci95'][0]:.3f}–{bz['pooled_auc_ci95'][1]:.3f}); per-target AUC {min(r['auc'] for r in brows):.3f}–{max(r['auc'] for r in brows):.3f}",
    f"- Paired ΔAUC vs Vina: **{bz['pooled_delta_boltz_minus_vina']:.3f}** (95% CI {bz['pooled_delta_ci95'][0]:.3f} to {bz['pooled_delta_ci95'][1]:.3f}); CI includes zero",
    f"- rho range {min(brho):.3f} to {max(brho):.3f} (no consistent direction); median ligand_iptm range {min(cmed):.3f}–{max(cmed):.3f}",
    "- Boltz-2 is higher overall in this nested subset but per-target varies (GyrA_Ec 0.367 vs Vina 0.878); no leading/complementary superiority claim.",
    "- Top-5 recall/EF not reported for Boltz n=20; q1 NA listed in S27b.",
    "",
    "## ECFP4 subset reference — WITHDRAWN (Addendum 12a)",
    "ECFP4 test-only was not recomputed on the balanced subset because subset actives overlap the temporal-train reference (68–94% for seven of eight comparison targets; TEM1 38%, LepB 100% n=3). The full-panel test-only result (median AUC 0.775, §3.3) remains the 2D reference and is not numerically compared to subset methods.",
    "",
    "## Single-chain sensitivity (206 attempted / 202 ok / 4 NA)",
]
for _, r in single.iterrows():
    lines.append(f"- {r['target']}: multichain {r['auc_multichain']:.3f} vs single-chain {r['auc_singlechain']:.3f}; Δ {r['delta_singlechain_minus_multichain']:+.3f} (95% CI {r['delta_ci_lo']:.3f} to {r['delta_ci_hi']:.3f}); n={int(r['n_analysis'])}")
lines += [
    "- Single-chain sensitivity changes only DiffDock receptor representation; Vina ranking is unchanged.",
    "",
    "## NA / coverage",
    f"- Main DiffDock NA: **{len(cov[cov['status'] != 'ok'])}** rows (DHFR 4, Gyrase 3, GyrA 4, GyrB 2); all in S27.",
    f"- Boltz-2 NA: **{len(bna[bna['status'] != 'ok'])}** row (GyrA_Ec dec_KOIREGUAEMHXSP-PYQAXFSVSA-N; attempts=2); S27b.",
    f"- Single-chain NA: **{len(scna)}** rows (GyrA); S28.",
    "",
    "## Output files",
    "- `metrics_per_target.csv`, `metrics_pooled.json`, `paired_bootstrap.json`",
    "- `coverage_na.csv`, `boltz2_coverage_na.csv`, `singlechain_coverage_na.csv`",
    "- `fig7_data.csv`, `table_S25...`, `table_S26...`, `table_S26b_boltz2_paired.csv`, `table_S27...`, `table_S27b_boltz2_coverage.csv`, `table_S28...`",
    "- `figures/fig7_method_axis_draft.png` / `.pdf`; final `figures/v3/fig7_method_axis.png` / `.pdf`",
    "- `manuscript/v3/method_axis_integration_draft.md`",
]
(OUT / "DRYRUN_ISSUES.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("wrote", OUT / "DRYRUN_ISSUES.md")
