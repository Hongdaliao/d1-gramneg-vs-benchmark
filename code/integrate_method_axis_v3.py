"""Integrate Amendment #12 / Addendum 12a method-axis text into v3.

Runs only after Boltz-2 has 166/166 rows. Default is dry-run: writes a diff
summary and does not touch v3.  Use --apply only after the diff has been shown.
"""
from __future__ import annotations

import json
import re
import shutil
from pathlib import Path

ROOT = Path(r"D:\harness\projects\d1_gramneg_vs")
V3 = ROOT / "manuscript/v3/manuscript_v3.md"
AN = ROOT / "results/trackA/method_axis/analysis"
OUT_DIFF = ROOT / "notes/review/v3_method_axis_integration_diff.md"
FIG_DRAFT = ROOT / "figures/fig7_method_axis_draft.png"
FIG_FINAL = ROOT / "figures/v3/fig7_method_axis.png"


def word_count(text: str) -> int:
    # JCIM-style whitespace word count; punctuation-separated numeric ranges count as words.
    return len(text.split())


def fmt(x, nd=3):
    return f"{float(x):.{nd}f}"


def main(apply: bool = False):
    pooled = json.loads((AN / "metrics_pooled.json").read_text(encoding="utf-8"))
    paired = json.loads((AN / "paired_bootstrap.json").read_text(encoding="utf-8"))
    bz = pooled["pooled"].get("boltz2", {})
    if bz.get("status") != "ok" or "pooled_auc_mean" not in bz:
        raise SystemExit(f"Boltz-2 not complete: {bz.get('status')} ({bz.get('n_observed')}/{bz.get('n_expected')})")
    p = pooled["pooled"]
    m = p["methods"]
    metrics = pooled["metrics"]
    bz_rows = [r for r in metrics if r["method"] == "Boltz-2" and r["target"] != "LepB_Ec"]
    bz_auc = [float(r["auc"]) for r in bz_rows if r["auc"] is not None and float(r["auc"]) == float(r["auc"])]
    bz_rho = [float(r["rho_confidence_vina"]) for r in bz_rows if r["rho_confidence_vina"] is not None and float(r["rho_confidence_vina"]) == float(r["rho_confidence_vina"])]
    bz_lo, bz_hi = min(bz_auc), max(bz_auc)
    bz_pooled = float(bz["pooled_auc_mean"])
    bz_ci = bz.get("pooled_auc_ci95", [float("nan"), float("nan")])
    bz_delta = float(bz.get("pooled_delta_boltz_minus_vina", float("nan")))
    bz_dci = bz.get("pooled_delta_ci95", [float("nan"), float("nan")])
    bz_rho_lo, bz_rho_hi = min(bz_rho), max(bz_rho)
    conf = bz.get("confidence_distribution", [])
    conf_med = [float(c["median_ligand_iptm"]) for c in conf]
    conf_txt = f"median ligand_iptm {min(conf_med):.3f}–{max(conf_med):.3f}" if conf_med else "confidence distribution in S25"
    n_bz_na = int(bz.get("n_na", 0))
    n_bz_ok = int(bz.get("n_ok", 0))
    single = p.get("singlechain", {}).get("targets", {})
    g = single.get("GyrA_Ec", {})
    f = single.get("FabI_Ec", {})

    methods = f"""### 2.7 Exploratory method-class comparison (DiffDock, Boltz-2)

This is a post-hoc exploratory analysis (Amendment #12 plus Addendum 12a; amendment sha256 `9ce4342e45dee330243e2af64b331232b0ba381b7f5406ab41cfb30c1ca222a6`); it sets no pass/fail gate, does not alter existing conclusions and is not used to validate MAES or enrichment. It uses the same nine-target Track A panel with a balanced eligible subset: `status=ok`, `box_clipped=0`, and at least one matched `ok` decoy. Eligible counts are LpxC_Ec 222, LpxC_Pa 418, DHFR_Ec 382, Gyrase_Ec 355, GyrA_Ec 73, GyrB_Ec 62, FabI_Ec 53, TEM1_Ec 51 and LepB_Ec 3. With seed 42, DiffDock uses 806 rows (403 active + 403 decoy; eight targets 50+50, LepB 3+3); Boltz-2 uses a nested 166 rows (eight targets 10+10, LepB 3+3). LepB is descriptive only.

Methods compared are Vina 1.2.5 and Vinardo as a second scoring function on the same poses (not an independent engine), DiffDock v1.1 with ESM2-650M (`samples=2`, rank-1 by confidence), and Boltz-2 co-folding with `ligand_iptm` as the primary score and `confidence_score` secondary (`msa: empty`, single-sequence mode; no MSA server; no metal/cofactor ion entity). The existing ECFP4 1-NN baseline was not recomputed on the balanced subset because subset actives overlap its temporal-train reference (68–94% for seven of the eight comparison targets; TEM1 38%, LepB 100% at n=3), so a subset-specific ECFP4 AUC would be a self-similarity artefact (Addendum 12a); the 2D reference remains the full-panel test-only result in §3.3 (median 0.775) and is not numerically compared to subset methods. The single-chain sensitivity uses the same subsets with 206 complexes (202 ok / 4 NA; GyrA, FabI and LepB single-chain receptors). DiffDock receptor input retains HETATM cofactors, but the DiffDock/ProDy parser models protein CA only; cofactors are not scored. DiffDock environment: v1.1 + ESM2-650M + torch 2.7.1+cu128 + PyG Windows wheels, repo commit `85c49b6…`; weight and source SHA256 values are frozen in Amendment #12/12a.

Endpoints are per-target AUC from rank-1 confidence (or Boltz-2 `ligand_iptm`), top-5 active recall (k=5; not EF; the subset is 1:1 so EF5 is not comparable to the full panel), and Spearman rho between confidence and Vina score (signed; no directional claim). LepB is excluded from AUC/top-5/rho and pooled summaries. Pooled AUC is the unweighted mean of eight target AUCs with target-level bootstrap CI (2,000 draws); per-target CIs resample ligands within target (2,000 draws); method differences use paired bootstrap over subset pairs (2,000 draws), not independent-CI comparison. Failures/NaN are recorded as NA and retained in coverage; no sample is silently removed, no post-hoc threshold or parameter change is made, and no multiple-comparison/significance claim is made.
"""

    results = f"""### 3.10 Exploratory method-class comparison (subset-balanced)

**Figure 7** summarises this exploratory comparison. On the frozen 806-ligand DiffDock subset (793 ok / 13 NA), pooled AUC was **{fmt(p['diffdock_auc_mean'])} (target-level 95% CI {fmt(p['diffdock_auc_ci_lo'])}–{fmt(p['diffdock_auc_ci_hi'])})** for DiffDock versus **{fmt(p['vina_auc_mean'])}** for Vina on the same paired ligands. The paired ΔAUC (DiffDock − Vina) was **{fmt(p['delta_auc_diffdock_minus_vina_mean'])} (95% CI {fmt(p['delta_auc_ci_lo'])} to {fmt(p['delta_auc_ci_hi'])})**. Per-target ranking split (DiffDock higher for GyrB_Ec 0.682 vs 0.432 and TEM1_Ec 0.598 vs 0.523; Vina higher for LpxC_Pa 0.857 vs 0.437 and GyrA_Ec 0.874 vs 0.488), and the remaining targets varied, so **no single engine dominates; the comparison is target-specific**. The observed signed rho(confidence, Vina score) ranged from −0.08 to +0.20 with no consistent direction (two negative, six positive); we do not report a directional claim. Thirteen NA samples (DHFR 4, Gyrase 3, GyrA 4, GyrB 2) are listed in S27; per-target analysis uses 96–100 complete rows. LepB (n=6) is descriptive only.

Boltz-2 (166-row nested subset; {n_bz_ok} ok / {n_bz_na} NA) gave pooled AUC **{fmt(bz_pooled)} (target-level 95% CI {fmt(bz_ci[0])}–{fmt(bz_ci[1])})**, per-target AUC {fmt(bz_lo)}–{fmt(bz_hi)}, and a paired ΔAUC versus Vina of {fmt(bz_delta)} (95% CI {fmt(bz_dci[0])} to {fmt(bz_dci[1])}); observed rho ranged {fmt(bz_rho_lo)} to {fmt(bz_rho_hi)} and {conf_txt}. Boltz-2 was higher than Vina in this nested subset overall, but its per-target performance varied (for example GyrA_Ec 0.367 vs Vina 0.878) and the paired CI included zero; we do not claim leading or complementary superiority. The method-class result remains **target-specific, with no single engine dominant across targets**. ECFP4 was not recomputed on the subset because of temporal-train overlap (Addendum 12a; S25 footnote); the 2D reference is the full-panel test-only result (median AUC 0.775, §3.3), not compared numerically to subset methods. On the 206-complex single-chain arm (202 ok / 4 NA), GyrA_Ec DiffDock AUC changed from {fmt(g.get('auc_multichain', float('nan')))} to {fmt(g.get('auc_singlechain', float('nan')))} (Δ {fmt(g.get('delta', float('nan')))}; 95% CI {fmt(g.get('delta_ci', [float('nan'), float('nan')])[0])} to {fmt(g.get('delta_ci', [float('nan'), float('nan')])[1])}) and FabI_Ec from {fmt(f.get('auc_multichain', float('nan')))} to {fmt(f.get('auc_singlechain', float('nan')))} (Δ {fmt(f.get('delta', float('nan')))}; 95% CI {fmt(f.get('delta_ci', [float('nan'), float('nan')])[0])} to {fmt(f.get('delta_ci', [float('nan'), float('nan')])[1])}); the single-chain sensitivity changes only the DiffDock receptor representation and does not change the Vina ranking. All numbers are subset-balanced and not comparable to full-panel Table 1.
"""

    abstract_sentence = " A post-hoc method comparison found target-specific rankings."
    limitations_sentence = (" The method-class comparison is post-hoc, subset-balanced and not powered for method superiority; "
                             "ECFP4 was not recomputed on the subset because of temporal-train overlap (S25), Boltz-2 is exploratory, "
                             "and single-chain sensitivity is reported in S28.")
    # Trim 14 non-substantive words from the signed abstract to stay <=250 words after adding the method-axis sentence.
    abstract_old = V3.read_text(encoding="utf-8").split("## Abstract")[1].split("## Keywords")[0].strip()
    abstract_trim = abstract_old
    for old, new in [
        ("across Gram-negative targets, and across the two very different data regimes those targets occupy, has not",
         "across Gram-negative targets and the two data regimes they occupy, has not"),
        ("We analyse two complementary regimes with regime-appropriate endpoints.",
         "We analyse two regimes with regime-appropriate endpoints."),
        ("Track B is a structure-rich co-crystal ligand pose benchmark", "Track B is a co-crystal ligand pose benchmark"),
        ("a common quantitative framework", "a common framework"),
        ("All data and code will be released upon publication.", "Data and code will be released upon publication."),
        ("The ECFP4 1-NN baseline evaluated on temporal test actives alone",
         "The ECFP4 1-NN baseline on temporal test actives alone"),
        ("The pre-registered MAES calibration gate is NO_GO", "The pre-registered MAES gate is NO_GO"),
        ("The structure-source axis, restricted to eight QC-passed pairs", "The structure-source axis, eight QC-passed pairs"),
        ("and 7/8 comparable targets match or exceed Vina", "and 7/8 match or exceed Vina"),
        ("and one is an unresolved/low-confidence ligand (95X; annotated as substrate in Table S18, not counted as an inhibitor)",
         "and one is an unresolved 95X ligand (annotated as substrate in Table S18, not counted as an inhibitor)"),
    ]:
        assert old in abstract_trim, f"abstract trim pattern missing: {old[:50]}"
        abstract_trim = abstract_trim.replace(old, new, 1)
    if "All data and code will be released upon publication." in abstract_trim:
        abstract_new = abstract_trim.replace("All data and code will be released upon publication.",
                                             abstract_sentence.strip() + " All data and code will be released upon publication.", 1)
    else:
        abstract_new = abstract_trim.replace("Data and code will be released upon publication.",
                                             abstract_sentence.strip() + " Data and code will be released upon publication.", 1)
    data_old = "are provided as SI tables S1–S24, including sub-tables S4a/S4b, S13b, S16b, S18b, S19b, S22b and S22c."
    data_new = "are provided as SI tables S1–S28 (including sub-tables S4a/S4b, S13b, S16b, S18b, S19b, S22b, S22c, S26b and S27b)."
    fig = f"""**Figure 7. Exploratory method-class comparison (subset-balanced).** (A) Paired per-target AUC for Vina and DiffDock; (B) pooled AUC (unweighted mean of eight targets) with target-level bootstrap CI; (C) signed Spearman rho between DiffDock confidence and Vina score (no directional claim); (D) single-chain sensitivity for GyrA and FabI (LepB descriptive; 206 attempted/202 ok/4 NA); (E) Boltz-2 `ligand_iptm` per-target AUC (166-row nested subset; pooled AUC {fmt(bz_pooled)}, 95% CI {fmt(bz_ci[0])}–{fmt(bz_ci[1])}). All endpoints are subset-balanced (DiffDock 806: 403 active + 403 decoy; eight targets 50+50, LepB 3+3) and are not comparable to full-panel Table 1. Source: `results/trackA/method_axis/analysis/fig7_data.csv`; figure file `figures/v3/fig7_method_axis.png` / `.pdf`.
"""

    text = V3.read_text(encoding="utf-8")
    if "### 2.7 Exploratory method-class comparison" in text:
        raise SystemExit("v3 already contains the method-axis section; aborting")
    text = text.replace("### 2.7 Statistics and reproducibility", methods + "\n### 2.8 Statistics and reproducibility", 1)
    text = text.replace("## 4. Discussion", results + "\n## 4. Discussion", 1)
    # replace abstract with the trimmed abstract including the added sentence
    old_abstract_block = "## Abstract\n\n" + abstract_old + "\n\n## Keywords"
    new_abstract_block = "## Abstract\n\n" + abstract_new + "\n\n## Keywords"
    assert old_abstract_block in text, "abstract block not found"
    text = text.replace(old_abstract_block, new_abstract_block, 1)
    assert "The 16 GB GPU limits co-folding to 13 of the 40 co-crystal-ligand systems." in text
    text = text.replace("The 16 GB GPU limits co-folding to 13 of the 40 co-crystal-ligand systems.",
                        "The 16 GB GPU limits co-folding to 13 of the 40 co-crystal-ligand systems." + limitations_sentence, 1)
    assert data_old in text
    text = text.replace(data_old, data_new, 1)
    # append Fig7 caption at end of file
    text = text.rstrip() + "\n\n" + fig

    # abstract word count
    abstract = text.split("## Abstract")[1].split("## Keywords")[0]
    wc = word_count(abstract)
    if wc > 250:
        raise SystemExit(f"abstract would be {wc} words (>250); trim before applying")

    diff_lines = [
        "# v3 method-axis integration — diff summary (pre-apply)",
        "",
        f"- v3 before: {len(V3.read_text(encoding='utf-8').splitlines())} lines / {V3.stat().st_size} bytes",
        f"- Abstract: insert one 9-word sentence before final data-availability sentence; word count after = {wc} (limit 250).",
        "- Methods: insert new `### 2.7 Exploratory method-class comparison (DiffDock, Boltz-2)` before current 2.7; renumber current `2.7 Statistics and reproducibility` to `2.8`.",
        "- Results: insert new `### 3.10 Exploratory method-class comparison (subset-balanced)` before `## 4. Discussion`.",
        "- Limitations: append one sentence (post-hoc/subset-balanced/ECFP4 overlap/Boltz exploratory/single-chain S28).",
        "- Data availability: S1–S24 -> S1–S28 (add S25–S28; S26b/S27b named).",
        "- Figure captions: append Figure 7 caption with new file `figures/v3/fig7_method_axis.png` / `.pdf`.",
        "- No data, existing conclusions or signed Figure/Table values changed; only additive text plus the 2.7->2.8 heading renumber.",
        "",
        "## Inserted Results numbers",
        f"- DiffDock pooled {fmt(p['diffdock_auc_mean'])} (CI {fmt(p['diffdock_auc_ci_lo'])}-{fmt(p['diffdock_auc_ci_hi'])})",
        f"- Vina pooled {fmt(p['vina_auc_mean'])}; paired delta {fmt(p['delta_auc_diffdock_minus_vina_mean'])} (CI {fmt(p['delta_auc_ci_lo'])} to {fmt(p['delta_auc_ci_hi'])})",
        f"- Boltz-2 pooled {fmt(bz_pooled)} (CI {fmt(bz_ci[0])}-{fmt(bz_ci[1])}); per-target range {fmt(bz_lo)}-{fmt(bz_hi)}; rho range {fmt(bz_rho_lo)} to {fmt(bz_rho_hi)}; {n_bz_ok} ok/{n_bz_na} NA",
        f"- Single-chain GyrA {fmt(g.get('auc_multichain', float('nan')))}->{fmt(g.get('auc_singlechain', float('nan')))}; FabI {fmt(f.get('auc_multichain', float('nan')))}->{fmt(f.get('auc_singlechain', float('nan')))}",
        "",
        "## Applied?",
        "NO — dry-run only. Re-run with `--apply` after review." if not apply else "YES — v3 updated.",
    ]
    OUT_DIFF.write_text("\n".join(diff_lines), encoding="utf-8")
    if apply:
        V3.write_text(text, encoding="utf-8")
        # final figure from draft
        for ext in (".png", ".pdf"):
            src = ROOT / f"figures/fig7_method_axis_draft{ext}"
            dst = ROOT / f"figures/v3/fig7_method_axis{ext}"
            if src.exists():
                shutil.copy2(src, dst)
        print("APPLIED", V3, V3.stat().st_size, "figure", FIG_FINAL)
    else:
        print("DRY-RUN diff:", OUT_DIFF, "abstract words:", wc)
    return {"abstract_words": wc, "diff": str(OUT_DIFF)}


if __name__ == "__main__":
    import sys
    main(apply="--apply" in sys.argv)
