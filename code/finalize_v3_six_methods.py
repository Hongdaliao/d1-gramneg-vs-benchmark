"""Final v3 integration: Amendment #12/#13/#14 five-method-class axis.

Run with --apply after reviewing the diff summary.  Only this script edits v3.
"""
from __future__ import annotations
import json, re, shutil
from pathlib import Path

ROOT = Path(r"D:\harness\projects\d1_gramneg_vs")
V3 = ROOT / "manuscript/v3/manuscript_v3.md"
AN = ROOT / "results/trackA/method_axis/analysis"
DIFF = ROOT / "notes/review/v3_six_method_integration_diff.md"


def f(x, nd=3):
    return f"{float(x):.{nd}f}"


def main(apply=False):
    p = json.loads((AN / "metrics_pooled_final.json").read_text(encoding="utf-8"))
    import pandas as pd
    pp = {r["method"]: r for r in json.loads((AN / "paired_pooled_final.json").read_text(encoding="utf-8"))}
    metrics = pd.read_csv(AN / "metrics_all_methods_final.csv")
    cov = pd.read_csv(AN / "coverage_all_methods_final.csv")
    old12 = (ROOT / "notes/preregistration_amendment12_method_axis_draft.md").read_bytes()
    old13 = (ROOT / "notes/preregistration_amendment13_drugclip_draft.md").read_bytes()
    old14 = (ROOT / "notes/preregistration_amendment14_unidock_draft.md").read_bytes()
    import hashlib
    sha = lambda b: hashlib.sha256(b).hexdigest()
    bz = p["methods"]["boltz2"]["pooled"]; dd = p["methods"]["diffdock_main"]["pooled"]
    vina = p["methods"]["vina"]["pooled"]; dclip = p["methods"]["drugclip"]["pooled"]; uni = p["methods"]["unidock"]["pooled"]
    common = {k: p["methods"][k]["common7"]["mean_auc"] for k in p["methods"]}
    def rng(method, key="auc"):
        d = metrics[metrics["method"] == method]
        d = d[d["target"].isin(["LpxC_Ec","LpxC_Pa","DHFR_Ec","Gyrase_Ec","GyrA_Ec","GyrB_Ec","FabI_Ec","TEM1_Ec"])]
        return float(d["auc"].min()), float(d["auc"].max())
    dd_r, bz_r, dc_r, uni_r = rng("diffdock_main"), rng("boltz2"), rng("drugclip"), rng("unidock")
    # single-chain values from earlier frozen analysis
    old = json.loads((AN / "metrics_pooled.json").read_text(encoding="utf-8"))["pooled"]
    g = old["singlechain"]["targets"]["GyrA_Ec"]; fab = old["singlechain"]["targets"]["FabI_Ec"]
    na_main = len(pd.read_csv(AN / "coverage_na.csv").query("status != 'ok'")) if (AN / "coverage_na.csv").exists() else 13

    methods = f"""### 2.7 Exploratory method-class comparison (Vina, Uni-Dock, DiffDock, Boltz-2, DrugCLIP, ECFP4 reference)

This is a post-hoc exploratory analysis (Amendment #12, #13 and #14; sha256 `{sha(old12)}`, `{sha(old13)}` and `{sha(old14)}`; Addendum 12a and 14a). It sets no pass/fail gate, does not alter existing conclusions and is not used to validate MAES or enrichment. It uses the same nine-target Track A panel with a balanced eligible subset (`status=ok`, `box_clipped=0`, and at least one matched `ok` decoy): DiffDock main 806 rows (403 active + 403 decoy; eight targets 50+50, LepB 3+3), Boltz-2 166 rows (eight targets 10+10, LepB 3+3), and Uni-Dock/DrugCLIP the same 806-ligand subset; seed 42. LepB is descriptive only.

Methods compared are Vina 1.2.5 and Vinardo as a second scoring function on the same poses (not an independent engine); **Uni-Dock v1.2.0 [19]** in WSL with the same receptor PDBQT, box, `--scoring vina --exhaustiveness 4 --seed 42 --num_modes 3` as Vina; **DiffDock v1.1 [21]** with ESM2-650M (`samples=2`, rank-1 by confidence); **Boltz-2** co-folding with `ligand_iptm` as the primary score and `confidence_score` secondary (`msa: empty`; no MSA server; no metal/cofactor ion entity); and **DrugCLIP [20]** with the official 6_folds model (CC BY-NC 4.0; pretrained weights; train-test overlap cannot be excluded; native-ligand 6 Å pocket encoding excludes metals/cofactors). ECFP4 is not a subset-scored method: its test-only reference was **not** recomputed on the balanced subset because subset actives overlap the temporal-train reference (Addendum 12a; 68–94% for seven of eight comparison targets, TEM1 38%, LepB 100% at n=3); the full-panel test-only result in §3.3 (median 0.775) is the 2D reference and is not numerically compared to subset methods. The single-chain sensitivity uses the same DiffDock subsets with 206 complexes (202 ok / 4 NA; GyrA, FabI and LepB). DiffDock receptor input retains HETATM cofactors, but the DiffDock/ProDy parser models protein CA only. DiffDock environment: v1.1 + ESM2-650M + torch 2.7.1+cu128 + PyG Windows wheels, repo commit `85c49b6…`; weight and source SHA256 values are frozen in Amendment #12/13/14.

Endpoints are per-target AUC from rank-1 confidence or `ligand_iptm`, top-5 active recall (k=5; not EF; the subset is 1:1 so EF5 is not comparable to the full panel), and Spearman rho between the method score and Vina score (signed; no directional claim). LepB is excluded from AUC/top-5/rho and pooled summaries. Pooled AUC is the unweighted mean of target AUCs with target-level bootstrap CI (2,000 draws); the main 8-target pooled set covers Vina, DiffDock, Boltz-2 and DrugCLIP; **Uni-Dock is pooled over seven targets (FabI excluded because 13/50 FabI actives are NA_positive_affinity and 1 has no model)**, and a common seven-target sensitivity set is reported for all methods. Per-target CIs resample ligands within target (2,000 draws); method differences use paired bootstrap over subset pairs, not independent-CI comparison. Failures/NaN are recorded as NA and retained in coverage; positive Uni-Dock affinities are treated as NA; TEM1_Ec boron actives use the pre-specified `_B2C` surrogate (22 complexes, Addendum 14a). No post-hoc threshold or parameter change is made, and no multiple-comparison/significance claim is made.
"""
    results = f"""### 3.10 Exploratory method-class comparison (subset-balanced)

**Figure 7** summarises this exploratory comparison. On the frozen 806-ligand DiffDock subset (793 ok / 13 NA), pooled AUC was **{f(dd['mean_auc'])} (target-level 95% CI {f(dd['ci95'][0])}–{f(dd['ci95'][1])})** for DiffDock versus **{f(vina['mean_auc'])}** for Vina on the same paired ligands; the paired ΔAUC (DiffDock − Vina) was **{f(pp['diffdock_main']['mean_delta_vs_vina'])} (95% CI {f(pp['diffdock_main']['ci_lo'])} to {f(pp['diffdock_main']['ci_hi'])})**. Per-target ranking varied (DiffDock higher for GyrB_Ec 0.682 vs 0.432 and TEM1_Ec 0.598 vs 0.523; Vina higher for LpxC_Pa 0.857 vs 0.437 and GyrA_Ec 0.874 vs 0.488), so no single engine dominates; the comparison is target-specific. Observed signed rho(confidence, Vina score) ranged from −0.08 to +0.20 with no consistent direction; 13 NA samples are listed in S27.

Boltz-2 (166-row nested subset; 165 ok / 1 NA) gave pooled AUC **{f(bz['mean_auc'])} (95% CI {f(bz['ci95'][0])}–{f(bz['ci95'][1])})**, per-target AUC {f(bz_r[0])}–{f(bz_r[1])}, and paired ΔAUC versus Vina of **{f(pp['boltz2']['mean_delta_vs_vina'])} (95% CI {f(pp['boltz2']['ci_lo'])} to {f(pp['boltz2']['ci_hi'])})**. Boltz-2 was higher than Vina in this nested subset overall, but its per-target performance varied (GyrA_Ec 0.367 vs Vina 0.878) and the paired CI included zero; it is exploratory and does not support a leading or complementary-superiority claim.

DrugCLIP (801 ok / 5 NA; all five NA are GyrA decoys with RDKit 3D generation failure) gave pooled 8-target AUC **{f(dclip['mean_auc'])} (95% CI {f(dclip['ci95'][0])}–{f(dclip['ci95'][1])})**, per-target AUC {f(dc_r[0])}–{f(dc_r[1])}, and paired ΔAUC versus Vina of **{f(pp['drugclip']['mean_delta_vs_vina'])} (95% CI {f(pp['drugclip']['ci_lo'])} to {f(pp['drugclip']['ci_hi'])})**. DrugCLIP uses pretrained weights (CC BY-NC 4.0), so train-test overlap cannot be excluded; the apparent performance is therefore descriptive only.

Uni-Dock (775 ok / 31 NA; 30 positive-affinity and one no-model; FabI 86 ok / 14 NA) gave a **7-target pooled AUC of {f(uni['mean_auc'])} (95% CI {f(uni['ci95'][0])}–{f(uni['ci95'][1])})**, per-target AUC {f(uni_r[0])}–{f(uni_r[1])}, and paired ΔAUC versus Vina of **{f(pp['unidock']['mean_delta_vs_vina'])} (95% CI {f(pp['unidock']['ci_lo'])} to {f(pp['unidock']['ci_hi'])})**. TEM1_Ec was fully covered after the 22-complex `_B2C` boron surrogate; FabI is excluded from the pooled Uni-Dock estimate because of selective positive-affinity missingness. The 7-target common sensitivity gave means of {f(common['vina'])} (Vina), {f(common['diffdock_main'])} (DiffDock), {f(common['boltz2'])} (Boltz-2), {f(common['drugclip'])} (DrugCLIP) and {f(common['unidock'])} (Uni-Dock).

ECFP4 was not recomputed on the subset because of temporal-train overlap (Addendum 12a; S25 footnote); the 2D reference is the full-panel test-only result (median AUC 0.775, §3.3), not compared numerically to subset methods. On the 206-complex single-chain arm (202 ok / 4 NA), GyrA_Ec DiffDock AUC changed from {f(g.get('auc_multichain'))} to {f(g.get('auc_singlechain'))} (Δ {f(g.get('delta'))}; 95% CI {f(g['delta_ci'][0])} to {f(g['delta_ci'][1])}) and FabI_Ec from {f(fab.get('auc_multichain'))} to {f(fab.get('auc_singlechain'))} (Δ {f(fab.get('delta'))}; 95% CI {f(fab['delta_ci'][0])} to {f(fab['delta_ci'][1])}); the sensitivity changes only the DiffDock receptor representation and does not change the Vina ranking. All numbers are subset-balanced and not comparable to full-panel Table 1.
"""
    abstract_sentence = "Five method classes gave target-specific rankings on balanced subsets."
    limitations_sentence = (" The method-class comparison is post-hoc, subset-balanced and not powered for method superiority; "
                             "Uni-Dock has positive-affinity/FabI missingness and needs the TEM1 boron surrogate, DrugCLIP is a pretrained model whose train-test overlap cannot be excluded, ECFP4 was not recomputed on the subset (S25), Boltz-2 is exploratory, and single-chain sensitivity is reported in S28.")
    data_old = "are provided as SI tables S1–S28 (including sub-tables S4a/S4b, S13b, S16b, S18b, S19b, S22b, S22c, S26b and S27b)."
    data_new = "are provided as SI tables S1–S31 (including sub-tables S4a/S4b, S13b, S16b, S18b, S19b, S22b, S22c, S26b, S27b, S31b, S31c and S31d)."
    fig = f"""**Figure 7. Exploratory method-class comparison (subset-balanced).** (A) Per-target AUC for Vina, Uni-Dock, DiffDock, Boltz-2 and DrugCLIP; (B) pooled AUC with target-level bootstrap CI (main 8 targets; Uni-Dock 7 targets); (C) paired per-target ΔAUC versus Vina; (D) DiffDock single-chain sensitivity for GyrA and FabI (LepB descriptive; 206 attempted/202 ok/4 NA); (E) common seven-target sensitivity means; (F) method coverage (ok rows). Boltz-2 and DrugCLIP are exploratory/pretrained and no leading or complementary-superiority claim is made. All endpoints are subset-balanced and not comparable to full-panel Table 1. Source: `results/trackA/method_axis/analysis/fig7_data.csv`; figure file `figures/v3/fig7_method_axis.png` / `.pdf`.
"""
    text = V3.read_text(encoding="utf-8")
    if "### 2.7 Exploratory method-class comparison (Vina, Uni-Dock" in text:
        raise SystemExit("v3 already final six-method version")
    # old line 47 outdated sentence
    old_sent = "Uni-Dock [19] was considered but could not be deployed on the available Windows host; DrugCLIP [20] was available as a contrastive-screening arm but was not run because of environment/compute constraints."
    new_sent = "Uni-Dock [19] was run in WSL (v1.2.0) and DrugCLIP [20] was run with the 6_folds model on the same 806-ligand subset; execution, pooling and coverage are described in §2.7 and S29/S30."
    assert old_sent in text; text = text.replace(old_sent, new_sent, 1)
    # replace Methods 2.7 block
    start = text.index("### 2.7 Exploratory method-class comparison")
    end = text.index("### 2.8 Statistics and reproducibility")
    text = text[:start] + methods + "\n" + text[end:]
    # replace Results 3.10 block
    start = text.index("### 3.10 Exploratory method-class comparison")
    end = text.index("## 4. Discussion")
    text = text[:start] + results + "\n" + text[end:]
    # Abstract sentence
    old_abs = "A post-hoc method comparison found target-specific rankings."
    assert old_abs in text
    text = text.replace(old_abs, abstract_sentence, 1)
    # Limitations sentence: replace old method-axis sentence
    old_lim = (" The method-class comparison is post-hoc, subset-balanced and not powered for method superiority; "
               "ECFP4 was not recomputed on the subset because of temporal-train overlap (S25), Boltz-2 is exploratory, "
               "and single-chain sensitivity is reported in S28.")
    assert old_lim in text
    text = text.replace(old_lim, limitations_sentence, 1)
    assert data_old in text
    text = text.replace(data_old, data_new, 1)
    # replace Fig7 caption (last occurrence)
    idx = text.rfind("**Figure 7. Exploratory method-class comparison")
    if idx != -1:
        text = text[:idx] + fig
    else:
        text = text.rstrip() + "\n\n" + fig
    # --- add DiffDock reference and renumber citations >=21 ---
    marker = "@@DIFFDOCKREF@@"
    # new methods text already has [21] literal; change it to marker before shifting
    text = text.replace("DiffDock v1.1 [21]", f"DiffDock v1.1 {marker}", 1)

    def shift_cites(s: str) -> str:
        def repl(m):
            a, b = m.group(1), m.group(2)
            if b is None:
                n = int(a)
                return f"[{n+1}]" if n >= 21 else m.group(0)
            aa, bb = int(a), int(b)
            if aa >= 21:
                return f"[{aa+1}-{bb+1}]"
            if bb >= 21:
                return f"[{aa}-20,22-{bb+1}]"
            return m.group(0)
        return re.sub(r"\[(\d+)(?:-(\d+))?\]", repl, s)

    body_start = 0
    refs_idx = text.index("## References")
    body = shift_cites(text[:refs_idx])
    refs_end = text.index("## Figure captions")
    refs = text[refs_idx:refs_end]
    # shift reference entries 21-38
    def shift_ref(m):
        n = int(m.group(1))
        return f"[{n+1}]" if n >= 21 else m.group(0)
    refs = re.sub(r"\[(\d+)\]", shift_ref, refs)
    new_ref = "[21] Corso G, Stärk H, Jing B, et al. DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking. International Conference on Learning Representations (ICLR), 2023. arXiv:2210.01776\n"
    # insert after old [20] line
    lines = refs.splitlines()
    out_lines, inserted = [], False
    for line in lines:
        out_lines.append(line)
        if line.startswith("[20]") and not inserted:
            out_lines.append(new_ref.rstrip()); inserted = True
    assert inserted
    refs = "\n".join(out_lines) + "\n\n"
    text = body + refs + text[refs_end:]
    text = text.replace(marker, "[21]", 1)
    # word count
    wc = len(text.split("## Abstract")[1].split("## Keywords")[0].split())
    if wc > 250:
        raise SystemExit(f"abstract {wc} > 250")
    lines = text.splitlines()
    diff = f"""# v3 six-method integration diff{' (APPLIED)' if apply else ' (DRY-RUN)'}

- before: 309 lines / 57577 bytes / sha 2c47fff6…
- after: {len(lines)} lines / {len(text.encode())} bytes
- Abstract: {wc} words (limit 250); replaced 7-word sentence with 9-word five-method sentence.
- Methods 2.7: replaced with five-method-class text (Vina/Vinardo, Uni-Dock, DiffDock, Boltz-2, DrugCLIP, ECFP4 reference); line-47 stale Uni-Dock/DrugCLIP sentence replaced.
- Results 3.10: replaced with six-method numbers; Boltz/DiffDock retained; new DrugCLIP/Uni-Dock/common-7.
- Limitations/Data availability S1–S31; Figure 7 A–F caption; DiffDock reference inserted as [21], old [21–38] shifted to [22–39].
- fixes included: E DiffDock ref, 'a audit' (already fixed), Boltz CI 0.585–0.840.
- No existing data/values changed beyond method-axis text and reference numbering.
"""
    DIFF.write_text(diff, encoding="utf-8")
    if apply:
        V3.write_text(text, encoding="utf-8")
        for ext in (".png", ".pdf"):
            s = ROOT / f"figures/v3/fig7_method_axis{ext}"
            if s.exists(): pass
        print("APPLIED", V3, V3.stat().st_size, "abstract", wc)
    else:
        print("DRY-RUN", DIFF, "abstract", wc)
    return wc


if __name__ == "__main__":
    import sys
    main(apply="--apply" in sys.argv)
