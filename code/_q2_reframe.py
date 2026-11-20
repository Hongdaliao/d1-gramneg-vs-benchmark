from pathlib import Path

p = Path(r"manuscript/main_skeleton.md")
t = p.read_text(encoding="utf-8")

repls = [
    (
        "*Two data regimes, multi-axis transferability: a membrane-aware benchmark of structure-based virtual screening against Gram-negative antibacterial targets*",
        "*Benchmarking structure-based virtual screening across Gram-negative antibacterial targets: target-specific enrichment, transfer axes and failure modes*",
    ),
    (
        "For membrane poses we develop **MAES**, a calibrated Membrane-Aware Evaluation Score that combines chemical validity (building on PoseBusters), membrane-environment consistency and conformational-state consistency, with method ability and target difficulty separated by item-response modelling; metal-dependent targets are treated as a pre-specified stratum rather than a claimed new finding.",
        "We evaluate a membrane-aware pose-quality score (MAES) as a pre-registered exploratory arm, but its calibration gate fails (AUROC 0.496, |rho| 0.075; 146 poses), so we report it descriptively only and instead centre the paper on target-specific enrichment, three transfer axes and three benchmark/tool failure modes (metal coordination, boron atom typing, DUD-E-style decoy bias), including the negative external-consistency and calibration results; metal-dependent targets are treated as a pre-specified stratum rather than a claimed new finding.",
    ),
    (
        "4. Contributions (framed as integration, not primacy): (i) a two-regime panel for Gram-negative antibacterial targets; (ii) joint quantification of three transfer axes on that panel (applying time-split/predicted-receptor protocols established elsewhere [REF-TIME-1, REF-PRED-1]); (iii) MAES, a membrane-aware pose-quality score building on pose-validity and confidence-calibration work [REF-PB, REF-CONF-1], reported **descriptively only** after the W12 calibration gate returned **NO-GO** (AUROC 0.496, |rho| 0.075; see §2.6); (iv) decision rules with explicit uncertainty + open release (upon publication).",
        "4. Contributions (framed as integration, not primacy): (i) a two-regime panel for Gram-negative antibacterial targets (**9-target Track A, 1,748 actives / 28,136 successfully docked ligands**; **Track B: 40 inhibitor cognate-redocking poses across 9 membrane targets**); (ii) per-target quantification of three transfer axes on that panel (temporal, cross-species, structure-source; R6), showing **axis- and target-dependent** transfer rather than a single transferability score; (iii) three reproducible benchmark/tool failure modes — metal coordination, boron atom typing, DUD-E-style decoy bias — separated from target biology; (iv) a negative-result-driven reporting protocol (pre-registered gates, per-target cluster CIs, 2D baselines, no post-hoc threshold/method/data changes), with MAES retained only as an exploratory/descriptive pose-quality analysis.",
    ),
    (
        "### 2.6 MAES: a membrane-aware pose-quality score (Fig 6)",
        "### 2.6 MAES as an exploratory pose-quality score: calibration gate failed (Fig 6)",
    ),
    (
        "- Definition; calibration curves (Brier/ECE); AUROC vs baselines; ablation; IRT ability/difficulty map.\n- **W12 calibration gate (formal; 2026-11-19): NO-GO.** 146 pose observations (vina 40 / vinardo 40 / vina_afdb 40 / boltz2 26; 146 unique pose keys, 40 systems; all `phi_source=reference_structure`): primary **V+M** AUROC **0.496** (wild cluster CI 0.318–0.672), |rho(MAES,RMSD)| **0.075** (cluster CI −0.302–0.161), Brier 0.145, ECE 0.124; both pre-registered criteria fail (AUROC≥0.65 and |rho|≥0.40). Full V+M+S is exploratory only. **MAES is therefore descriptive only, not used for ranking or as a calibrated probability; no post-hoc threshold/method change.** Q2 fallback: the two-regime benchmark + three-axis transfer (R6) + failure-mode analyses carry the paper.",
        "- MAES is reported as an exploratory/descriptive pose-quality analysis only: the pre-registered W12 calibration gate fails on 146 pose observations (AUROC 0.496, wild-CI 0.318–0.672; |rho(MAES,RMSD)| 0.075, cluster-CI −0.302–0.161; Brier 0.145, ECE 0.124; 19 eval positives), so we make no probability, ranking or MAES-as-contribution claim, and present the calibration curves as a negative-result diagnostic (full V+M+S exploratory only).",
    ),
    (
        "| Fig 6 | MAES calibration + IRT | WP4 | **descriptive only; W12 formal NO-GO (AUROC 0.496, |rho| 0.075)** |",
        "| Fig 6 | **Calibration-gate negative result for MAES (reliability/Brier/ECE; AUROC 0.496; |rho| 0.075; 146 poses)** — exploratory/descriptive only, not a probability score or ranking figure | WP4 | **NO_GO; retained as negative result** |",
    ),
    (
        '| "MAES is calibrated and reusable" | AUROC/Brier/ECE + ablation + cross-target validation | **REJECTED — W12 formal NO-GO (AUROC 0.496, |rho| 0.075); descriptive only, Q2 fallback** |',
        '| "MAES provides a probability/ranking score" | AUROC/Brier/ECE + ablation + cross-target validation | **REJECTED — W12 formal NO-GO (AUROC 0.496, |rho| 0.075); descriptive only, Q2 fallback** |',
    ),
]
for old, new in repls:
    if old not in t:
        print("MISSING:", old[:120])
    else:
        t = t.replace(old, new)
p.write_text(t, encoding="utf-8")
print("written", p)
