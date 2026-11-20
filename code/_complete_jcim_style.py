from pathlib import Path

p = Path(r"manuscript/v3/manuscript_v3_jcim_style.md")
t = p.read_text(encoding="utf-8")
repls = [
("### Data and code availability\n\nCurated data, decoy manifests, docking scores, analysis scripts, figure-generation code and intermediate tables will be released upon publication. The MAES calibration manifest is identified by sha256 e4ae5c23461474f773a016ec415e9da8e5dd6ff10bbef55eda73c726cb4fa3ca. SI tables S1–S31 document the panel, enrichments, MAES-gate observations and related analyses. The repository also includes the raw score tables and coverage files used for every pooled estimate.",
 "### Data and software release\n\nThe released data, decoy manifests, docking scores, analysis scripts, figure-generation code and intermediate tables are described in the back-matter Data and code availability statement; SI tables S1–S31 accompany the submission."),
("paired comparisons resample matched subset pairs rather than comparing independent confidence intervals. No multiple-comparison correction or significance claim is made; per-target estimates are the primary descriptive unit.",
 "paired comparisons resample matched subset pairs rather than comparing independent confidence intervals. No multiple-comparison correction or significance claim is made; per-target estimates are the primary descriptive unit. We report minimum detectable effect sizes (MDES) rather than post-hoc power: Track A detects approximately 0.93–0.99 × between-target SD; Track B effective n_eff is 20–35 and detects only ≥15–25 percentage-point acceptable-rate differences, or MAES differences of about 0.10–0.16."),
("The single-chain arm changed DiffDock AUC from 0.488 to 0.667 on GyrA and from 0.469 to 0.444 on FabI; it changes only the DiffDock receptor representation, not the Vina ranking.",
 "The single-chain arm changed DiffDock AUC from 0.488 to 0.667 on GyrA (Δ 0.178; 95% CI 0.021–0.362) and from 0.469 to 0.444 on FabI (Δ −0.025; 95% CI −0.188 to 0.146); it changes only the DiffDock receptor representation, not the Vina ranking. No single engine dominated across targets."),
("AC3 failed under both SD definitions (maximum CV 28.3% population SD and 34.64% sample SD, MTORC1), while AC4, AC6 and AC7 passed;",
 "AC3 failed under both SD definitions (maximum CV 28.28%, reported as 28.3% under population SD, and 34.64% under sample SD, MTORC1), while AC4, AC6 and AC7 passed;"),
("The benchmark response is to report per-target estimates, observed decoy ratios, coverage, receptor QC and pre-registered gates, and to retain negative results rather than selecting a single score.",
 "The benchmark response is to report per-target estimates, observed decoy ratios, coverage, receptor QC and pre-registered gates, and to retain negative results rather than selecting a single score. No post-hoc threshold, method or observation was changed after these results were observed."),
]
for old, new in repls:
    if old not in t:
        print("MISSING:", old[:100])
    else:
        t = t.replace(old, new, 1)
p.write_text(t, encoding="utf-8")
print("patched", p.stat().st_size, len(t.split()), "words")
