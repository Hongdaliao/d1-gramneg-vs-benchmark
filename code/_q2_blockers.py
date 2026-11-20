from pathlib import Path

# Fig 6 caption + deprecation note
Path("figures/fig6_w12_calibration_formal_caption.md").write_text(
    """# Fig 6 — calibration-gate negative result (final)

**Official file:** `figures/fig6_w12_calibration_formal.png` (300 dpi; PDF: `figures/fig6_w12_calibration_formal.pdf`).

**Source:** W12 formal calibration gate, manifest `results/maes/w12/manifest_features.csv` sha256 `e4ae5c23461474f773a016ec415e9da8e5dd6ff10bbef55eda73c726cb4fa3ca` (146 pose observations; 40 systems; vina 40 / vinardo 40 / vina_afdb 40 / boltz2 26; 19 evaluation positives).

**Primary V+M (S neutral):** AUROC 0.496 (wild cluster CI 0.318–0.672), |ρ(MAES,RMSD)| 0.075 (cluster CI −0.302–0.161), Brier 0.145, ECE 0.124 → **NO_GO** (criteria AUROC≥0.65 and |ρ|≥0.40 not met). Full V+M+S is exploratory only.

**Caption:** Calibration-gate negative result for MAES. Reliability (left) and formal metrics (right) for the pre-registered W12 gate (146 pose observations; V+M primary). Both pre-registered criteria fail; MAES is retained only as an exploratory/descriptive pose-quality analysis, not as a calibrated probability or ranking score.

**Deprecated:** `figures/fig6_maes_calibration_draft.*` (n=36/AUROC 0.586, preliminary 2026-11-17) must not be used for Fig 6.
""",
    encoding="utf-8",
)
Path("figures/fig6_maes_calibration_draft_DEPRECATED.md").write_text(
    """# DEPRECATED — old Fig 6 draft (n=36)

This is the preliminary n=36 Vina-only Fig 6 draft (AUROC 0.586, 2026-11-17). It is **superseded** by the W12 formal calibration gate (146 poses; AUROC 0.496; |ρ| 0.075; NO_GO).

Use `figures/fig6_w12_calibration_formal.png` / `.pdf` and `figures/fig6_w12_calibration_formal_caption.md` instead. Do not cite or embed the draft figure in the manuscript.
""",
    encoding="utf-8",
)

# main skeleton Fig 6 row: add official source path
p = Path("manuscript/main_skeleton.md")
t = p.read_text(encoding="utf-8")
old = "| Fig 6 | **Calibration-gate negative result for MAES (reliability/Brier/ECE; AUROC 0.496; |ρ| 0.075; 146 poses)** — exploratory/descriptive only, not a probability score or ranking figure | WP4 | **NO_GO; retained as negative result** |"
new = "| Fig 6 | **Calibration-gate negative result for MAES (reliability/Brier/ECE; AUROC 0.496; |ρ| 0.075; 146 poses)** — exploratory/descriptive only, not a probability score or ranking figure | `figures/fig6_w12_calibration_formal.png` / `.pdf` | **NO_GO; retained as negative result; old draft deprecated** |"
assert old in t
p.write_text(t.replace(old, new), encoding="utf-8")

# results_trackb R2 title/ref
p = Path("manuscript/results_trackb_draft.md")
t = p.read_text(encoding="utf-8")
t = t.replace("## R2 Primary MAES evaluation (preliminary first-pass; W12 formal NO-GO; Fig 6a–b; Table 2)",
              "## R2 Exploratory MAES evaluation (calibration gate NO-GO; Fig 6; Table 2)")
t = t.replace("## R1 Analysis set and data quality (Fig 5a; Table 2)",
              "## R1 Analysis set and data quality (Fig 5; Table 2)")
p.write_text(t, encoding="utf-8")

# discussion TODO / [[FILL]]
p = Path("manuscript/discussion_draft.md")
t = p.read_text(encoding="utf-8")
t = t.replace("> 状态：骨架（等 anchor / 多方法 / 三轴结果填入 `[[FILL]]`）。数字以 `figures/trackA_enrichment_full.json`、`manuscript/results_trackA_draft.md` v1.1、`manuscript/results_trackb_draft.md` 为准。",
              "> 状态：D3b/R6 与 W12 已填；数字以 `figures/trackA_enrichment_full.json`、`manuscript/results_trackA_draft.md`、`manuscript/results_trackb_draft.md` 为准。")
t = t.replace("## [[FILL]] 清单", "## 内部待办（不投稿）")
t = t.replace("- 时间/结构来源轴的正式指标仍在补（novel-scaffold 178/178 已可关联、AFDB 臂在跑）。",
              "- 时间/结构来源轴的正式指标已完成（R6；`figures/task35_validation.json` 12/12）。")
p.write_text(t, encoding="utf-8")

# methods_maes: remove probability claim and unify W12 sample
p = Path("manuscript/methods_maes_draft.md")
t = p.read_text(encoding="utf-8")
t = t.replace("MAES converts these three axes into a probability that a pose is acceptable.",
              "MAES combines these three axes into a descriptive pose-quality score (not a calibrated probability; see W12 NO-GO).")
t = t.replace("> v1.2 changes (W12 formal): **W12 calibration gate returned NO-GO** (n=146 poses; primary V+M AUROC 0.496, |ρ| 0.075); MAES is reported **descriptively only**, not as a probability or ranking score. v1.1 changes (Lead G0 decision): MAES main analysis runs on **Vina docking poses** (40 inhibitor systems, full coverage); Boltz-2 (≤400 tokens) is a **secondary arm** reported separately; release wording changed to \"will be released upon publication\" (P9 gate); limitation sentences added below (§4.8 word target kept 300–400).",
              "> v1.2 changes (W12 formal): **W12 calibration gate returned NO-GO** (n=146 poses; primary V+M AUROC 0.496, |ρ| 0.075); MAES is reported **descriptively only**, not as a probability or ranking score. W12 calibration data are **146 pose observations across 4 methods/arms (vina 40, vinardo 40, vina_afdb 40, boltz2 26)**; the earlier Vina-only 40-system arm is historical and not the W12 sample. v1.1 changes (Lead G0 decision): release wording changed to \"will be released upon publication\" (P9 gate); limitation sentences added below (§4.8 word target kept 300–400).")
p.write_text(t, encoding="utf-8")

print("Fig6 caption/deprecation and manuscript minor fixes applied")
