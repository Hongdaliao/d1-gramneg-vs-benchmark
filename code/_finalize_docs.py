from pathlib import Path

# --- package README / checklist size and page corrections ---
for f in (r"manuscript/submission_jcim/README_submission_jcim.md",
          r"manuscript/submission_jcim/SUBMISSION_CHECKLIST_JCIM.md"):
    p = Path(f); t = p.read_text(encoding="utf-8")
    t = t.replace("17 pages", "15 pages")
    t = t.replace("54 sheets", "50 sheets")
    t = t.replace("README + 53 tables", "README + 49 table sheets")
    p.write_text(t, encoding="utf-8")

# --- review changelog / mapping final stats ---
p = Path(r"notes/review/jcim_style_changelog.md"); t = p.read_text(encoding="utf-8")
t = t.replace("current: 53,139 B / sha `00622cdce3925785006a97634ba935ddb2495a6ab0934cca0d06f0a35d997cb2` / 7,539 words.",
              "current: 53,155 B / sha `f3144889f3368a5c37205a0e7448dab640a3859751df418156b7738fbcebc0b3` / 7,541 words.")
t = t.replace("- **M7** Removed duplicate sentences:", "- **M7 residual:** CPU-cost duplicate merged; interface-bound caveat restored once in Track B. Removed duplicate sentences:")
p.write_text(t, encoding="utf-8")

p = Path(r"notes/review/jcim_style_rewrite_mapping.md"); t = p.read_text(encoding="utf-8")
t = t.replace("Draft total: **7,539 words** vs v3 **8,436 words** (−10.6%; inside ±15%).",
              "Draft total: **7,541 words** vs v3 **8,436 words** (−10.6%; inside ±15%).")
t = t.replace("Total words 7,539", "Total words 7,541")
p.write_text(t, encoding="utf-8")

# --- v3 changelog append ---
p = Path(r"manuscript/v3/V3_CHANGELOG.md"); t = p.read_text(encoding="utf-8")
t += """

## JCIM-style replacement (2026-11-20)

- Replaced the previous signed v3 with the E-approved JCIM-style submission (`manuscript_v3.md`).
- Structure: Abstract 239 words; Introduction 3 paragraphs; Results 5 subsections (3/4/3/3/3); Discussion 9 paragraphs; Methods 12 subsections after Discussion; no standalone Limitations/Conclusion.
- References renumbered by first appearance (39/39); Tables 1–4 and Figures 1–7 retained.
- **One numerical direction correction** relative to the previous signed v3 L119: LpxC Zn control is **0.220 (Zn-retained) → 0.222 (Zn-stripped)**; the old wording had the direction reversed.
- Old signed v3 archived at `manuscript/v3/archive/manuscript_v3_signed_7c7033_20261120.md` (60,167 B / sha `7c703373…`).
- Submission package rebuilt from the new baseline.
"""
p.write_text(t, encoding="utf-8")
print("docs updated")
