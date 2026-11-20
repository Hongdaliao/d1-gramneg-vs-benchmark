from pathlib import Path

for f in (r"manuscript/submission_jcim/README_submission_jcim.md",
          r"manuscript/submission_jcim/SUBMISSION_CHECKLIST_JCIM.md"):
    p = Path(f); t = p.read_text(encoding="utf-8")
    t = t.replace("f3144889f3368a5c37205a0e7448dab640a3859751df418156b7738fbcebc0b3",
                  "d6a318d3b94999172a6d434f376718739886fd53f19a721290c8de6aa808a1a8")
    t = t.replace("f3144889", "d6a318d3")
    t = t.replace("53,155 B", "52,711 B").replace("53155 B", "52711 B")
    t = t.replace("7,541 words", "7,552 words")
    t = t.replace("58062 B", "126918 B").replace("215311 B", "299404 B").replace("15 pages", "16 pages")
    t = t.replace("16 files", "17 files")
    t = t.replace("Fig1–7 PNG+PDF (300 dpi) + TOC graphic", "Fig1–7 PNG+PDF (300 dpi) + TOC (PNG/PDF/TIFF)")
    t = t.replace("Fig1–7 + TOC", "Fig1–7 + TOC (PNG/PDF/TIFF)")
    p.write_text(t, encoding="utf-8"); print("updated", p)

p = Path(r"manuscript/v3/V3_CHANGELOG.md"); t = p.read_text(encoding="utf-8")
t = t.replace("53,155 B", "52,711 B").replace("7,541 words", "7,552 words").replace("f3144889", "d6a318d3")
t += """

## JCIM production audit fixes (2026-11-20)

- Table 2 header repaired to 7 columns and given explicit widths (DOCX/PDF checked on page 6).
- Added `## Supporting Information` before References.
- Back matter: Acknowledgments (Funding + AI disclosure), Competing interests, Data and Software Availability (Datasets./Software.).
- Figure captions: removed source tails and backticks.
- TOC graphic replaced by exact 3.25 × 1.75 in @300 dpi image; DOCX last page carries image + "For Table of Contents Only"; package keeps PNG/PDF/TIFF.
- References: ref14 article number 141; ref36 Science article eads9530; ref27 ICML PMLR 70:1321–1330.
- Title page: version line and ORCID line removed; ORCID stays in the submission system.
"""
p.write_text(t, encoding="utf-8"); print("changelog updated")
