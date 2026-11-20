from pathlib import Path

# Review changelog
p = Path(r"notes/review/jcim_style_changelog.md"); t = p.read_text(encoding="utf-8")
t = t.replace("current: 53,155 B / sha `f3144889f3368a5c37205a0e7448dab640a3859751df418156b7738fbcebc0b3` / 7,541 words.",
              "current: 52,711 B / sha `d6a318d3b94999172a6d434f376718739886fd53f19a721290c8de6aa808a1a8` / 7,552 words.")
t += "\n\n## Production audit batch (2026-11-20)\n\n- Table 2 corrected to 7 columns, explicit widths, verified on rendered page 6.\n- Added Supporting Information section before References.\n- Back matter expanded: Acknowledgments (Funding + AI), Competing interests, Data and Software Availability (Datasets./Software.).\n- Figure captions stripped of Source tails/backticks.\n- TOC graphic: exact 3.25×1.75 in @300 dpi (975×525) image in DOCX last page + `figures/toc_graphic.tif` (also PNG/PDF).\n- Reference fields: ref14 141; ref36 eads9530; ref27 PMLR 70:1321–1330.\n- Title page version/ORCID lines removed; ORCID retained only in submission system.\n- Cover letter: no dual submission/all-author approval/no preprint; date and corresponding signature added.\n"
p.write_text(t, encoding="utf-8")

# Mapping doc
p = Path(r"notes/review/jcim_style_rewrite_mapping.md"); t = p.read_text(encoding="utf-8")
t = t.replace("7,541 words", "7,552 words").replace("7,541", "7,552")
t = t.replace("7,539 words", "7,552 words").replace("7,539", "7,552")
p.write_text(t, encoding="utf-8")
print("review docs updated")
