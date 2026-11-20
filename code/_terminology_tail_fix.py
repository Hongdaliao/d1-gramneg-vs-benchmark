from pathlib import Path

p = Path(r"manuscript/v3/manuscript_v3.md")
t = p.read_text(encoding="utf-8")
repls = [
    ("subset actives overlap the temporal-train reference (Addendum 12a; 68–94% for seven of eight comparison targets, TEM1 38%, LepB 100% at n=3)",
     "subset actives overlap the temporal-train reference (68–94% for seven of eight comparison targets, TEM1 38%, LepB 100% at n=3; S25)"),
    ("TEM1_Ec boron actives use the pre-specified `_B2C` surrogate (22 complexes, Addendum 14a).",
     "TEM1_Ec boron actives use the pre-specified `_B2C` surrogate (22 complexes; input substitution pre-specified, S24)."),
    ("ECFP4 was not recomputed on the subset because of temporal-train overlap (Addendum 12a; S25 footnote);",
     "ECFP4 was not recomputed on the subset because of overlap with the temporal-train reference (S25 footnote);"),
    ("`figures/v3/fig6_w12_calibration.png` / `.pdf`", "`figures/v3/fig6_maes_calibration.png` / `.pdf`"),
]
for old, new in repls:
    if old not in t:
        print("MISSING", old[:100])
    else:
        t = t.replace(old, new, 1)
p.write_text(t, encoding="utf-8")
print("written", p.stat().st_size)
