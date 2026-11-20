from pathlib import Path

p = Path(r"manuscript/v3/manuscript_v3_jcim_style.md")
t = p.read_text(encoding="utf-8")

# m7 residual: merge the two CPU-cost sentences into one.
old = "The practical cost is asymmetric: CPU-only Vina remains inexpensive, while GPU arms require scheduling, model-weight provenance and deterministic environment disclosure. Vina remains cheap on CPU, whereas DiffDock and Boltz-2 require GPU time, and Boltz-2's nested subset is deliberately small."
new = "The practical cost is asymmetric: CPU-only Vina remains inexpensive, whereas DiffDock and Boltz-2 require GPU time, and Boltz-2's nested subset is deliberately small; GPU arms add scheduling, model-weight provenance and deterministic environment disclosure."
if old not in t:
    print("MISSING CPU duplicate")
else:
    t = t.replace(old, new, 1)

# minor: restore interface-bound caveat once.
old2 = "AcrB contributes 34 of the 40 main-panel poses, so apparent ligand-class effects are also target-composition effects; the class medians are descriptive and not tested."
new2 = "Only 3 of 40 poses are interface-bound. AcrB contributes 34 of the 40 main-panel poses, so apparent ligand-class effects are also target-composition effects; the class medians are descriptive and not tested."
if old2 not in t:
    print("MISSING interface caveat anchor")
else:
    t = t.replace(old2, new2, 1)

# final version line for the replacement file
t = t.replace("Manuscript version: v3 JCIM-style (2026-11-20)",
              "Manuscript version: v3-submission (2026-11-20)", 1)
p.write_text(t, encoding="utf-8")
print("patched", p.stat().st_size, len(t.split()), "words")
