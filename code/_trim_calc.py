from pathlib import Path
old = Path(r"manuscript/v3/manuscript_v3.md").read_text(encoding="utf-8").split("## Abstract")[1].split("## Keywords")[0].strip()
repls = [
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
]
t = old
print("base", len(t.split()))
for i, (a, b) in enumerate(repls):
    if a not in t:
        print("missing", i, a[:40]); continue
    before = len(t.split()); t = t.replace(a, b, 1); after = len(t.split())
    print(i, "save", before-after, "now", after)
print("final", len(t.split()), "with sentence", len((t + " A post-hoc method comparison found no single engine dominated.").split()))
