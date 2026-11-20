from pathlib import Path

p = Path(r"manuscript/v3/manuscript_v3.md")
t = p.read_text(encoding="utf-8")
repls = [
    ("multimodal pretraining [26-28]", "multimodal pretraining [27,28]"),
    ("hit identification [31-33]", "hit identification [32,33]"),
    ("time-split docking evaluations [37,38] and AI-docking context [39]", "docking evaluation context [38,39]"),
    ("physics checks [22-24]", "physics checks [23,24]"),
    ("DiffDock higher for GyrB_Ec 0.682 vs 0.432", "DiffDock higher for GyrB_Ec 0.682 vs 0.420"),
    ("Vina higher for LpxC_Pa 0.857 vs 0.437 and GyrA_Ec 0.874 vs 0.488", "Vina higher for LpxC_Pa 0.857 vs 0.437 and GyrA_Ec 0.869 vs 0.488"),
    ("Vina on the same paired ligands;", "Vina on the same 806-ligand subset (pairing is used only for the paired ΔAUC);"),
    ("+0.089 (95% CI -0.120 to 0.302)", "+0.089 (95% CI -0.120 to 0.303)"),
    ("S26b, S27b, S31b, S31c and S31d", "S26b, S27b, S31b, S31c, S31d and S31e"),
    ("Source: `results/trackA/method_axis/analysis/fig7_data.csv`; figure file `figures/v3/fig7_method_axis.png` / `.pdf`.",
     "Source: `results/trackA/method_axis/analysis/fig7_data.csv`, `table_S31c_method_coverage.csv` (panel F) and `table_S31e_common7_sensitivity.csv` (panel E); figure file `figures/v3/fig7_method_axis.png` / `.pdf`."),
]
for old, new in repls:
    if old not in t:
        print("MISSING", old[:120])
    else:
        t = t.replace(old, new, 1)
# blank line between refs [20] and [21]
lines = t.splitlines()
out = []
for line in lines:
    out.append(line)
    if line.startswith("[20] Jia Y"):
        out.append("")
t = "\n".join(out) + ("\n" if t.endswith("\n") else "")
p.write_text(t, encoding="utf-8")
print("v3 E fixes applied", p.stat().st_size)
