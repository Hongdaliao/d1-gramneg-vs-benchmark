from pathlib import Path

# P4 wording
p = Path("manuscript/methods_maes_draft.md")
t = p.read_text(encoding="utf-8")
t = t.replace("MAES combines these three axes into a descriptive pose-quality score (not a calibrated probability; see W12 NO-GO).",
              "MAES combines these three axes into a descriptive pose-quality score; a candidate probability score was pre-registered, but its W12 gate failed (AUROC 0.496, |ρ| 0.075), so it is reported descriptively only.")
p.write_text(t, encoding="utf-8")

# P6 drop fill-map section (internal ingest keys, not submission text)
p = Path("manuscript/results_trackb_draft.md")
lines = p.read_text(encoding="utf-8").splitlines()
out = []
skip = False
for line in lines:
    if line.strip() == "### Fill map (ingest JSON keys)":
        skip = True
        continue
    if skip:
        if line.strip().startswith("|"):
            continue
        if line.strip() == "":
            continue
        if line.strip().startswith("**"):
            skip = False
            out.append(line)
            continue
        skip = False
    out.append(line)
p.write_text("\n".join(out) + "\n", encoding="utf-8")

# P7 figures README note
readme = Path("figures/README.md")
old = readme.read_text(encoding="utf-8") if readme.exists() else "# Figures\n"
old += """
## Fig 6 (MAES calibration gate)

- **Use in manuscript:** `figures/fig6_w12_calibration_formal.png` / `.pdf`; caption `figures/fig6_w12_calibration_formal_caption.md` (W12 final: n=146, AUROC 0.496, |ρ| 0.075, NO_GO).
- **History only — do not use in manuscript:** `figures/fig6_maes_calibration_draft.*` (old n=36/AUROC 0.586 preliminary draft; marked deprecated in `fig6_maes_calibration_draft_DEPRECATED.md`).
"""
readme.write_text(old, encoding="utf-8")
print("P4/P6/P7 fixed")
