"""Rebuild the JCIM submission package from final v3 and S1-S31 tables."""
from __future__ import annotations

import csv
import re
import shutil
from pathlib import Path

from docx import Document
from docx.shared import Inches, Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH
from openpyxl import Workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter
from PIL import Image

ROOT = Path(r"D:\harness\projects\d1_gramneg_vs")
V3 = ROOT / "manuscript/v3/manuscript_v3.md"
PKG = ROOT / "manuscript/submission_jcim"
SI = ROOT / "manuscript/si"
OLD = PKG / "manuscript_v3_acs.docx"
FIG = ROOT / "figures"
TOC_PNG = Path(r"D:\harness\projects\d1_manuscript_figs\figures\cell_0011_toc_graphic_jcim_exact.png")
TOC_PDF = Path(r"D:\harness\projects\d1_manuscript_figs\figures\cell_0011_toc_graphic_jcim_exact.pdf")

TITLE = "Benchmarking structure-based virtual screening across Gram-negative antibacterial targets: target-specific enrichment, co-crystal ligand pose limitations and failure modes"
RUNNING = "Gram-negative SBVS benchmark and failure modes"
AUTHORS = "Yutao Zhao^1, Hongda Liao^2,*"
AFFIL = "^1 Southern University of Science and Technology, Shenzhen, Guangdong, China; ^2 Shenzhen Bay Laboratory, Shenzhen, Guangdong, China"
CORR = "Corresponding author: Hongda Liao, Shenzhen Bay Laboratory, Shenzhen, Guangdong, China; E-mail: hongdaliao@outlook.com"
ORCID = ""


def add_runs(par, text):
    """Render **bold** fragments without destroying the rest of the text."""
    parts = re.split(r"(\*\*.*?\*\*)", text)
    for part in parts:
        if not part:
            continue
        bold = part.startswith("**") and part.endswith("**")
        run = par.add_run(part[2:-2] if bold else part)
        run.bold = bold


def add_table(doc, rows):
    if not rows:
        return
    ncol = max(len(r) for r in rows)
    table = doc.add_table(rows=len(rows), cols=ncol)
    table.style = "Table Grid"
    for i, row in enumerate(rows):
        for j in range(ncol):
            cell = table.rows[i].cells[j]
            cell.text = row[j] if j < len(row) else ""
            if i == 0:
                for p in cell.paragraphs:
                    for r in p.runs:
                        r.bold = True
    if ncol == 7:
        widths = [Inches(1.15), Inches(0.45), Inches(0.70), Inches(1.45),
                  Inches(1.45), Inches(0.50), Inches(0.45)]
        table.autofit = False
        for row in table.rows:
            for cell, width in zip(row.cells, widths):
                cell.width = width
                for par in cell.paragraphs:
                    for run in par.runs:
                        run.font.size = Pt(8)


def build_docx():
    lines = V3.read_text(encoding="utf-8").splitlines()
    doc = Document()
    sec = doc.sections[0]
    sec.top_margin = sec.bottom_margin = sec.left_margin = sec.right_margin = Inches(1)
    style = doc.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(10.5)
    title = doc.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = title.add_run(TITLE)
    run.bold = True
    run.font.size = Pt(15)
    for line in [f"Running title: {RUNNING}",
                 f"Authors: {AUTHORS}", f"Affiliations: {AFFIL}", CORR]:
        p = doc.add_paragraph(line)
        p.paragraph_format.space_after = Pt(2)

    i = 0
    frontmatter_prefixes = (
        "Manuscript version:", "**Running title:**", "**Authors:**",
        "**Affiliations:**", "**Corresponding author:**", "**ORCID:**",
    )
    while i < len(lines):
        line = lines[i].rstrip()
        if not line:
            i += 1
            continue
        if line.startswith(frontmatter_prefixes):
            i += 1
            continue
        if line.startswith("# "):  # title already added
            i += 1
            continue
        if line.startswith("## "):
            doc.add_heading(line[3:].strip(), level=1)
            i += 1
            continue
        if line.startswith("### "):
            doc.add_heading(line[4:].strip(), level=2)
            i += 1
            continue
        if line.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-{2,}:?", c) for c in cells):
                    rows.append(cells)
                i += 1
            add_table(doc, rows)
            continue
        # ordinary paragraph
        if line.startswith("**") and line.endswith("**") and len(line) < 180:
            p = doc.add_paragraph(); add_runs(p, line); 
        else:
            p = doc.add_paragraph(); add_runs(p, line)
        p.paragraph_format.space_after = Pt(4)
        i += 1
    if TOC_PNG.exists():
        doc.add_page_break()
        pic = doc.add_paragraph()
        pic.alignment = WD_ALIGN_PARAGRAPH.CENTER
        pic.add_run().add_picture(str(TOC_PNG), width=Inches(3.25), height=Inches(1.75))
        cap = doc.add_paragraph("For Table of Contents Only")
        cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    doc.save(PKG / "manuscript_v3_acs.docx")
    return PKG / "manuscript_v3_acs.docx"


def build_cover():
    text = """**Cover letter — Journal of Chemical Information and Modeling**

Dear Editor,

We submit our manuscript **"Benchmarking structure-based virtual screening across Gram-negative antibacterial targets: target-specific enrichment, co-crystal ligand pose limitations and failure modes"** for consideration as an Article in the *Journal of Chemical Information and Modeling*.

This manuscript has not been submitted elsewhere, all authors have approved the submitted version, and no preprint has been posted.

**What the paper reports.** We organise public Gram-negative antibacterial targets into two complementary regimes with regime-appropriate endpoints: a potency-driven nine-target panel (1,748 curated actives; 28,136 successfully docked ligands) evaluated by enrichment, and a structure-rich co-crystal ligand pose benchmark (40 poses across five of nine membrane-target entries) evaluated by pose geometry and interface recovery. We report three benchmark/tool failure modes — metal coordination, missing boron atom typing in Vina 1.2.5, and DUD-E-style synthetic-decoy bias — together with per-target estimates, observed decoy ratios and cluster-aware confidence intervals.

**Method-class comparison.** A post-hoc, subset-balanced comparison adds five method classes (Vina, Uni-Dock, DiffDock, Boltz-2 and DrugCLIP) plus the full-panel ECFP4 2D reference. No single engine dominates; all subset numbers are descriptive and not comparable to the full panel. DrugCLIP is a pretrained model (CC BY-NC 4.0) and train-test overlap cannot be excluded; Uni-Dock has positive-affinity/FabI missingness and a TEM1 boron surrogate.

**Negative results are retained.** A pre-registered MAES calibration gate is NO_GO on the decisive |rho| criterion (upper CI 0.161 < 0.40) and inconclusive on AUROC (19 positives; complete-case n = 132, AUROC 0.434); a pre-registered LIT-PCBA external-consistency check does not robustly reach AC2 (2/9 = 22.2%; 3/10 with TP53) and fails AC3 under both SD definitions (28.3%/34.64%). No thresholds, methods or observations were changed after observing these results.

**Data and code.** All curated data, decoy manifests, docking scores, analysis scripts and figure-generation code will be released upon publication; no proprietary or institution-restricted data were used and no wet-lab experiments were performed.

**Authorship and conflicts.** Yutao Zhao (first author; Southern University of Science and Technology, Shenzhen, Guangdong, China; ORCID 0000-0002-9917-3264) and Hongda Liao (corresponding author; Shenzhen Bay Laboratory, Shenzhen, Guangdong, China; hongdaliao@outlook.com; ORCID 0009-0001-6861-7600). The authors declare no competing interests.

Thank you for considering our manuscript.

Date: 2026-11-20

Sincerely,
Hongda Liao
Shenzhen Bay Laboratory
hongdaliao@outlook.com
"""
    md = PKG / "cover_letter_jcim.md"; md.write_text(text, encoding="utf-8")
    doc = Document()
    for block in text.split("\n\n"):
        p = doc.add_paragraph(); add_runs(p, block.replace("\n", " "))
        p.paragraph_format.space_after = Pt(6)
    doc.save(PKG / "cover_letter_jcim.docx")


def build_xlsx():
    files = sorted(SI.glob("table_S*.csv"))
    wb = Workbook()
    ws = wb.active
    ws.title = "README"
    ws.append(["sheet", "source_csv", "n_rows", "n_cols", "note"])
    used = {"README"}
    for f in files:
        name = f.name.replace("table_", "").replace(".csv", "")
        name = re.sub(r"[^A-Za-z0-9_]", "_", name)[:31]
        base = name
        k = 2
        while name in used:
            name = (base[:28] + f"_{k}")[:31]; k += 1
        used.add(name)
        sheet = wb.create_sheet(name)
        with f.open(encoding="utf-8-sig", newline="") as fh:
            rows = list(csv.reader(fh))
        for row in rows:
            sheet.append(row)
        for cell in sheet[1]:
            cell.font = Font(bold=True)
        sheet.freeze_panes = "A2"
        if sheet.max_column:
            for j in range(1, min(sheet.max_column, 20) + 1):
                sheet.column_dimensions[get_column_letter(j)].width = 16
        ws.append([name, f.name, sheet.max_row, sheet.max_column,
                   "final v3 submission table" if name != "README" else ""])
    wb.save(PKG / "SI_tables_S1-S31.xlsx")
    return len(wb.sheetnames)


def copy_figures():
    dest = PKG / "figures"; dest.mkdir(exist_ok=True)
    pairs = {
        "fig1": (FIG / "fig1_data_landscape.png", FIG / "fig1_data_landscape.pdf"),
        "fig2": (FIG / "v3/fig2_two_regime_panel.png", FIG / "v3/fig2_two_regime_panel.pdf"),
        "fig3": (FIG / "v3/fig3_trackA_enrichment.png", FIG / "v3/fig3_trackA_enrichment.pdf"),
        "fig4": (FIG / "v3/fig4_three_axis_transfer.png", FIG / "v3/fig4_three_axis_transfer.pdf"),
        "fig5": (FIG / "v3/fig5_failure_modes.png", FIG / "v3/fig5_failure_modes.pdf"),
        "fig6": (FIG / "v3/fig6_maes_calibration.png", FIG / "v3/fig6_maes_calibration.pdf"),
        "fig7": (FIG / "v3/fig7_method_axis.png", FIG / "v3/fig7_method_axis.pdf"),
    }
    for name, (png, pdf) in pairs.items():
        if png.exists(): shutil.copy2(png, dest / f"{name}.png")
        if pdf.exists(): shutil.copy2(pdf, dest / f"{name}.pdf")
    if TOC_PNG.exists():
        shutil.copy2(TOC_PNG, dest / "toc_graphic.png")
        Image.open(TOC_PNG).convert("RGB").save(
            dest / "toc_graphic.tif", format="TIFF", dpi=(300, 300), compression="tiff_deflate")
    if TOC_PDF.exists():
        shutil.copy2(TOC_PDF, dest / "toc_graphic.pdf")


def main():
    docx = build_docx()
    build_cover()
    nsheets = build_xlsx()
    copy_figures()
    print("docx", docx.stat().st_size, "xlsx sheets", nsheets)
    print("package", [p.name for p in sorted(PKG.iterdir())])


if __name__ == "__main__":
    main()
