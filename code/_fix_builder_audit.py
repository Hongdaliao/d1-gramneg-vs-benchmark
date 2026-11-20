from pathlib import Path

p = Path(r"code/maes/build_submission_jcim_final.py")
t = p.read_text(encoding="utf-8")

# front matter constants: no version line and no ORCID on the title page
t = t.replace('CORR = "Corresponding author: Hongda Liao, Shenzhen Bay Laboratory, Shenzhen, Guangdong, China; E-mail: hongdaliao@outlook.com; ORCID 0009-0001-6861-7600"',
              'CORR = "Corresponding author: Hongda Liao, Shenzhen Bay Laboratory, Shenzhen, Guangdong, China; E-mail: hongdaliao@outlook.com"')
t = t.replace('ORCID = "ORCID: Yutao Zhao 0000-0002-9917-3264; Hongda Liao 0009-0001-6861-7600"',
              'ORCID = ""')
t = t.replace('FIG = ROOT / "figures"', 'FIG = ROOT / "figures"\nTOC_PNG = Path(r"D:\\harness\\projects\\d1_manuscript_figs\\figures\\cell_0011_toc_graphic_jcim_exact.png")\nTOC_PDF = Path(r"D:\\harness\\projects\\d1_manuscript_figs\\figures\\cell_0011_toc_graphic_jcim_exact.pdf")')
t = t.replace('    for line in [f"Manuscript version: v3-submission (2026-11-20)", f"Running title: {RUNNING}",\n                 f"Authors: {AUTHORS}", f"Affiliations: {AFFIL}", CORR, ORCID]:',
              '    for line in [f"Running title: {RUNNING}",\n                 f"Authors: {AUTHORS}", f"Affiliations: {AFFIL}", CORR]:')

# table 2: seven columns, tighter widths and 8 pt for legibility
old = '''    for i, row in enumerate(rows):
        for j in range(ncol):
            cell = table.rows[i].cells[j]
            cell.text = row[j] if j < len(row) else ""
            if i == 0:
                for p in cell.paragraphs:
                    for r in p.runs:
                        r.bold = True
'''
new = '''    for i, row in enumerate(rows):
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
'''
if old not in t:
    print("MISSING table width anchor")
else:
    t = t.replace(old, new, 1)

# copy exact TOC and fix fig6 source path
old_copy = '''        "fig6": (FIG / "v3/fig6_w12_calibration.png", FIG / "v3/fig6_w12_calibration.pdf"),'''
new_copy = '''        "fig6": (FIG / "v3/fig6_maes_calibration.png", FIG / "v3/fig6_maes_calibration.pdf"),'''
t = t.replace(old_copy, new_copy, 1)
old_toc = '''    if (dest / "toc_graphic.png").exists():
        pass'''
new_toc = '''    if TOC_PNG.exists():
        shutil.copy2(TOC_PNG, dest / "toc_graphic.png")
        Image.open(TOC_PNG).convert("RGB").save(
            dest / "toc_graphic.tif", format="TIFF", dpi=(300, 300), compression="tiff_deflate")
    if TOC_PDF.exists():
        shutil.copy2(TOC_PDF, dest / "toc_graphic.pdf")'''
t = t.replace(old_toc, new_toc, 1)

# insert TOC page at the end of DOCX
old_end = '''    doc.save(PKG / "manuscript_v3_acs.docx")
    return PKG / "manuscript_v3_acs.docx"'''
new_end = '''    if TOC_PNG.exists():
        doc.add_page_break()
        pic = doc.add_paragraph()
        pic.alignment = WD_ALIGN_PARAGRAPH.CENTER
        pic.add_run().add_picture(str(TOC_PNG), width=Inches(3.25), height=Inches(1.75))
        cap = doc.add_paragraph("For Table of Contents Only")
        cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    doc.save(PKG / "manuscript_v3_acs.docx")
    return PKG / "manuscript_v3_acs.docx"'''
t = t.replace(old_end, new_end, 1)

# add PIL import
if 'from PIL import Image' not in t:
    t = t.replace('from openpyxl.utils import get_column_letter', 'from openpyxl.utils import get_column_letter\nfrom PIL import Image')

p.write_text(t, encoding="utf-8")
print("builder patched")
