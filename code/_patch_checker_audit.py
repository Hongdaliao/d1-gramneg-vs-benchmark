from pathlib import Path

p = Path(r"notes/review/_check_jcim_style.py")
t = p.read_text(encoding="utf-8")
t = t.replace(
'''required = ["Abstract", "Keywords", "Introduction", "Results", "Discussion", "Methods",
            "Author contributions", "Data and code availability", "Use of artificial intelligence", "References", "Figure captions"]''',
'''required = ["Abstract", "Keywords", "Introduction", "Results", "Discussion", "Methods",
            "Author contributions", "Acknowledgments", "Competing interests",
            "Data and Software Availability", "Supporting Information", "References", "Figure captions"]''')
t = t.replace(
'''checks["Methods data-software heading"] = "### Data and software release" in t
checks["backmatter data availability only"] = t.count("## Data and code availability") == 1''',
'''checks["Methods data-software heading"] = "### Data and software release" in t
checks["backmatter data availability only"] = t.count("## Data and Software Availability") == 1
checks["no old data availability heading"] = "## Data and code availability" not in t''')
t = t.replace(
'''checks["m7 title version"] = "Manuscript version: v3-submission (2026-11-20)" in t''',
'''checks["m7 title clean"] = "Manuscript version:" not in t and "**ORCID:**" not in t''')
t = t.replace(
'''checks["m8 6-418 exact range"] = "scored actives per target range from 6 to 418 after exclusion" in t''',
'''checks["m8 6-418 exact range"] = "scored actives per target range from 6 to 418 after exclusion" in t
checks["Table 2 seven columns"] = table_rows(2) >= 0 and len([l for l in t[t.index("**Table 2."):].splitlines() if l.startswith("|")][0].split("|")) - 2 == 7
checks["captions no backticks"] = "`" not in t.split("## Figure captions")[1]
checks["captions no source tails"] = "Source:" not in t.split("## Figure captions")[1]''')
p.write_text(t, encoding="utf-8")
print("checker patched")
