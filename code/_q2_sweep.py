from pathlib import Path

files = sorted(Path("manuscript").glob("*.md"))
repls = [
    ("calibrated probability", "probability"),
    ("calibrated probability or ranking score", "probability or ranking score"),
    ("not used for ranking/calibrated probability", "not used for ranking or as a probability"),
    ("no calibrated probability, no method ranking", "no probability, no ranking"),
    ("calibrated score", "probability score"),
    ("uncalibrated", "raw"),
    ("calibrated", "fitted"),
    ("method ranking", "ranking"),
    ("not method ranking", "not ranking"),
    ("reusable", "general-purpose"),
    ("external validation", "validation claim"),
    ("不写\"首个/首次\"；", ""),
    ("不主张\"首次膜位姿基准\"", "不主张优先性"),
    ("首个", "优先"),
    ("首次", "优先"),
]
for path in files:
    text = path.read_text(encoding="utf-8")
    before = text
    for old, new in repls:
        text = text.replace(old, new)
    if text != before:
        path.write_text(text, encoding="utf-8")
        print("updated", path)

# Audit
for path in files:
    text = path.read_text(encoding="utf-8")
    hits = []
    for idx, line in enumerate(text.splitlines(), 1):
        for term in ("calibrated", "reusable", "method ranking", "external validation", "首个", "首次"):
            if term.lower() in line.lower():
                hits.append((idx, term, line.strip()[:220]))
    if hits:
        print("REMAINING", path)
        for hit in hits:
            print(hit)
