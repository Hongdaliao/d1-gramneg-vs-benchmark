import re
from pathlib import Path
t = Path(r"manuscript/v3/manuscript_v3.md").read_text(encoding="utf-8")
a = t.split("## Abstract")[1].split("## Keywords")[0].strip()
print("split", len(a.split()))
print("regex", len(re.findall(r"[A-Za-z0-9][A-Za-z0-9'’\-]*", a)))
print(a)
