from pathlib import Path
import re

ROOT = Path(r"D:\harness\projects\d1_gramneg_vs")
p = ROOT / "manuscript/v3/manuscript_v3_jcim_style.md"
text = p.read_text(encoding="utf-8")
body, rest = text.split("## References", 1)
refs_part, caps_part = rest.split("## Figure captions", 1)

# parse old references
entries = {}
for line in refs_part.splitlines():
    m = re.match(r"\[(\d+)\]\s+(.*)", line)
    if m:
        entries[int(m.group(1))] = m.group(2).strip()

# parse citations and first-occurrence order
def expand(tok):
    out = []
    for part in tok.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        elif part:
            out.append(int(part))
    return out

order = []
for m in re.finditer(r"\[([0-9,\-]+)\]", body):
    for n in expand(m.group(1)):
        if n not in order:
            order.append(n)
missing = [i for i in range(1, 40) if i not in order]
if missing:
    raise SystemExit(f"missing old refs in draft: {missing}")
mapping = {old: i + 1 for i, old in enumerate(order)}

def compress(nums):
    nums = sorted(set(nums)); chunks = []
    i = 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        if j - i >= 2:
            chunks.append(f"{nums[i]}-{nums[j]}")
        else:
            chunks.extend(str(n) for n in nums[i:j + 1])
        i = j + 1
    return ",".join(chunks)

def repl(m):
    nums = [mapping[n] for n in expand(m.group(1))]
    return "[" + compress(nums) + "]"

body2 = re.sub(r"\[([0-9,\-]+)\]", repl, body)
refs_lines = ["## References", ""]
for old in order:
    refs_lines.append(f"[{mapping[old]}] {entries[old]}")
    refs_lines.append("")
new_text = body2.rstrip() + "\n\n" + "\n".join(refs_lines).rstrip() + "\n\n## Figure captions" + caps_part
p.write_text(new_text, encoding="utf-8")
print("renumbered", len(order), "refs; new first/last", mapping[order[0]], mapping[order[-1]])
