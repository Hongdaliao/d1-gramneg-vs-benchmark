#!/usr/bin/env python
"""LIT-PCBA anchor integrity check (Lead request 2026-11-18).

For every set/subseed score file under results/litpcba/:
  - expected rows == manifest rows (manifest_s<sub>.csv), n_act/n_inact match
  - no NUL bytes / blank lines; status distribution; duplicate ligand ids
  - SHA256 of score file
Outputs: results/litpcba/anchor_integrity_check.json
Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\_b_anchor_integrity.py
"""
import csv, glob, hashlib, json, os, time

ROOT = r"D:\harness\projects\d1_gramneg_vs"
OUT = os.path.join(ROOT, "results", "litpcba")
SUBS = [42, 43, 44]

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def check_file(path, man_path):
    issues = []
    raw = open(path, "rb").read()
    if b"\x00" in raw:
        issues.append("NUL_bytes")
    try:
        rows = list(csv.DictReader(open(path, encoding="utf-8-sig")))
    except Exception as e:
        return {"path": path, "rows": None, "issues": issues + [f"csv_error:{type(e).__name__}"], "sha256": sha(path)}
    man = list(csv.DictReader(open(man_path, encoding="utf-8-sig"))) if os.path.exists(man_path) else []
    n_man = len(man)
    n_act_man = sum(1 for r in man if str(r["label"]) == "1")
    n_act = sum(1 for r in rows if str(r["label"]) == "1")
    ids = [r.get("ligand_id") for r in rows]
    if len(rows) != n_man:
        issues.append(f"row_mismatch:{len(rows)}vs{n_man}")
    if n_act != n_act_man:
        issues.append(f"act_mismatch:{n_act}vs{n_act_man}")
    if len(set(ids)) != len(ids):
        issues.append("duplicate_ligand_ids")
    status = {}
    for r in rows:
        status[r.get("status", "?")] = status.get(r.get("status", "?"), 0) + 1
    return {"path": path, "rows": len(rows), "n_act": n_act, "n_man": n_man, "status": status,
            "issues": issues, "sha256": sha(path)}

def main():
    rep = {"generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"), "sets": {}}
    bad = 0
    for sp in sorted(glob.glob(os.path.join(OUT, "*", "scores_s*_v42.csv"))):
        d = os.path.dirname(sp)
        s = os.path.basename(d)
        m = os.path.basename(sp).replace("scores_s", "").replace("_v42.csv", "")
        man = os.path.join(d, f"manifest_s{m}.csv")
        res = check_file(sp, man)
        rep["sets"].setdefault(s, {})[f"sub{m}"] = res
        if res.get("issues"):
            bad += 1
    rep["summary"] = {"n_files": sum(len(v) for v in rep["sets"].values()),
                      "n_with_issues": bad,
                      "n_sets": len(rep["sets"]),
                      "missing_subs": {s: [str(x) for _ in [0] for x in SUBS
                                           if not os.path.exists(os.path.join(OUT, s, f"scores_s{x}_v42.csv"))]
                                       for s in rep["sets"]}}
    p = os.path.join(OUT, "anchor_integrity_check.json")
    json.dump(rep, open(p, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(json.dumps(rep["summary"], indent=1))
    print("wrote", p)

if __name__ == "__main__":
    main()
