#!/usr/bin/env python
"""Parse Uni-Dock raw outputs into the frozen scores schema.

Raw layout: results/trackA/method_axis/unidock/out/<target>/<ligand_id>_out.pdbqt
  * first "REMARK VINA RESULT:" = best affinity (kcal/mol, lower is better)
  * number of MODEL records = n_models
Outputs:
  results/trackA/method_axis/unidock/scores/unidock_<target>.csv
      target,ligand_id,label,affinity,n_models,status   (status ok | failed_missing | failed_no_models)
  results/trackA/method_axis/unidock/scores/unidock_scores_806.csv (combined)
  results/trackA/method_axis/unidock/scores/unidock_na_coverage.csv
  results/trackA/method_axis/unidock/scores/parse_summary.json (sha256 + counts)
Safe to re-run as WSL jobs finish; no files overwritten outside unidock/scores.
"""
import argparse
import csv
import hashlib
import json
import os
import re
from collections import Counter

ROOT = r"D:\harness\projects\d1_gramneg_vs"
BASE = os.path.join(ROOT, "results", "trackA", "method_axis")
U = os.path.join(BASE, "unidock")
OUT = os.path.join(U, "out")
SURROGATE = os.path.join(U, "out_surrogate")
SCORES = os.path.join(U, "scores")
TARGETS = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec",
           "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]
SUBSETS = os.path.join(BASE, "subsets_diffdock.csv")
RESULT_RE = re.compile(r"^REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)")


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_pdbqt(path):
    best = None
    n_models = 0
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("MODEL"):
                n_models += 1
            elif best is None:
                m = RESULT_RE.match(line)
                if m:
                    best = float(m.group(1))
    return best, n_models


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--allow-partial", action="store_true", help="write combined even if some targets missing")
    args = ap.parse_args()
    os.makedirs(SCORES, exist_ok=True)
    subs = [r for r in csv.DictReader(open(SUBSETS, encoding="utf-8-sig"))]
    all_rows = []
    summary = {"targets": {}, "at": __import__("time").strftime("%Y-%m-%dT%H:%M:%S")}
    for t in TARGETS:
        tdir = os.path.join(OUT, t)
        rows = []
        files = {}
        if os.path.isdir(tdir):
            for f in os.listdir(tdir):
                if f.endswith("_out.pdbqt"):
                    files[f[:-len("_out.pdbqt")]] = os.path.join(tdir, f)
        for r in [x for x in subs if x["target"] == t]:
            lid = r["ligand_id"]
            p = files.get(lid)
            src = "native"
            if not p and os.path.isdir(os.path.join(SURROGATE, t)):
                sp = os.path.join(SURROGATE, t, f"{lid}_out.pdbqt")
                if os.path.exists(sp):
                    p, src = sp, "surrogate_B2C"
            if not p:
                rows.append({"target": t, "ligand_id": lid, "label": r["label"],
                             "affinity": "", "n_models": 0, "status": "failed_missing", "source": ""})
                continue
            aff, n = parse_pdbqt(p)
            if aff is None or n == 0:
                rows.append({"target": t, "ligand_id": lid, "label": r["label"],
                             "affinity": "", "n_models": n, "status": "failed_no_models", "source": src})
            elif aff > 0:
                # Amendment #14 rule: affinity > 0 = NA (unconverged); keep raw value in SI column
                rows.append({"target": t, "ligand_id": lid, "label": r["label"],
                             "affinity": f"{aff:.3f}", "n_models": n,
                             "status": "NA_positive_affinity", "source": src})
            else:
                rows.append({"target": t, "ligand_id": lid, "label": r["label"],
                             "affinity": f"{aff:.3f}", "n_models": n, "status": "ok", "source": src})
        out = os.path.join(SCORES, f"unidock_{t}.csv")
        with open(out, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=["target", "ligand_id", "label", "affinity", "n_models", "status", "source"])
            w.writeheader()
            w.writerows(rows)
        ok = sum(1 for x in rows if x["status"] == "ok")
        summary["targets"][t] = {"rows": len(rows), "ok": ok, "na": len(rows) - ok,
                                 "csv": os.path.relpath(out, ROOT), "csv_sha256": sha256(out)}
        all_rows.extend(rows)
    agg = os.path.join(SCORES, "unidock_scores_806.csv")
    with open(agg, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["target", "ligand_id", "label", "affinity", "n_models", "status", "source"])
        w.writeheader()
        w.writerows(all_rows)
    na = [r for r in all_rows if r["status"] != "ok"]
    na_csv = os.path.join(SCORES, "unidock_na_coverage.csv")
    with open(na_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["target", "ligand_id", "label", "affinity", "n_models", "status", "source"])
        w.writeheader()
        w.writerows(na)
    pos = [r for r in all_rows if r["status"] == "NA_positive_affinity"]
    pos_csv = os.path.join(SCORES, "unidock_positive_affinity.csv")
    with open(pos_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["target", "ligand_id", "label", "affinity", "n_models", "status", "source"])
        w.writeheader(); w.writerows(pos)
    from collections import Counter as _C
    summary["status_counts"] = dict(_C(r["status"] for r in all_rows))
    summary["source_counts"] = dict(_C(r.get("source", "") for r in all_rows))
    summary["combined"] = {"csv": os.path.relpath(agg, ROOT), "sha256": sha256(agg),
                           "rows": len(all_rows),
                           "ok": sum(1 for x in all_rows if x["status"] == "ok"),
                           "na": len(na), "na_positive_affinity": len(pos),
                           "na_csv": os.path.relpath(na_csv, ROOT), "na_csv_sha256": sha256(na_csv),
                           "positive_affinity_csv": os.path.relpath(pos_csv, ROOT),
                           "positive_affinity_csv_sha256": sha256(pos_csv),
                           "rule": "Amendment #14: affinity > 0 recorded as NA (unconverged); raw value kept in SI; excluded from AUC/ranking"}
    json.dump(summary, open(os.path.join(SCORES, "parse_summary.json"), "w"), indent=1, ensure_ascii=False)
    print(json.dumps({t: {k: v[k] for k in ("rows", "ok", "na")} for t, v in summary["targets"].items()}, indent=1))
    print("combined:", summary["combined"]["rows"], "ok:", summary["combined"]["ok"], "na:", summary["combined"]["na"])
    print("sha:", summary["combined"]["sha256"])


if __name__ == "__main__":
    main()
