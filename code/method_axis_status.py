#!/usr/bin/env python
"""Programmatic status/coverage summary for Amendment #12 method-axis runs.

Reads scores.csv files only (source of truth): no hand counting.
Usage: python code\\trackA\\method_axis_status.py [--na-list]
"""
import argparse
import csv
import json
import os
from collections import Counter, OrderedDict

ROOT = r"D:\harness\projects\d1_gramneg_vs"
BASE = os.path.join(ROOT, "results", "trackA", "method_axis")
TARGETS = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec",
           "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]


def rows(path):
    if not os.path.exists(path):
        return []
    return list(csv.DictReader(open(path, encoding="utf-8-sig")))


def summarize_arm(arm):
    out = OrderedDict()
    tot = Counter()
    for t in TARGETS:
        rs = rows(os.path.join(BASE, "diffdock", arm, t, "scores.csv"))
        ok = sum(1 for r in rs if r.get("status") == "ok")
        na = len(rs) - ok
        dup = [k for k, v in Counter(r["complex_name"] for r in rs).items() if v > 1]
        done = os.path.exists(os.path.join(BASE, "diffdock", arm, t, "DONE.json"))
        out[t] = {"attempted": len(rs), "ok": ok, "na": na, "done": done, "dups": dup}
        tot["attempted"] += len(rs)
        tot["ok"] += ok
        tot["na"] += na
    return out, tot


def summarize_boltz():
    rs = rows(os.path.join(BASE, "boltz2_scores.csv"))
    out = OrderedDict()
    for t in TARGETS:
        sub = [r for r in rs if r.get("target") == t]
        ok = sum(1 for r in sub if r.get("status") == "ok")
        dup = [k for k, v in Counter(r["complex_name"] for r in sub).items() if v > 1]
        out[t] = {"attempted": len(sub), "ok": ok, "na": len(sub) - ok, "dups": dup}
    tot = {"attempted": len(rs), "ok": sum(1 for r in rs if r.get("status") == "ok"),
           "na": sum(1 for r in rs if r.get("status") != "ok")}
    return out, tot


def na_rows(arm):
    out = []
    for t in TARGETS:
        p = os.path.join(BASE, "diffdock", arm, t, "scores.csv")
        done_p = os.path.join(BASE, "diffdock", arm, t, "DONE.json")
        fin = ""
        if os.path.exists(done_p):
            try:
                fin = json.load(open(done_p, encoding="utf-8")).get("finished_at", "")
            except Exception:
                pass
        for r in rows(p):
            if r.get("status") != "ok":
                out.append({"arm": arm, "target": t, "complex_name": r.get("complex_name", ""),
                            "ligand_id": r.get("ligand_id", ""), "label": r.get("label", ""),
                            "attempts": int(r.get("retry_used") or 0) + 1,
                            "failed_at": fin})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--na-list", action="store_true")
    args = ap.parse_args()
    print("=== DiffDock main ===")
    s, tot = summarize_arm("main")
    for t, d in s.items():
        print(f"  {t:10s} attempted={d['attempted']:3d} ok={d['ok']:3d} na={d['na']:2d} "
              f"done={str(d['done']):5s} dups={len(d['dups'])}")
    print(f"  TOTAL attempted={tot['attempted']} ok={tot['ok']} na={tot['na']} "
          f"coverage={100.0*tot['ok']/max(1,tot['attempted']):.1f}%")
    print("=== DiffDock singlechain ===")
    s2, tot2 = summarize_arm("singlechain")
    for t, d in s2.items():
        if d["attempted"]:
            print(f"  {t:10s} attempted={d['attempted']:3d} ok={d['ok']:3d} na={d['na']:2d} "
                  f"done={str(d['done']):5s} dups={len(d['dups'])}")
    print(f"  TOTAL attempted={tot2['attempted']} ok={tot2['ok']} na={tot2['na']}")
    print("=== Boltz-2 ===")
    s3, tot3 = summarize_boltz()
    for t, d in s3.items():
        if d["attempted"]:
            print(f"  {t:10s} attempted={d['attempted']:3d} ok={d['ok']:3d} na={d['na']:2d} dups={len(d['dups'])}")
    print(f"  TOTAL attempted={tot3['attempted']} ok={tot3['ok']} na={tot3['na']}")
    if args.na_list:
        print("=== NA rows (DiffDock main) ===")
        for r in na_rows("main"):
            print("  {target:10s} {complex_name:55s} attempts={attempts} failed_at={failed_at}".format(**r))
        print("=== NA rows (DiffDock singlechain) ===")
        for r in na_rows("singlechain"):
            print("  {target:10s} {complex_name:55s} attempts={attempts} failed_at={failed_at}".format(**r))


if __name__ == "__main__":
    main()
