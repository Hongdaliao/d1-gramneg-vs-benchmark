#!/usr/bin/env python
"""Recompute LIT-PCBA anchor metrics from final score files (post-run, Amendment A2).

Reads results/litpcba/<set>/scores_s{42,43,44}_v42.csv + manifests, recomputes per-set
metrics (EF5/EF1/BEDROC/AUC), AC1-AC3, writes:
  results/litpcba/anchor_metrics.csv        (primary = subseed 42; AC3 = CV% across 3 subseeds)
  results/litpcba/anchor_AC.json
Notes ADRB2 as chain-validation / tiny sample (4 actives) => descriptive only.
Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\_b_anchor_recompute.py
"""
import csv, json, math, os, statistics, sys, time

ROOT = r"D:\harness\projects\d1_gramneg_vs"
OUT = os.path.join(ROOT, "results", "litpcba")
SUBS = [42, 43, 44]
TINY = {"ADRB2", "OPRK1", "ESR1_ago", "PPARG"}  # <20 validation actives => descriptive only

def metrics(rows):
    data = [(float(r["score"]), 1 if str(r["label"]) == "1" else 0, r["ligand_id"])
            for r in rows if r["status"] == "ok" and r.get("score") not in ("", None)]
    data.sort(key=lambda x: (x[0], x[2]))
    N = len(data); nact = sum(1 for x in data if x[1] == 1); nneg = N - nact
    if nact == 0 or nneg == 0:
        return {"N": N, "n_act": nact, "auc": None, "ef5": None, "ef1": None, "bedroc20": None}
    ranks = {}; i = 0
    while i < N:
        j = i
        while j + 1 < N and data[j + 1][0] == data[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks[data[k][2]] = avg
        i = j + 1
    sr = sum(ranks[x[2]] for x in data if x[1] == 1)
    auc = 1 - ((sr - nact * (nact + 1) / 2.0) / (nact * nneg))
    def ef(f):
        k = max(1, int(f * N))
        hits = sum(1 for x in data[:k] if x[1] == 1)
        return round((hits / nact) / (k / N), 4)
    alpha = 20.0; ra = nact / N
    s = sum(math.exp(-alpha * i / N) for i, x in enumerate(data) if x[1] == 1)
    denom = (1 - math.exp(-alpha)) / (math.exp(alpha / N) - 1) + ra * alpha / (1 - math.exp(-alpha))
    return {"N": N, "n_act": nact, "n_neg": nneg, "auc": round(auc, 4), "ef5": ef(0.05),
            "ef1": ef(0.01), "bedroc20": round(s / (nact * denom), 4)}

def main():
    sets = sorted({os.path.basename(d) for d in
                   [os.path.dirname(p) for p in __import__("glob").glob(os.path.join(OUT, "*", "scores_s*_v42.csv"))]})
    per, agg = [], []
    for s in sets:
        efs = []
        for sub in SUBS:
            sp = os.path.join(OUT, s, f"scores_s{sub}_v42.csv")
            if not os.path.exists(sp):
                per.append({"set": s, "subseed": sub, "N": 0, "n_act": 0, "status": "missing"})
                continue
            rows = list(csv.DictReader(open(sp, encoding="utf-8-sig")))
            m = metrics(rows)
            m.update({"set": s, "subseed": sub, "status": "done"})
            per.append(m)
            if m.get("ef5") is not None:
                efs.append(m["ef5"])
        r42 = [r for r in per if r["set"] == s and r["subseed"] == 42 and r.get("status") == "done"]
        a42 = r42[0] if r42 else {}
        mean_ef = sum(efs) / len(efs) if efs else None
        cv = round(statistics.pstdev(efs) / mean_ef * 100, 2) if efs and mean_ef else None
        n_act = a42.get("n_act")
        agg.append({"set": s, "n_act_val": n_act, "N_primary": a42.get("N"),
                    "auc": a42.get("auc"), "ef5": a42.get("ef5"), "ef1": a42.get("ef1"),
                    "bedroc20": a42.get("bedroc20"),
                    "ef5_mean_3sub": round(mean_ef, 4) if mean_ef is not None else None,
                    "ef5_cv_pct": cv,
                    "ef5_ge2": bool(a42.get("ef5") is not None and a42["ef5"] >= 2),
                    "ef1_ge2": bool(a42.get("ef1") is not None and a42["ef1"] >= 2),
                    "n_act_ge150": bool(n_act is not None and n_act >= 150),
                    "descriptive_only": s in TINY,
                    "note": "chain validation / tiny sample (4 actives), descriptive only" if s == "ADRB2" else "",
                    "status": "done" if a42.get("N") else "failed"})
    pp = os.path.join(OUT, "anchor_metrics_by_subseed.csv")
    with open(pp, "w", newline="", encoding="utf-8-sig") as f:
        cols = ["set", "subseed", "status", "N", "n_act", "n_neg", "auc", "ef5", "ef1", "bedroc20"]
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(per)
    pa = os.path.join(OUT, "anchor_metrics.csv")
    with open(pa, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(agg[0].keys())); w.writeheader(); w.writerows(agg)
    # AC2 per Amendment A2: primary = sets with n_val>=20; proportion EF5>=2 >=30%
    elig = [r for r in agg if (r["n_act_val"] or 0) >= 20 and r["status"] == "done"]
    ge2 = [r for r in elig if r["ef5"] is not None and r["ef5"] >= 2]
    ac = {"AC1_attempted": len(agg), "AC1_completed": sum(1 for r in agg if r["status"] == "done"),
          "AC2_eligible_sets_n_val_ge20": len(elig), "AC2_ef5_ge2": len(ge2),
          "AC2_proportion": round(len(ge2) / len(elig), 3) if elig else None,
          "AC2_pass_ge_30pct": bool(elig and len(ge2) / len(elig) >= 0.30),
          "AC2_set_list_eligible": [r["set"] for r in elig],
          "AC2_sets_ef5_ge2": [r["set"] for r in ge2],
          "AC3_ef5_cv_pct_max": max([r["ef5_cv_pct"] for r in agg if r["ef5_cv_pct"] is not None], default=None),
          "AC3_pass_le_10pct": (max([r["ef5_cv_pct"] for r in agg if r["ef5_cv_pct"] is not None], default=999) <= 10),
          "ef1_reported_sets": [r["set"] for r in agg if r["n_act_ge150"]],
          "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
          "protocol": "Amendment A2: 20:1 cap1000 inactives, 3 inactive subseeds, Vina seed 42, exhaustiveness 4, box 22^3; AC2 on sets with n_val>=20",
          "notes": {"ADRB2": "chain validation / tiny sample (4 actives), descriptive only"}}
    json.dump(ac, open(os.path.join(OUT, "anchor_AC.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(json.dumps(ac, indent=1, ensure_ascii=False))
    print("wrote", pa)

if __name__ == "__main__":
    main()
