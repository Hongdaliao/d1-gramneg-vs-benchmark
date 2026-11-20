#!/usr/bin/env python
"""LIT-PCBA full anchor runner (Amendment A2, task-31).

A2 protocol:
  * inactives/active = 20:1, cap 1,000/set: n_inact = min(20*n_val_actives, 1000, pool)
  * 3 inactive subsample seeds (42/43/44); Vina seed 42 single (stability 3-seed only on ADRB2/VDR/ALDH1)
  * Vina 1.2.5, exhaustiveness 4, num_modes 9, box 22^3, 16x cpu1
  * 15/15 sets attempted; failures recorded
Outputs:
  results/litpcba/<set>/scores_s<subseed>_v42.csv
  results/litpcba/litpcba_scores_all.csv
  results/litpcba/anchor_metrics.csv (AC1-AC3 per set: EF5/EF1/BEDROC/AUC + subsample variance)
  results/litpcba/full_run.log
Resumable: skips any (set,subseed) whose score CSV already has rows.
Run detached/background:
  tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\_b_litpcba_full.py
"""
import csv, json, math, os, sys, time

ROOT = r"D:\harness\projects\d1_gramneg_vs"
sys.path.insert(0, os.path.join(ROOT, "code", "trackA"))
import _b_litpcba as L

OUT = os.path.join(ROOT, "results", "litpcba")
LOG = os.path.join(OUT, "full_run.log")
SUBS = [42, 43, 44]
VSEED = 42
EXH = 4

def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    try:
        open(LOG, "a", encoding="utf-8").write(line + "\n")
    except Exception:
        pass

def counts(s, subseed):
    d = os.path.join(OUT, s)
    man = os.path.join(d, f"manifest_s{subseed}.csv")
    if os.path.exists(man):
        rows = list(csv.DictReader(open(man, encoding="utf-8-sig")))
        return sum(1 for r in rows if r["label"] == "1"), sum(1 for r in rows if r["label"] == "0")
    return None, None

def run_one(s, subseed):
    d = os.path.join(OUT, s)
    sp = os.path.join(d, f"scores_s{subseed}_v{VSEED}.csv")
    if os.path.exists(sp):
        rows = list(csv.DictReader(open(sp, encoding="utf-8-sig")))
        if len(rows) > 0:
            log(f"{s} sub={subseed}: scores exist ({len(rows)} rows), skip")
            return sp
    log(f"{s} sub={subseed}: prep")
    L.cmd_prep(s, subseed, 12, 20, 1000)
    log(f"{s} sub={subseed}: dock (exh{EXH}, vseed {VSEED})")
    L.cmd_dock(s, subseed, VSEED, 16, 1, EXH)
    return sp

def metrics_from_scores(rows):
    """rows: ligand_id,label,score,status. Returns AUC/EF1/EF5/BEDROC."""
    data = []
    for r in rows:
        if r["status"] != "ok" or r.get("score") in ("", None):
            continue
        data.append((float(r["score"]), 1 if str(r["label"]) == "1" else 0, r["ligand_id"]))
    data.sort(key=lambda x: (x[0], x[2]))
    N = len(data); nact = sum(1 for x in data if x[1] == 1); nneg = N - nact
    if nact == 0 or nneg == 0:
        return {"N": N, "n_act": nact, "auc": None, "ef1": None, "ef5": None, "bedroc20": None}
    ranks = {}
    i = 0
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
    def ef(frac):
        k = max(1, int(frac * N))
        hits = sum(1 for x in data[:k] if x[1] == 1)
        return round((hits / nact) / (k / N), 4)
    # BEDROC a=20
    alpha = 20.0
    ra = nact / N
    s = sum(math.exp(-alpha * i / N) for i, x in enumerate(data) if x[1] == 1)
    denom = (1 - math.exp(-alpha)) / (math.exp(alpha / N) - 1 + 1e-12) + ra * alpha / (1 - math.exp(-alpha))
    bedroc = round(s / (nact * denom), 4) if nact and denom else None
    return {"N": N, "n_act": nact, "n_neg": nneg, "auc": round(auc, 4),
            "ef5": ef(0.05), "ef1": ef(0.01), "bedroc20": bedroc}

def main():
    log("=== LIT-PCBA full anchor start (Amendment A2) ===")
    panel = L.load_panel()
    sets = list(panel.keys())
    per = []
    for s in sets:
        rows_all = []
        for sub in SUBS:
            try:
                sp = run_one(s, sub)
                rows = list(csv.DictReader(open(sp, encoding="utf-8-sig")))
                rows_all.append((sub, rows))
            except Exception as e:
                log(f"{s} sub={sub} FAILED {type(e).__name__}: {str(e)[:200]}")
                rows_all.append((sub, []))
        # per-set metrics: per-subseed + mean/min/max/VAR of EF5
        efs = []
        for sub, rows in rows_all:
            m = metrics_from_scores(rows)
            m.update({"set": s, "subseed": sub, "vina_seed": VSEED, "exhaustiveness": EXH})
            per.append(m)
            if m.get("ef5") is not None:
                efs.append(m["ef5"])
        log(f"{s}: done subsets; EF5 values={efs}")
    # write per-subseed metrics
    p = os.path.join(OUT, "anchor_metrics_by_subseed.csv")
    if per:
        cols = ["set", "subseed", "vina_seed", "exhaustiveness", "N", "n_act", "n_neg", "auc", "ef5", "ef1", "bedroc20"]
        with open(p, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(per)
    # aggregate per set: primary = subseed 42; AC3 = CV% of EF5 across 3 subseeds
    agg = []
    for s in sets:
        r42 = [r for r in per if r["set"] == s and r["subseed"] == 42]
        efs = [r["ef5"] for r in per if r["set"] == s and r.get("ef5") is not None]
        a42 = r42[0] if r42 else {}
        mean_ef = sum(efs) / len(efs) if efs else None
        cv = None
        if efs and mean_ef:
            import statistics
            cv = round(statistics.pstdev(efs) / mean_ef * 100, 2)
        n_act = a42.get("n_act") if a42 else None
        agg.append({"set": s, "n_act_val": n_act, "N_primary": a42.get("N"), "auc": a42.get("auc"),
                    "ef5": a42.get("ef5"), "ef1": a42.get("ef1"), "bedroc20": a42.get("bedroc20"),
                    "ef5_mean_3sub": round(mean_ef, 4) if mean_ef is not None else None,
                    "ef5_cv_pct": cv, "ef5_ge2": (a42.get("ef5") is not None and a42.get("ef5") >= 2),
                    "ef1_ge2": (a42.get("ef1") is not None and a42.get("ef1") >= 2),
                    "n_act_ge150": (n_act is not None and n_act >= 150),
                    "status": "done" if a42.get("N") else "failed"})
    pa = os.path.join(OUT, "anchor_metrics.csv")
    with open(pa, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(agg[0].keys())); w.writeheader(); w.writerows(agg)
    # AC2: sets with n_act_val>=20; proportion EF5>=2
    elig = [r for r in agg if (r["n_act_val"] or 0) >= 20 and r["status"] == "done"]
    ge2 = [r for r in elig if r["ef5"] is not None and r["ef5"] >= 2]
    ac = {"AC1_attempted": len(agg), "AC1_completed": sum(1 for r in agg if r["status"] == "done"),
          "AC2_eligible_sets": len(elig), "AC2_ef5_ge2": len(ge2),
          "AC2_proportion": round(len(ge2) / len(elig), 3) if elig else None,
          "AC2_pass_ge_30pct": (len(ge2) / len(elig) >= 0.30) if elig else False,
          "AC3_ef5_cv_pct_max": max([r["ef5_cv_pct"] for r in agg if r["ef5_cv_pct"] is not None], default=None),
          "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
          "protocol": "Amendment A2: 20:1 cap1000, 3 inactive subseeds, Vina seed 42, exhaustiveness 4, box 22^3"}
    json.dump(ac, open(os.path.join(OUT, "anchor_AC.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    log("AC: " + json.dumps(ac, ensure_ascii=False))
    log("=== full anchor done ===")

if __name__ == "__main__":
    main()
