#!/usr/bin/env python
"""Box-size three-tier comparison (Amendment #11): box22 / box26 / box28 on LpxC_Ec subset.

box26 = main panel (base), box22 = informative failure (physical truncation), box28 = frozen
equivalence arm (Amendment #11, 2026-11-19). Same 1,000 ligands, same center, exh4, seed 42.
Outputs: results/trackA/sensitivity/box/three_tier.csv + box28_summary.json
Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\_b_box28_compare.py
"""
import csv, json, math, os, random, statistics

ROOT = r"D:\harness\projects\d1_gramneg_vs"
TRACKA = os.path.join(ROOT, "results", "trackA", "full")
SENS = os.path.join(ROOT, "results", "trackA", "sensitivity")
SEED = 42
BOOT = 2000

def load(p):
    return {r["ligand_id"]: r for r in csv.DictReader(open(p, encoding="utf-8-sig"))}

def ef_auc(pairs):
    data = sorted(pairs, key=lambda x: (x[0], x[2]))
    N = len(data); nact = sum(1 for x in data if x[1] == 1); nneg = N - nact
    if nact == 0 or nneg == 0:
        return None, None
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
    k = max(1, int(0.05 * N))
    hits = sum(1 for x in data[:k] if x[1] == 1)
    return round(auc, 4), round((hits / nact) / (k / N), 4)

def spearman(xs, ys):
    n = len(xs)
    def ranks(v):
        order = sorted(range(n), key=lambda i: v[i]); r = [0] * n; i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = ranks(xs), ranks(ys)
    mx = sum(rx) / n; my = sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else None

def compare(ids, base, alt):
    pb, pa = [], []
    for lid in ids:
        b = base.get(lid); a = alt.get(lid)
        if not b or not a or b.get("status") != "ok" or a.get("status") != "ok":
            continue
        if b.get("score") in ("", None) or a.get("score") in ("", None):
            continue
        lab = 1 if str(b["label"]) == "1" else 0
        pb.append((float(b["score"]), lab, lid)); pa.append((float(a["score"]), lab, lid))
    if len(pb) < 20:
        return None
    auc_b, ef_b = ef_auc(pb); auc_a, ef_a = ef_auc(pa)
    rho = spearman([p[0] for p in pb], [p[0] for p in pa])
    rng = random.Random(SEED)
    bmap = {p[2]: p[0] for p in pb}; amap = {p[2]: p[0] for p in pa}; labs = {p[2]: p[1] for p in pb}
    ids_all = list(bmap)
    d_ef = []
    for _ in range(BOOT):
        smp = [rng.choice(ids_all) for _ in ids_all]
        xb = [(bmap[l], labs[l], str(i)) for i, l in enumerate(smp)]
        xa = [(amap[l], labs[l], str(i)) for i, l in enumerate(smp)]
        _, eb = ef_auc(xb); _, ea = ef_auc(xa)
        if eb is not None and ea is not None:
            d_ef.append(ea - eb)
    d_ef.sort()
    ci = (round(d_ef[int(0.025 * len(d_ef))], 4), round(d_ef[int(0.975 * len(d_ef)) - 1], 4)) if d_ef else (None, None)
    rho_ok = rho is not None and rho >= 0.90
    dauc_ok = abs(auc_a - auc_b) <= 0.05
    ef_ci_ok = ci[0] is not None and ci[0] <= 0 <= ci[1]
    return {"n": len(pb), "auc_base": auc_b, "auc_alt": auc_a, "dAUC": round(auc_a - auc_b, 4),
            "ef5_base": ef_b, "ef5_alt": ef_a, "dEF5": round(ef_a - ef_b, 4), "dEF5_ci95": ci,
            "spearman_rho": round(rho, 4) if rho is not None else None,
            "rho_ge_0.90": rho_ok, "abs_dAUC_le_0.05": dauc_ok, "dEF5_in_CI95": ef_ci_ok,
            "verdict": "PASS" if (rho_ok and dauc_ok and ef_ci_ok) else "FAIL"}

def main():
    ids = [r["ligand_id"] for r in csv.DictReader(open(os.path.join(SENS, "subset_LpxC_Ec.csv"), encoding="utf-8-sig"))]
    base = load(os.path.join(TRACKA, "LpxC_Ec", "LpxC_Ec_scores.csv"))
    arms = {"box22": os.path.join(SENS, "box22", "LpxC_Ec", "scores.csv"),
            "box28": os.path.join(SENS, "box28", "LpxC_Ec", "scores.csv")}
    rows = []
    for name, p in arms.items():
        if not os.path.exists(p):
            rows.append({"arm": name, "status": "missing"})
            continue
        r = compare(ids, base, load(p))
        rows.append({"arm": name, **(r or {"status": "too_few"})})
    # per-arm scores vs box26 table
    p = os.path.join(SENS, "box", "three_tier.csv")
    os.makedirs(os.path.join(SENS, "box"), exist_ok=True)
    cols = ["arm", "n", "auc_base", "auc_alt", "dAUC", "ef5_base", "ef5_alt", "dEF5", "dEF5_ci95",
            "spearman_rho", "rho_ge_0.90", "abs_dAUC_le_0.05", "dEF5_in_CI95", "verdict", "status"]
    with open(p, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    box28 = next((r for r in rows if r["arm"] == "box28"), {})
    summary = {"base": "box26 (main panel)", "reference": "LpxC_Ec 1,000-ligand subset (all actives + sampled decoys, seed42)",
               "exhaustiveness": 4, "vina_seed": 42, "box22": next((r for r in rows if r["arm"] == "box22"), {}),
               "box28": box28,
               "amendment": "#11 (FROZEN 2026-11-19): box28 equivalence arm; PASS -> box26 stays main, box22 = informative failure (physical truncation)",
               "conclusion": ("box28 PASS: box26 retained as main; box22 recorded as informative failure"
                              if box28.get("verdict") == "PASS" else
                              ("box28 FAIL: report to Lead; absolute metrics flagged box-dependent (no main-panel change without Lead decision)"
                               if box28 else "box28 pending"))}
    json.dump(summary, open(os.path.join(SENS, "box", "box28_summary.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(json.dumps(summary, indent=1, ensure_ascii=False))
    print("wrote", p)

if __name__ == "__main__":
    main()
