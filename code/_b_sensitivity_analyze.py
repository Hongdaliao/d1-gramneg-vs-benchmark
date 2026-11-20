#!/usr/bin/env python
"""Analyze frozen sensitivity arms #7a (exhaustiveness 4 vs 8) and #8 (box 22 vs 26 Å).

Inputs:
  #7a: main panel exh4 scores (results/trackA/full/<target>/<target>_scores.csv) restricted to
       the 1,000-ligand subset (results/trackA/sensitivity/subset_<target>.csv) vs
       exh8 arm scores (results/trackA/sensitivity/exh8/<target>/scores.csv)
  #8 : LpxC_Ec subset: main panel box26/exh4 vs box22/exh4 arm
       (results/trackA/sensitivity/box22/LpxC_Ec/scores.csv)
Criteria (Amendment #7a/#8, FROZEN): Spearman rho >= 0.90, |dAUC| <= 0.05,
EF5 difference within bootstrap 95% CI (compound-level bootstrap, 2000 resamples).
Outputs: results/trackA/sensitivity/exh/summary.json + criteria.csv
         results/trackA/sensitivity/box/summary.json + criteria.csv
Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\_b_sensitivity_analyze.py
"""
import csv, json, math, os, random, statistics, sys

ROOT = r"D:\harness\projects\d1_gramneg_vs"
TRACKA = os.path.join(ROOT, "results", "trackA", "full")
SENS = os.path.join(ROOT, "results", "trackA", "sensitivity")
SEED = 42
BOOT = 2000

def load_scores(path):
    return {r["ligand_id"]: r for r in csv.DictReader(open(path, encoding="utf-8-sig"))}

def subset_ids(target):
    p = os.path.join(SENS, f"subset_{target}.csv")
    return [r["ligand_id"] for r in csv.DictReader(open(p, encoding="utf-8-sig"))]

def ef_auc(pairs):
    """pairs: list of (score,label,active_key). lower score better."""
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
    ef5 = (hits / nact) / (k / N)
    return round(auc, 4), round(ef5, 4)

def spearman(xs, ys):
    n = len(xs)
    def ranks(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0] * n
        i = 0
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
    """base/alt: dict ligand_id -> row (score). Returns metrics + bootstrap CI of dEF5."""
    pairs_b, pairs_a = [], []
    rng = random.Random(SEED)
    for lid in ids:
        b = base.get(lid); a = alt.get(lid)
        if not b or not a:
            continue
        if b.get("score") in ("", None) or a.get("score") in ("", None):
            continue
        if b.get("status") != "ok" or a.get("status") != "ok":
            continue
        lab = 1 if str(b["label"]) == "1" else 0
        pairs_b.append((float(b["score"]), lab, lid))
        pairs_a.append((float(a["score"]), lab, lid))
    if len(pairs_b) < 20:
        return {"n": len(pairs_b), "error": "too_few"}
    auc_b, ef_b = ef_auc(pairs_b)
    auc_a, ef_a = ef_auc(pairs_a)
    rho = spearman([p[0] for p in pairs_b], [p[0] for p in pairs_a])
    # compound-level bootstrap of dEF5 (resample ligands)
    d_ef = []
    labels = {p[2]: p[1] for p in pairs_b}
    base_map = {p[2]: p[0] for p in pairs_b}
    alt_map = {p[2]: p[0] for p in pairs_a}
    ids_all = list(base_map.keys())
    for _ in range(BOOT):
        sample = [rng.choice(ids_all) for _ in ids_all]
        pb = [(base_map[l], labels[l], str(i)) for i, l in enumerate(sample)]
        pa = [(alt_map[l], labels[l], str(i)) for i, l in enumerate(sample)]
        ab, eb = ef_auc(pb); aa, ea = ef_auc(pa)
        if eb is not None and ea is not None:
            d_ef.append(ea - eb)
    d_ef.sort()
    ci = (round(d_ef[int(0.025 * len(d_ef))], 4), round(d_ef[int(0.975 * len(d_ef)) - 1], 4)) if d_ef else (None, None)
    return {"n": len(pairs_b),
            "auc_base": auc_b, "auc_alt": auc_a, "dAUC": round(auc_a - auc_b, 4),
            "ef5_base": ef_b, "ef5_alt": ef_a, "dEF5": round(ef_a - ef_b, 4),
            "dEF5_ci95": ci, "spearman_rho": round(rho, 4) if rho is not None else None}

def main():
    out = {}
    crit = []
    # #7a: exh8 arms (base = main panel exh4, same box)
    for t in ["LpxC_Ec", "DHFR_Ec"]:
        ids = subset_ids(t)
        base = load_scores(os.path.join(TRACKA, t, f"{t}_scores.csv"))
        alt = load_scores(os.path.join(SENS, "exh8", t, "scores.csv"))
        r = compare(ids, base, alt)
        rho_ok = r.get("spearman_rho") is not None and r["spearman_rho"] >= 0.90
        dauc_ok = r.get("dAUC") is not None and abs(r["dAUC"]) <= 0.05
        ci = r.get("dEF5_ci95") or (None, None)
        ef_ci_ok = (ci[0] is not None and ci[0] <= 0 <= ci[1])
        verdict = "PASS" if (rho_ok and dauc_ok and ef_ci_ok) else "FAIL"
        out[t] = {"arm": "7a exh4_vs_exh8", **r, "criteria": {"rho_ge_0.90": rho_ok, "abs_dAUC_le_0.05": dauc_ok,
                                                              "dEF5_in_CI95": ef_ci_ok, "verdict": verdict}}
        crit.append({"arm": "exh", "target": t, **out[t]["criteria"], **{k: r.get(k) for k in ("n", "auc_base", "auc_alt", "dAUC", "ef5_base", "ef5_alt", "dEF5", "dEF5_ci95", "spearman_rho")}})
        print(t, "7a", json.dumps(out[t], ensure_ascii=False), flush=True)
    os.makedirs(os.path.join(SENS, "exh"), exist_ok=True)
    json.dump(out, open(os.path.join(SENS, "exh", "summary.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    with open(os.path.join(SENS, "exh", "criteria.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(crit[0].keys())); w.writeheader(); w.writerows(crit)
    # #8: box22 vs box26 for LpxC_Ec subset
    out2 = {}
    crit2 = []
    ids = subset_ids("LpxC_Ec")
    base = load_scores(os.path.join(TRACKA, "LpxC_Ec", "LpxC_Ec_scores.csv"))
    alt = load_scores(os.path.join(SENS, "box22", "LpxC_Ec", "scores.csv"))
    r = compare(ids, base, alt)
    rho_ok = r.get("spearman_rho") is not None and r["spearman_rho"] >= 0.90
    dauc_ok = r.get("dAUC") is not None and abs(r["dAUC"]) <= 0.05
    ci = r.get("dEF5_ci95") or (None, None)
    ef_ci_ok = (ci[0] is not None and ci[0] <= 0 <= ci[1])
    verdict = "PASS" if (rho_ok and dauc_ok and ef_ci_ok) else "FAIL"
    out2["LpxC_Ec"] = {"arm": "8 box22_vs_box26_exh4", **r, "criteria": {"rho_ge_0.90": rho_ok, "abs_dAUC_le_0.05": dauc_ok,
                                                                        "dEF5_in_CI95": ef_ci_ok, "verdict": verdict}}
    crit2.append({"arm": "box", "target": "LpxC_Ec", **out2["LpxC_Ec"]["criteria"], **{k: r.get(k) for k in ("n", "auc_base", "auc_alt", "dAUC", "ef5_base", "ef5_alt", "dEF5", "dEF5_ci95", "spearman_rho")}})
    print("box", json.dumps(out2["LpxC_Ec"], ensure_ascii=False), flush=True)
    os.makedirs(os.path.join(SENS, "box"), exist_ok=True)
    json.dump(out2, open(os.path.join(SENS, "box", "summary.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    with open(os.path.join(SENS, "box", "criteria.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(crit2[0].keys())); w.writeheader(); w.writerows(crit2)
    print("done")

if __name__ == "__main__":
    main()
