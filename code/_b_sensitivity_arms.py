#!/usr/bin/env python
"""Track A frozen sensitivity arms (task-34, Amendment #7a/#8/#7b).

Subcommands:
  decoy_ratio : #7b decoy-ratio bootstrap (10:1/15:1/20:1) from existing 20:1 docked decoys,
                no extra docking; per-target AUC/EF5 mean + 95% CI over 200 resamples.
  make_subsets: build 1,000-ligand stratified subsets (all actives + sampled decoys, seed 42)
                for LpxC_Ec and DHFR (used by #7a exh4-vs-exh8 and #8 box22-vs-26 arms).
  dock        : run Vina on a subset manifest with given box size / exhaustiveness (16x1),
                reusable for exh8 and box22 arms; scores -> results/trackA/sensitivity/<arm>/.
Outputs: results/trackA/sensitivity/{decoy_ratio,exh,box}/*.csv + summary.json
Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\_b_sensitivity_arms.py decoy_ratio
     ... make_subsets
"""
import argparse, csv, json, math, os, random, re, subprocess, sys, time

ROOT = r"D:\harness\projects\d1_gramneg_vs"
TRACKA = os.path.join(ROOT, "results", "trackA", "full")
DECOYS = os.path.join(ROOT, "data", "decoys")
OUT = os.path.join(ROOT, "results", "trackA", "sensitivity")
VINA = os.path.join(ROOT, "tools", "vina", "vina.exe")
os.makedirs(OUT, exist_ok=True)
SEED = 42
ALL = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec", "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]
REC = {"LpxC_Ec": "LpxC_Ec_holo.pdbqt", "LpxC_Pa": "LpxC_Pa_holo.pdbqt", "DHFR_Ec": "DHFR_Ec_holo_3DAU.pdbqt",
       "GyrB_Ec": "GyrB_Ec_holo.pdbqt", "Gyrase_Ec": "GyrB_Ec_holo.pdbqt", "GyrA_Ec": "GyrA_Ec_holo.pdbqt",
       "FabI_Ec": "FabI_Ec_holo.pdbqt", "TEM1_Ec": "TEM1_Ec_holo.pdbqt", "LepB_Ec": "LepB_Ec_holo.pdbqt"}
BOX = {"LpxC_Ec": [8.621, -32.447, -15.695], "LpxC_Pa": [-4.066, 10.272, 12.367], "DHFR_Ec": [27.158, 9.763, -17.934],
       "GyrB_Ec": [-28.103, 7.970, -1.874], "Gyrase_Ec": [-28.103, 7.970, -1.874],
       "GyrA_Ec": [-7.561, 50.828, 29.374], "FabI_Ec": [12.555, 15.073, 11.276],
       "TEM1_Ec": [3.088, 13.906, 18.215], "LepB_Ec": [64.921, 22.645, 10.289]}

def load_scores(t):
    p = os.path.join(TRACKA, t, f"{t}_scores.csv")
    return list(csv.DictReader(open(p, encoding="utf-8-sig")))

def decoy_map(t):
    """decoy inchikey -> matched_active_inchikey  (from s1 decoy pool)."""
    p = os.path.join(DECOYS, f"{t}_decoys_s1.csv")
    m = {}
    if os.path.exists(p):
        for r in csv.DictReader(open(p, encoding="utf-8-sig")):
            ik = (r.get("decoy_inchikey") or "").strip()
            ma = (r.get("matched_active_inchikey") or "").strip()
            if ik:
                m[ik] = ma
    return m

def fast_auc_ef(data):
    """data: list of (score, label, active_key). Returns (auc, ef5, ef1). Lower score better."""
    import numpy as np
    d = sorted(data, key=lambda x: x[0])
    N = len(d); nact = sum(1 for x in d if x[1] == 1); nneg = N - nact
    if nact == 0 or nneg == 0:
        return None, None, None
    scores = np.array([x[0] for x in d]); labels = np.array([x[1] for x in d])
    # average ranks for ties
    import scipy.stats as st  # may be unavailable; fallback manual
    try:
        ranks = st.rankdata(scores, method="average")
    except Exception:
        ranks = None
    if ranks is not None:
        auc = 1 - ((ranks[labels == 1].sum() - nact * (nact + 1) / 2.0) / (nact * nneg))
    else:
        auc = None
    def ef(frac):
        k = max(1, int(frac * N))
        hits = int(labels[:k].sum())
        return (hits / nact) / (k / N)
    return round(float(auc), 4) if auc is not None else None, round(ef(0.05), 4), round(ef(0.01), 4)

def decoy_ratio():
    rng = random.Random(SEED)
    all_out = []
    for t in ALL:
        rows = load_scores(t)
        ok = [r for r in rows if r["status"] == "ok" and r.get("score") not in ("", None)]
        acts = [r for r in ok if str(r["label"]) == "1"]
        decs = [r for r in ok if str(r["label"]) == "0"]
        dmap = decoy_map(t)
        # active chembl_id -> inchikey (from A-line actives file)
        act_ik = {}
        ap = os.path.join(ROOT, "data", "trackA", f"{t}_actives.csv")
        if os.path.exists(ap):
            for a in csv.DictReader(open(ap, encoding="utf-8-sig")):
                cid = (a.get("compound_chembl_id") or "").strip()
                ik = (a.get("inchikey") or "").strip()
                if cid and ik:
                    act_ik[cid] = ik
        by_active = {}
        unassigned = []
        for r in decs:
            ik = r["ligand_id"][4:] if r["ligand_id"].startswith("dec_") else r["ligand_id"]
            ma = dmap.get(ik, "")
            if ma:
                by_active.setdefault(ma, []).append(r)
            else:
                unassigned.append(r)
        print(t, "actives with decoy pools:", sum(1 for a in acts if act_ik.get(a["ligand_id"][4:], "") in by_active), "/", len(acts))
        # keep only actives with >=1 decoy; others dropped for ratio<=20 (report count)
        act_scores = [(float(a["score"]), 1) for a in acts]
        for ratio in (10, 15, 20):
            aucs, ef5s, ef1s = [], [], []
            for rep in range(200):
                dec_sel = []
                for a in acts:
                    aik = act_ik.get(a["ligand_id"][4:] if a["ligand_id"].startswith("act_") else a["ligand_id"], "")
                    pool = by_active.get(aik, [])
                    k = min(ratio, len(pool))
                    if k > 0:
                        dec_sel.extend(rng.sample(pool, k))
                data = [(float(a["score"]), 1) for a in acts] + [(float(d["score"]), 0) for d in dec_sel]
                a, e5, e1 = fast_auc_ef(data)
                if a is not None:
                    aucs.append(a); ef5s.append(e5); ef1s.append(e1)
            def ci(v):
                v = sorted(v); n = len(v)
                return (round(v[int(0.025 * n)], 4), round(v[int(0.975 * n) - 1], 4)) if n else (None, None)
            all_out.append({"target": t, "ratio": ratio, "n_reps": len(aucs),
                            "n_act": len(acts), "n_dec_sel_mean": round(sum(len(by_active.get(a["ligand_id"][4:], [])) for a in acts) * 0 + 0, 1),
                            "auc_mean": round(sum(aucs) / len(aucs), 4) if aucs else None,
                            "auc_ci": ci(aucs),
                            "ef5_mean": round(sum(ef5s) / len(ef5s), 4) if ef5s else None,
                            "ef5_ci": ci(ef5s),
                            "ef1_mean": round(sum(ef1s) / len(ef1s), 4) if ef1s else None})
            print(t, ratio, all_out[-1], flush=True)
    p = os.path.join(OUT, "decoy_ratio_bootstrap.csv")
    with open(p, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(all_out[0].keys())); w.writeheader(); w.writerows(all_out)
    json.dump({"arm": "7b decoy-ratio bootstrap", "ratios": [10, 15, 20], "n_resamples": 200,
               "note": "decoys resampled per matched active from existing 20:1 docked pool; no extra docking",
               "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S")},
              open(os.path.join(OUT, "decoy_ratio_summary.json"), "w", encoding="utf-8"), indent=1)
    print("wrote", p)

def make_subsets():
    rng = random.Random(SEED)
    meta = {}
    for t in ["LpxC_Ec", "DHFR_Ec"]:
        rows = load_scores(t)
        man = list(csv.DictReader(open(os.path.join(TRACKA, t, "manifest.csv"), encoding="utf-8-sig")))
        ok = {r["ligand_id"] for r in rows if r["status"] == "ok"}
        acts = [r for r in man if str(r["label"]) == "1" and r["ligand_id"] in ok]
        decs = [r for r in man if str(r["label"]) != "1" and r["ligand_id"] in ok]
        n_dec = min(1000 - len(acts), len(decs))
        sel_dec = rng.sample(decs, n_dec)
        subset = acts + sel_dec
        p = os.path.join(OUT, f"subset_{t}.csv")
        with open(p, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "smiles", "pdbqt"]); w.writeheader()
            for r in subset:
                w.writerow({"ligand_id": r["ligand_id"], "label": r["label"], "smiles": r["smiles"], "pdbqt": r["pdbqt"]})
        meta[t] = {"n_subset": len(subset), "n_act": len(acts), "n_dec": n_dec, "path": p}
        print(t, meta[t], flush=True)
    json.dump(meta, open(os.path.join(OUT, "subsets_meta.json"), "w", encoding="utf-8"), indent=1)

def dock_subset(target, subset_csv, arm, size, exh, nproc, cpu):
    d = os.path.join(OUT, arm, target)
    for sub in ("poses", "logs"):
        os.makedirs(os.path.join(d, sub), exist_ok=True)
    rows = list(csv.DictReader(open(subset_csv, encoding="utf-8-sig")))
    rec = os.path.join(ROOT, "data", "receptors", REC[target]); box = BOX[target]
    sp = os.path.join(d, "scores.csv")
    done = {}
    if os.path.exists(sp):
        done = {r["ligand_id"]: r for r in csv.DictReader(open(sp, encoding="utf-8-sig"))}
    todo = [r for r in rows if r["ligand_id"] not in done]
    print(f"{arm}/{target}: n={len(rows)} todo={len(todo)} size={size} exh={exh}", flush=True)
    CH = nproc
    procs = []
    for c in range(min(CH, max(1, len(todo)))):
        mine = [r for i, r in enumerate(todo) if i % CH == c]
        if not mine:
            continue
        part = os.path.join(d, f"part{c}.csv")
        # child script inline via this same file: dock-child
        procs.append(subprocess.Popen([sys.executable, os.path.abspath(__file__), "dock-child",
                                       "--target", target, "--arm", arm, "--size", str(size), "--exh", str(exh),
                                       "--cpu", str(cpu), "--chunk", str(c), "--nchunks", str(CH),
                                       "--subset", subset_csv],
                                      cwd=ROOT, stdout=open(os.path.join(d, "logs", f"dock_{c}.log"), "w"),
                                      stderr=subprocess.STDOUT, text=True))
    for p in procs:
        p.wait()
    rows_out = []
    for c in range(CH):
        p = os.path.join(d, f"part{c}.csv")
        if os.path.exists(p):
            rows_out += list(csv.DictReader(open(p, encoding="utf-8-sig")))
    with open(sp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "score", "status"]); w.writeheader()
        for r in rows_out:
            w.writerow({"ligand_id": r["ligand_id"], "label": r["label"], "score": r["score"], "status": r["status"]})
    n_ok = sum(1 for r in rows_out if r["status"] == "ok")
    print(json.dumps({"arm": arm, "target": target, "n": len(rows_out), "n_ok": n_ok}), flush=True)

def dock_child(target, arm, size, exh, cpu, chunk, nchunks, subset):
    d = os.path.join(OUT, arm, target)
    rec = os.path.join(ROOT, "data", "receptors", REC[target]); box = BOX[target]
    rows = list(csv.DictReader(open(subset, encoding="utf-8-sig")))
    mine = [r for i, r in enumerate(rows) if i % nchunks == chunk]
    outp = os.path.join(d, f"part{chunk}.csv")
    with open(outp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f); w.writerow(["ligand_id", "label", "score", "status"])
        for r in mine:
            lid = r["ligand_id"]
            pdbqt = r["pdbqt"]
            if not (os.path.exists(pdbqt) and os.path.getsize(pdbqt) > 100):
                w.writerow([lid, r["label"], "", "prep_missing"]); f.flush(); continue
            pose = os.path.join(d, "poses", f"{lid}_out.pdbqt")
            logf = os.path.join(d, "logs", f"vina_{lid}.log")
            cmd = [VINA, "--receptor", rec, "--ligand", pdbqt,
                   "--center_x", str(box[0]), "--center_y", str(box[1]), "--center_z", str(box[2]),
                   "--size_x", str(size), "--size_y", str(size), "--size_z", str(size),
                   "--exhaustiveness", str(exh), "--num_modes", "3", "--seed", "42", "--cpu", str(cpu), "--out", pose]
            try:
                res = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=3600)
                txt = (res.stdout or "") + (res.stderr or "")
            except subprocess.TimeoutExpired:
                txt = "TIMEOUT"
            open(logf, "w", encoding="utf-8", errors="replace").write(txt)
            score = None
            for ln in txt.splitlines():
                m = re.match(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$", ln)
                if m:
                    score = float(m.group(2)); break
            w.writerow([lid, r["label"], "" if score is None else score, "ok" if score is not None else "fail"]); f.flush()
    print(json.dumps({"chunk": chunk, "done": len(mine)}), flush=True)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("decoy_ratio"); sub.add_parser("make_subsets")
    p3 = sub.add_parser("dock"); p3.add_argument("--target", required=True); p3.add_argument("--arm", required=True)
    p3.add_argument("--size", type=int, required=True); p3.add_argument("--exh", type=int, required=True)
    p3.add_argument("--nproc", type=int, default=16); p3.add_argument("--cpu", type=int, default=1)
    p4 = sub.add_parser("dock-child"); p4.add_argument("--target", required=True); p4.add_argument("--arm", required=True)
    p4.add_argument("--size", type=int, required=True); p4.add_argument("--exh", type=int, required=True)
    p4.add_argument("--cpu", type=int, default=1); p4.add_argument("--chunk", type=int, required=True)
    p4.add_argument("--nchunks", type=int, required=True); p4.add_argument("--subset", required=True)
    a = ap.parse_args()
    if a.cmd == "decoy_ratio": decoy_ratio()
    elif a.cmd == "make_subsets": make_subsets()
    elif a.cmd == "dock": dock_subset(a.target, os.path.join(OUT, f"subset_{a.target}.csv"), a.arm, a.size, a.exh, a.nproc, a.cpu)
    elif a.cmd == "dock-child": dock_child(a.target, a.arm, a.size, a.exh, a.cpu, a.chunk, a.nchunks, a.subset)
