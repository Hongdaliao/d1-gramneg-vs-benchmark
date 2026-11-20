#!/usr/bin/env python
"""Vinardo `--score_only` rescoring of Track A top poses (task-30 fallback).

Rationale (approved by Lead, D7): Uni-Dock infeasible on this Windows sandbox
(no PyPI wheel; GitHub releases ship no binaries; WSL access denied). Fallback =
Vina 1.2.5 built-in alternative scoring function **Vinardo** applied with
`--score_only` to the SAME docked top poses -> "second scoring function, not an
independent docking engine" (Methods wording), not merged into AI-vs-physics claims.

Inputs: results/trackA/full/<target>/poses/<ligand>_out.pdbqt (top pose)
        receptors data/receptors/<target>_holo.pdbqt, box per target (22 or 26 A)
Outputs: results/trackA/uni_dock/<target>/scores_vinardo.csv
         results/trackA/multi_method_comparison.csv (Vina vs Vinardo AUC/EF5/rho)
Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\_b_vinardo_rescore.py [--nproc 8]
"""
import argparse, csv, json, math, os, re, subprocess, sys, time

ROOT = r"D:\harness\projects\d1_gramneg_vs"
TRACKA = os.path.join(ROOT, "results", "trackA", "full")
OUT = os.path.join(ROOT, "results", "trackA", "uni_dock")
VINA = os.path.join(ROOT, "tools", "vina", "vina.exe")
os.makedirs(OUT, exist_ok=True)
REC = {"LpxC_Ec": "LpxC_Ec_holo.pdbqt", "LpxC_Pa": "LpxC_Pa_holo.pdbqt", "DHFR_Ec": "DHFR_Ec_holo_3DAU.pdbqt",
       "Gyrase_Ec": "GyrB_Ec_holo.pdbqt", "GyrB_Ec": "GyrB_Ec_holo.pdbqt", "GyrA_Ec": "GyrA_Ec_holo.pdbqt",
       "FabI_Ec": "FabI_Ec_holo.pdbqt", "TEM1_Ec": "TEM1_Ec_holo.pdbqt", "LepB_Ec": "LepB_Ec_holo.pdbqt"}
BOX = {"LpxC_Ec": [8.621, -32.447, -15.695], "LpxC_Pa": [-4.066, 10.272, 12.367], "DHFR_Ec": [27.158, 9.763, -17.934],
       "GyrB_Ec": [-28.103, 7.970, -1.874], "Gyrase_Ec": [-28.103, 7.970, -1.874],
       "GyrA_Ec": [-7.561, 50.828, 29.374], "FabI_Ec": [12.555, 15.073, 11.276],
       "TEM1_Ec": [3.088, 13.906, 18.215], "LepB_Ec": [64.921, 22.645, 10.289]}
SIZE = {"LpxC_Ec": 26, "GyrB_Ec": 26, "Gyrase_Ec": 26, "LepB_Ec": 26}
TARGETS = ["LpxC_Ec", "LpxC_Pa", "DHFR_Ec", "Gyrase_Ec"]

def top_model_file(pose):
    """Vina --score_only rejects multi-MODEL PDBQT -> write first MODEL to a temp file."""
    txt = open(pose, encoding="utf-8", errors="replace").read()
    out = pose + ".top.pdbqt"
    if "MODEL" in txt:
        lines = txt.splitlines()
        keep, inside = [], False
        for ln in lines:
            if ln.strip().startswith("MODEL"):
                inside = True; continue
            if ln.strip().startswith("ENDMDL"):
                break
            if inside:
                keep.append(ln)
        open(out, "w", encoding="utf-8").write("\n".join(keep) + "\n")
    else:
        open(out, "w", encoding="utf-8").write(txt)
    return out

def score_one(pose, rec, box, size):
    pose = top_model_file(pose)
    cmd = [VINA, "--receptor", rec, "--ligand", pose,
           "--center_x", str(box[0]), "--center_y", str(box[1]), "--center_z", str(box[2]),
           "--size_x", str(size), "--size_y", str(size), "--size_z", str(size),
           "--scoring", "vinardo", "--score_only"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=600)
        txt = (r.stdout or "") + (r.stderr or "")
    except subprocess.TimeoutExpired:
        return None, "timeout"
    m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?\d+(?:\.\d+)?)", txt)
    if m:
        return float(m.group(1)), "ok"
    return None, "score_fail"

def child(target, chunk, nchunks):
    d = os.path.join(TRACKA, target)
    rows = list(csv.DictReader(open(os.path.join(d, f"{target}_scores.csv"), encoding="utf-8-sig")))
    ok = [r for r in rows if r["status"] == "ok"]
    mine = [r for i, r in enumerate(ok) if i % nchunks == chunk]
    rec = os.path.join(ROOT, "data", "receptors", REC[target]); box = BOX[target]; size = SIZE.get(target, 22)
    od = os.path.join(OUT, target); os.makedirs(od, exist_ok=True)
    outp = os.path.join(od, f"part{chunk}.csv")
    done = {}
    if os.path.exists(outp):
        try:
            for r0 in csv.DictReader(open(outp, encoding="utf-8-sig")):
                if r0.get("status") == "ok":
                    done[r0["ligand_id"]] = r0
        except Exception:
            done = {}
    with open(outp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f); w.writerow(["ligand_id", "label", "score_vinardo", "status"])
        for r0 in done.values():
            w.writerow([r0["ligand_id"], r0["label"], r0["score_vinardo"], r0["status"]])
        f.flush()
        for r in mine:
            if r["ligand_id"] in done:
                continue
            pose = os.path.join(d, "poses", r["ligand_id"] + "_out.pdbqt")
            if not (os.path.exists(pose) and os.path.getsize(pose) > 100):
                w.writerow([r["ligand_id"], r["label"], "", "pose_missing"]); f.flush(); continue
            s, st = score_one(pose, rec, box, size)
            w.writerow([r["ligand_id"], r["label"], "" if s is None else s, st]); f.flush()
    print(json.dumps({"target": target, "chunk": chunk, "done": len(mine)}), flush=True)

def launch(target, nproc):
    os.makedirs(os.path.join(OUT, target), exist_ok=True)
    d = os.path.join(TRACKA, target)
    args = [["child", "--target", target, "--chunk", str(c), "--nchunks", str(nproc)] for c in range(nproc)]
    procs = []
    for c, a in enumerate(args):
        logf = os.path.join(OUT, target, f"rescore_{c}.log")
        fh = open(logf, "w")
        procs.append((subprocess.Popen([sys.executable, os.path.abspath(__file__)] + a, cwd=ROOT,
                                       stdout=fh, stderr=subprocess.STDOUT, text=True), fh))
    for pr, fh in procs:
        pr.wait(); fh.close()
    rows = []
    for c in range(nproc):
        p = os.path.join(OUT, target, f"part{c}.csv")
        if os.path.exists(p):
            rows += list(csv.DictReader(open(p, encoding="utf-8-sig")))
    sp = os.path.join(OUT, target, "scores_vinardo.csv")
    with open(sp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "score_vinardo", "status"]); w.writeheader(); w.writerows(rows)
    n_ok = sum(1 for r in rows if r["status"] == "ok")
    print(json.dumps({"target": target, "n": len(rows), "n_ok": n_ok}), flush=True)
    return rows

def metrics(pairs, higher_better):
    """pairs: [(score,label)]. Returns AUC/EF5/EF1 (lower better unless higher_better)."""
    data = [(s, l) for s, l in pairs if s is not None]
    data.sort(key=lambda x: (-x[0] if higher_better else x[0]))
    N = len(data); nact = sum(1 for s, l in data if l == 1); nneg = N - nact
    if nact == 0 or nneg == 0:
        return None
    pos = [s for s, l in data if l == 1]; neg = [s for s, l in data if l == 0]
    wins = 0.0
    for ps in pos:
        for ns in neg:
            if (ps > ns) if higher_better else (ps < ns):
                wins += 1
            elif ps == ns:
                wins += 0.5
    auc = wins / (nact * nneg)
    def ef(f):
        k = max(1, int(f * N))
        hits = sum(1 for s, l in data[:k] if l == 1)
        return (hits / nact) / (k / N)
    return {"N": N, "n_act": nact, "auc": round(auc, 4), "ef5": round(ef(0.05), 4), "ef1": round(ef(0.01), 4)}

def main(nproc):
    rows_comp = []
    for t in TARGETS:
        print(f"=== rescoring {t} with Vinardo --score_only ===", flush=True)
        vrows = launch(t, nproc)
        # merge with Vina scores
        vina = {r["ligand_id"]: r for r in csv.DictReader(open(os.path.join(TRACKA, t, f"{t}_scores.csv"), encoding="utf-8-sig"))}
        pairs_v, pairs_vin = [], []
        for r in vrows:
            if r["status"] != "ok" or r["score_vinardo"] in ("", None):
                continue
            lid = r["ligand_id"]; lab = 1 if str(r["label"]) == "1" else 0
            if vina.get(lid, {}).get("score") in ("", None):
                continue
            pairs_v.append((float(vina[lid]["score"]), lab))
            pairs_vin.append((float(r["score_vinardo"]), lab))
        m_v = metrics(pairs_v, higher_better=False)
        m_vin = metrics(pairs_vin, higher_better=False)  # Vinardo affinity: lower better
        # Spearman between scores
        try:
            from scipy.stats import spearmanr as spr
            rho = float(spr([p[0] for p in pairs_v], [p[0] for p in pairs_vin]).statistic)
        except Exception:
            rho = None
        rows_comp.append({"target": t, "method": "vina", **m_v, "spearman_vs_vina": 1.0})
        rows_comp.append({"target": t, "method": "vinardo_score_only", **m_vin, "spearman_vs_vina": round(rho, 4) if rho is not None else None})
        print(t, "vina:", m_v, "| vinardo:", m_vin, "| rho:", rho, flush=True)
    p = os.path.join(ROOT, "results", "trackA", "multi_method_comparison.csv")
    with open(p, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["target", "method", "N", "n_act", "auc", "ef5", "ef1", "spearman_vs_vina"]); w.writeheader(); w.writerows(rows_comp)
    print("wrote", p, flush=True)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("run"); p1.add_argument("--nproc", type=int, default=8)
    p2 = sub.add_parser("child"); p2.add_argument("--target", required=True); p2.add_argument("--chunk", type=int, required=True); p2.add_argument("--nchunks", type=int, required=True)
    a = ap.parse_args()
    if a.cmd == "run": main(a.nproc)
    else: child(a.target, a.chunk, a.nchunks)
