#!/usr/bin/env python
"""Supplemental docking for ligands rescued after the main run (prep_missing -> docked).

Targets rows in <target>_scores.csv with status=prep_missing whose PDBQT now exists;
docks only those (same Vina settings), then merges scores/status back into the canonical CSV.
Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\_b_supplement.py --target LpxC_Ec --nproc 12 --exhaustiveness 4
"""
import argparse, csv, json, os, re, subprocess, sys, time

ROOT = r"D:\harness\projects\d1_gramneg_vs"
OUTROOT = os.path.join(ROOT, "results", "trackA", "full")
VINA = os.path.join(ROOT, "tools", "vina", "vina.exe")
REC = {"LpxC_Ec": "LpxC_Ec_holo.pdbqt", "LpxC_Pa": "LpxC_Pa_holo.pdbqt", "DHFR_Ec": "DHFR_Ec_holo_3DAU.pdbqt",
       "GyrB_Ec": "GyrB_Ec_holo.pdbqt", "Gyrase_Ec": "GyrB_Ec_holo.pdbqt", "GyrA_Ec": "GyrA_Ec_holo.pdbqt",
       "FabI_Ec": "FabI_Ec_holo.pdbqt", "TEM1_Ec": "TEM1_Ec_holo.pdbqt", "LepB_Ec": "LepB_Ec_holo.pdbqt"}
BOX = {"LpxC_Ec": [8.621, -32.447, -15.695], "LpxC_Pa": [-4.066, 10.272, 12.367], "DHFR_Ec": [27.158, 9.763, -17.934],
       "GyrB_Ec": [-28.103, 7.970, -1.874], "Gyrase_Ec": [-28.103, 7.970, -1.874],
       "GyrA_Ec": [-7.561, 50.828, 29.374], "FabI_Ec": [12.555, 15.073, 11.276],
       "TEM1_Ec": [3.088, 13.906, 18.215], "LepB_Ec": [64.921, 22.645, 10.289]}
SIZE = {"LpxC_Ec": 26, "GyrB_Ec": 26, "Gyrase_Ec": 26, "LepB_Ec": 26}
SEED = 42

def launch(args_list, logs):
    procs = []
    for args, log in zip(args_list, logs):
        os.makedirs(os.path.dirname(log), exist_ok=True)
        lf = open(log, "w", encoding="utf-8", errors="replace")
        procs.append((subprocess.Popen([sys.executable, os.path.abspath(__file__)] + args,
                                       stdout=lf, stderr=subprocess.STDOUT, cwd=ROOT), lf))
    for p, lf in procs:
        p.wait(); lf.close()

def child(target, chunk, nchunks, exh):
    d = os.path.join(OUTROOT, target)
    scores = list(csv.DictReader(open(os.path.join(d, f"{target}_scores.csv"), encoding="utf-8-sig")))
    todo = [r for r in scores if r["status"] == "prep_missing" and os.path.exists(os.path.join(d, "ligands", r["ligand_id"] + ".pdbqt"))]
    mine = [r for i, r in enumerate(todo) if i % nchunks == chunk]
    rec = os.path.join(ROOT, "data", "receptors", REC[target]); box = BOX[target]; size = SIZE.get(target, 22)
    out = []
    for r in mine:
        lig = os.path.join(d, "ligands", r["ligand_id"] + ".pdbqt")
        outp = os.path.join(d, "poses", r["ligand_id"] + "_out.pdbqt")
        logf = os.path.join(d, "logs", "sup_" + r["ligand_id"] + ".log")
        cmd = [VINA, "--receptor", rec, "--ligand", lig, "--center_x", str(box[0]), "--center_y", str(box[1]),
               "--center_z", str(box[2]), "--size_x", str(size), "--size_y", str(size), "--size_z", str(size),
               "--exhaustiveness", str(exh), "--seed", str(SEED), "--num_modes", "3", "--cpu", "1", "--out", outp]
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
        out.append({"ligand_id": r["ligand_id"], "label": r["label"], "score": "" if score is None else score,
                    "status": "ok(supplement)" if score is not None else "vina_fail", "seed": SEED, "exhaustiveness": exh})
    with open(os.path.join(d, f"supp_part{chunk}.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "score", "status", "seed", "exhaustiveness"])
        w.writeheader(); w.writerows(out)
    print(json.dumps({"chunk": chunk, "n": len(out)}), flush=True)

def main_cmd(target, nproc, exh):
    d = os.path.join(OUTROOT, target)
    scores_p = os.path.join(d, f"{target}_scores.csv")
    scores = list(csv.DictReader(open(scores_p, encoding="utf-8-sig")))
    todo = [r for r in scores if r["status"] == "prep_missing" and os.path.exists(os.path.join(d, "ligands", r["ligand_id"] + ".pdbqt"))]
    print(json.dumps({"target": target, "n_scores": len(scores), "n_rescued_to_dock": len(todo)}), flush=True)
    if not todo:
        return
    args = [["child", "--target", target, "--chunk", str(c), "--nchunks", str(nproc), "--exhaustiveness", str(exh)] for c in range(nproc)]
    logs = [os.path.join(d, "logs", f"supp_{c}.log") for c in range(nproc)]
    t0 = time.time(); launch(args, logs)
    merged = {r["ligand_id"]: r for r in scores}
    for c in range(nproc):
        p = os.path.join(d, f"supp_part{c}.csv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p, encoding="utf-8-sig")):
                merged[r["ligand_id"]] = {**merged.get(r["ligand_id"], {}), **r}
    cols = ["ligand_id", "label", "score", "status", "seed", "exhaustiveness"]
    with open(scores_p, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows([merged[k] for k in merged])
    ok = sum(1 for r in merged.values() if str(r.get("status", "")).startswith("ok"))
    print(json.dumps({"supplemented": len(todo), "elapsed_s": round(time.time()-t0, 1), "n_ok_total": ok,
                      "n_rows": len(merged)}, indent=1), flush=True)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("run"); p1.add_argument("--target", required=True); p1.add_argument("--nproc", type=int, default=12); p1.add_argument("--exhaustiveness", type=int, default=4)
    p2 = sub.add_parser("child"); p2.add_argument("--target", required=True); p2.add_argument("--chunk", type=int, required=True); p2.add_argument("--nchunks", type=int, required=True); p2.add_argument("--exhaustiveness", type=int, default=4)
    a = ap.parse_args()
    if a.cmd == "run":
        main_cmd(a.target, a.nproc, a.exhaustiveness)
    else:
        child(a.target, a.chunk, a.nchunks, a.exhaustiveness)
