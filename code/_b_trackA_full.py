#!/usr/bin/env python
"""Track A full Vina docking driver (B-line, W3).

Subcommands
-----------
manifest : build per-target ligand manifest (actives all + s1 decoys all)
           results/trackA/full/<target>/manifest.csv (ligand_id,label,smiles,pdbqt)
prep     : parallel SMILES->3D(SDF)->PDBQT (meeko), N child processes (no mp pipes)
test     : throughput test on first N ligands: strategy A (16 workers x cpu1) vs B (8 x cpu2)
dock     : full docking with chosen strategy; writes <target>_scores.csv + summary.json + progress
report   : aggregate per-target progress into results/trackA/full/progress.csv

Design notes
------------
* No multiprocessing (sandbox blocks named pipes): parallelism = subprocess.Popen children,
  each child is this same script with `--child-*` args, stdout/stderr -> log files, parent waits.
* Vina: tools/vina/vina.exe 1.2.5, num_modes=3, seed=42, box 22^3 (A-line handoff).
* Failure rows kept with status; resume-safe (existing pdbqt / existing scores reused).
"""
import argparse, csv, json, os, re, subprocess, sys, time, datetime

ROOT = r"D:\harness\projects\d1_gramneg_vs"
OUTROOT = os.path.join(ROOT, "results", "trackA", "full")
VINA = os.path.join(ROOT, "tools", "vina", "vina.exe")
MEEPREP = os.path.join(ROOT, "tools", "envs", "boltz", "Scripts", "mk_prepare_ligand.exe")
ACT_DIR = os.path.join(ROOT, "data", "trackA")
DEC_DIR = os.path.join(ROOT, "data", "decoys")
REC_DIR = os.path.join(ROOT, "data", "receptors")
SEED = 42

# box centers per A-line authoritative handoff (notes/pilot_handoff.md + docking handoff)
BOXES = {
    "LpxC_Ec": [8.621, -32.447, -15.695],
    "LpxC_Pa": [-4.066, 10.272, 12.367],
    "DHFR_Ec": [27.158, 9.763, -17.934],
}
RECEPTORS = {  # target -> holo pdbqt (Zn-retained for LpxC!)
    "LpxC_Ec": "LpxC_Ec_holo.pdbqt",
    "LpxC_Pa": "LpxC_Pa_holo.pdbqt",
    "DHFR_Ec": "DHFR_Ec_holo_3DAU.pdbqt",
}
ACTIVE_FILES = {
    "LpxC_Ec": "LpxC_Ec_actives.csv", "LpxC_Pa": "LpxC_Pa_actives.csv", "DHFR_Ec": "DHFR_Ec_actives.csv",
    "Gyrase_Ec": "Gyrase_Ec_actives.csv", "GyrA_Ec": "GyrA_Ec_actives.csv", "GyrB_Ec": "GyrB_Ec_actives.csv",
    "FabI_Ec": "FabI_Ec_actives.csv", "TEM1_Ec": "TEM1_Ec_actives.csv", "LepB_Ec": "LepB_Ec_actives.csv",
}
DECOY_FILES = {
    "LpxC_Ec": "LpxC_Ec_decoys_s1.csv", "LpxC_Pa": "LpxC_Pa_decoys_s1.csv", "DHFR_Ec": "DHFR_Ec_decoys_s1.csv",
    "Gyrase_Ec": "Gyrase_Ec_decoys_s1.csv", "GyrA_Ec": "GyrA_Ec_decoys_s1.csv", "GyrB_Ec": "GyrB_Ec_decoys_s1.csv",
    "FabI_Ec": "FabI_Ec_decoys_s1.csv", "TEM1_Ec": "TEM1_Ec_decoys_s1.csv", "LepB_Ec": "LepB_Ec_decoys_s1.csv",
}

def now():
    return datetime.datetime.now().isoformat(timespec="seconds")

def tdir(target):
    d = os.path.join(OUTROOT, target)
    os.makedirs(os.path.join(d, "ligands"), exist_ok=True)
    os.makedirs(os.path.join(d, "poses"), exist_ok=True)
    os.makedirs(os.path.join(d, "logs"), exist_ok=True)
    return d

def read_actives(target):
    p = os.path.join(ACT_DIR, ACTIVE_FILES[target])
    out = []
    for r in csv.DictReader(open(p, encoding="utf-8-sig")):
        smi = (r.get("canonical_smiles") or "").strip()
        if not smi:
            continue
        cid = (r.get("compound_chembl_id") or r.get("inchikey") or "").strip() or f"a{len(out)}"
        out.append({"ligand_id": "act_" + cid, "label": 1, "smiles": smi})
    return out

def read_decoys(target):
    p = os.path.join(DEC_DIR, DECOY_FILES[target])
    out = []
    for r in csv.DictReader(open(p, encoding="utf-8-sig")):
        smi = (r.get("decoy_canonical_smiles") or "").strip()
        if not smi:
            continue
        cid = (r.get("decoy_chembl_id") or r.get("decoy_inchikey") or "").strip() or f"d{len(out)}"
        out.append({"ligand_id": "dec_" + cid, "label": 0, "smiles": smi})
    return out

def build_manifest(target, limit=None):
    d = tdir(target)
    acts = read_actives(target)
    decs = read_decoys(target)
    rows = acts + decs
    if limit:
        na = max(1, int(limit * len(acts) / max(1, len(acts) + len(decs))))
        # keep actives proportion, sample decoys
        rows = acts[:max(1, min(len(acts), na))] + decs[:max(1, limit - min(len(acts), na))]
    for r in rows:
        r["pdbqt"] = os.path.join(d, "ligands", r["ligand_id"] + ".pdbqt")
        r["status"] = "pending"
    mpath = os.path.join(d, "manifest.csv" if not limit else f"manifest_sub{limit}.csv")
    with open(mpath, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "smiles", "pdbqt", "status"])
        w.writeheader(); w.writerows(rows)
    print(json.dumps({"target": target, "n_act": len(acts), "n_dec": len(decs), "n_rows": len(rows),
                      "manifest": mpath}, indent=1), flush=True)
    return mpath

def prep_child(manifest, chunk, nchunks):
    rows = list(csv.DictReader(open(manifest, encoding="utf-8-sig")))
    mine = [r for i, r in enumerate(rows) if i % nchunks == chunk]
    from rdkit import Chem
    from rdkit.Chem import AllChem
    ok = fail = skip = 0
    for r in mine:
        if os.path.exists(r["pdbqt"]) and os.path.getsize(r["pdbqt"]) > 100:
            skip += 1; continue
        try:
            m = Chem.MolFromSmiles(r["smiles"])
            if m is None:
                fail += 1; continue
            m = Chem.AddHs(m)
            p = AllChem.ETKDGv3(); p.randomSeed = SEED
            if AllChem.EmbedMolecule(m, p) != 0:
                AllChem.EmbedMolecule(m, useRandomCoords=True, randomSeed=SEED)
            try:
                AllChem.MMFFOptimizeMolecule(m, maxIters=200)
            except Exception:
                pass
            sdf = r["pdbqt"].replace(".pdbqt", ".sdf")
            Chem.MolToMolFile(m, sdf)
            res = subprocess.run([MEEPREP, "-i", sdf, "-o", r["pdbqt"]], capture_output=True,
                                 text=True, encoding="utf-8", errors="replace", timeout=300)
            if res.returncode == 0 and os.path.exists(r["pdbqt"]):
                ok += 1
            else:
                fail += 1
        except Exception:
            fail += 1
    print(json.dumps({"chunk": chunk, "prep_ok": ok, "skip": skip, "fail": fail}), flush=True)

def launch_children(args_list, tag):
    procs = []
    d = tdir(tag) if False else None
    for i, args in enumerate(args_list):
        log = args.pop("_log")
        os.makedirs(os.path.dirname(log), exist_ok=True) if os.path.dirname(log) else None
        with open(log, "w", encoding="utf-8", errors="replace") as lf:
            p = subprocess.Popen([sys.executable, os.path.abspath(__file__)] + args,
                                 stdout=lf, stderr=subprocess.STDOUT, cwd=ROOT)
            procs.append(p)
    for p in procs:
        p.wait()

def cmd_prep(target, manifest, nproc):
    d = tdir(target)
    args = []
    for c in range(nproc):
        args.append({"args": ["prep-child", "--manifest", manifest, "--chunk", str(c), "--nchunks", str(nproc)],
                     "_log": os.path.join(d, "logs", f"prep_{c}.log")})
    t0 = time.time()
    _launch([a["args"] for a in args], [a["_log"] for a in args])
    rows = list(csv.DictReader(open(manifest, encoding="utf-8-sig")))
    n_ok = sum(1 for r in rows if os.path.exists(r["pdbqt"]) and os.path.getsize(r["pdbqt"]) > 100)
    print(json.dumps({"prep_done": n_ok, "n": len(rows), "elapsed_s": round(time.time()-t0, 1)}), flush=True)
    return n_ok

def _launch(args_list, logs):
    procs = []
    for args, log in zip(args_list, logs):
        os.makedirs(os.path.dirname(log), exist_ok=True)
        lf = open(log, "w", encoding="utf-8", errors="replace")
        procs.append((subprocess.Popen([sys.executable, os.path.abspath(__file__)] + args,
                                       stdout=lf, stderr=subprocess.STDOUT, cwd=ROOT), lf))
    for p, lf in procs:
        p.wait(); lf.close()

def vina_one(rec, lig, box, cpu, out_pdbqt, logf, exh=8, seed=SEED):
    if os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 500:
        try:
            txt = open(logf, encoding="utf-8", errors="replace").read()
            for ln in txt.splitlines():
                m = re.match(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$", ln)
                if m:
                    return float(m.group(2)), "ok(reuse)"
        except Exception:
            pass
    cmd = [VINA, "--receptor", rec, "--ligand", lig, "--center_x", str(box[0]), "--center_y", str(box[1]),
           "--center_z", str(box[2]), "--size_x", "22", "--size_y", "22", "--size_z", "22",
           "--exhaustiveness", str(exh), "--seed", str(seed), "--num_modes", "3", "--cpu", str(cpu), "--out", out_pdbqt]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=3600)
        txt = (res.stdout or "") + (res.stderr or "")
    except subprocess.TimeoutExpired:
        return None, "timeout"
    open(logf, "w", encoding="utf-8", errors="replace").write(txt)
    for ln in txt.splitlines():
        m = re.match(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$", ln)
        if m:
            return float(m.group(2)), "ok"
    return None, "fail"

def dock_child(target, manifest, chunk, nchunks, cpu, box_csv):
    d = tdir(target)
    rec = os.path.join(REC_DIR, RECEPTORS[target])
    box = [float(x) for x in box_csv.split(",")]
    rows = list(csv.DictReader(open(manifest, encoding="utf-8-sig")))
    mine = [r for i, r in enumerate(rows) if i % nchunks == chunk]
    outp = os.path.join(d, f"scores_part{chunk}.csv")
    done = 0
    with open(outp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f); w.writerow(["ligand_id", "label", "score", "status", "seed"])
        for r in mine:
            if not (os.path.exists(r["pdbqt"]) and os.path.getsize(r["pdbqt"]) > 100):
                w.writerow([r["ligand_id"], r["label"], "", "prep_missing", SEED]); f.flush(); continue
            out_pdbqt = os.path.join(d, "poses", r["ligand_id"] + "_out.pdbqt")
            logf = os.path.join(d, "logs", "vina_" + r["ligand_id"] + ".log")
            score, status = vina_one(rec, r["pdbqt"], box, cpu, out_pdbqt, logf)
            w.writerow([r["ligand_id"], r["label"], "" if score is None else score, status, SEED]); f.flush()
            done += 1
    print(json.dumps({"chunk": chunk, "docked": done}), flush=True)

def cmd_dock(target, manifest, nproc, cpu):
    d = tdir(target)
    box = BOXES[target]
    args, logs = [], []
    for c in range(nproc):
        args.append(["dock-child", "--target", target, "--manifest", manifest, "--chunk", str(c),
                     "--nchunks", str(nproc), "--cpu", str(cpu), "--box", ",".join(map(str, box))])
        logs.append(os.path.join(d, "logs", f"dock_{c}.log"))
    t0 = time.time()
    _launch(args, logs)
    # aggregate
    parts = [os.path.join(d, f"scores_part{c}.csv") for c in range(nproc)]
    rows = []
    for p in parts:
        if os.path.exists(p):
            rows += list(csv.DictReader(open(p, encoding="utf-8-sig")))
    out = os.path.join(d, f"{target}_scores.csv")
    with open(out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "score", "status", "seed"])
        w.writeheader(); w.writerows(rows)
    n_ok = sum(1 for r in rows if r["status"].startswith("ok"))
    n_act = sum(1 for r in rows if r["label"] == "1"); n_dec = len(rows) - n_act
    summary = {"target": target, "date": "2026-11-17", "n_rows": len(rows), "n_act": n_act, "n_dec": n_dec,
               "n_ok": n_ok, "n_fail": len(rows) - n_ok,
               "fail_status": {s: sum(1 for r in rows if r["status"] == s)
                               for s in set(r["status"] for r in rows) if not s.startswith("ok")},
               "box": box, "box_size": 22, "receptor": RECEPTORS[target], "vina": "1.2.5",
               "exhaustiveness": 8, "seed": SEED, "num_modes": 3,
               "workers": nproc, "cpu_per_worker": cpu,
               "elapsed_s": round(time.time() - t0, 1), "source": "B-line _b_trackA_full.py"}
    json.dump(summary, open(os.path.join(d, "summary.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1), flush=True)

def progress():
    rows = []
    for t in ["LpxC_Ec", "LpxC_Pa", "DHFR_Ec", "Gyrase_Ec", "GyrA_Ec", "GyrB_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]:
        p = os.path.join(OUTROOT, t, "summary.json")
        if os.path.exists(p):
            s = json.load(open(p, encoding="utf-8"))
            rows.append({"target": t, "n_rows": s.get("n_rows"), "n_ok": s.get("n_ok"),
                         "elapsed_s": s.get("elapsed_s"), "status": "done"})
        else:
            m = os.path.join(OUTROOT, t, "manifest.csv")
            n = len(list(csv.DictReader(open(m, encoding="utf-8-sig")))) if os.path.exists(m) else 0
            rows.append({"target": t, "n_rows": n, "n_ok": 0, "elapsed_s": "", "status": "not_started"})
    p = os.path.join(OUTROOT, "progress.csv")
    with open(p, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["target", "n_rows", "n_ok", "elapsed_s", "status"]); w.writeheader(); w.writerows(rows)
    print(json.dumps(rows, indent=1), flush=True)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("manifest"); p1.add_argument("--target", required=True); p1.add_argument("--limit", type=int)
    p2 = sub.add_parser("prep"); p2.add_argument("--target", required=True); p2.add_argument("--manifest", required=True); p2.add_argument("--nproc", type=int, default=10)
    p3 = sub.add_parser("dock"); p3.add_argument("--target", required=True); p3.add_argument("--manifest", required=True); p3.add_argument("--nproc", type=int, default=8); p3.add_argument("--cpu", type=int, default=2)
    p4 = sub.add_parser("prep-child"); p4.add_argument("--manifest", required=True); p4.add_argument("--chunk", type=int, required=True); p4.add_argument("--nchunks", type=int, required=True)
    p5 = sub.add_parser("dock-child"); p5.add_argument("--target", required=True); p5.add_argument("--manifest", required=True); p5.add_argument("--chunk", type=int, required=True); p5.add_argument("--nchunks", type=int, required=True); p5.add_argument("--cpu", type=int, default=2); p5.add_argument("--box", required=True)
    sub.add_parser("progress")
    a = ap.parse_args()
    if a.cmd == "manifest":
        build_manifest(a.target, a.limit)
    elif a.cmd == "prep":
        cmd_prep(a.target, a.manifest, a.nproc)
    elif a.cmd == "dock":
        cmd_dock(a.target, a.manifest, a.nproc, a.cpu)
    elif a.cmd == "prep-child":
        prep_child(a.manifest, a.chunk, a.nchunks)
    elif a.cmd == "dock-child":
        dock_child(a.target, a.manifest, a.chunk, a.nchunks, a.cpu, a.box)
    elif a.cmd == "progress":
        progress()
