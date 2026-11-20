#!/usr/bin/env python
"""LIT-PCBA external anchor driver (task-31, frozen protocol v1.1 + Amendment A1).

Protocol (data/litpcba/protocol_freeze.txt):
  * test actives = official <set>_actives_val.csv
  * negatives    = official inactive_V stratified subsample <set>_inactives_s{42,43,44}.csv
  * receptor     = data/litpcba/receptors/<set>_<pdb>.pdbqt (A-line prepared)
  * box          = panel.csv grid_center, 22^3
  * Vina         = 1.2.5, exhaustiveness 8, top-9 modes, seeds 42/43/44 (Vina seed replicates)
  * metrics      = EF5% primary (EF1% only if n_act_val>=150); BEDROC(a=20); AUC; AC1-AC7

Subcommands: manifest|prep|dock|metrics
Run (chain validation, ADRB2 single subsample seed + single Vina seed):
  tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\_b_litpcba.py prep --set ADRB2 --subseed 42
  ... dock --set ADRB2 --subseed 42 --vseed 42 --cpu 16
"""
import argparse, csv, json, math, os, re, subprocess, sys, time

ROOT = r"D:\harness\projects\d1_gramneg_vs"
LIT = os.path.join(ROOT, "data", "litpcba")
OUT = os.path.join(ROOT, "results", "litpcba")
VINA = os.path.join(ROOT, "tools", "vina", "vina.exe")
MEEPREP = os.path.join(ROOT, "tools", "envs", "boltz", "Scripts", "mk_prepare_ligand.exe")
os.makedirs(OUT, exist_ok=True)

def now():
    return time.strftime("%Y-%m-%dT%H:%M:%S")

def set_dir(s):
    d = os.path.join(OUT, s)
    for sub in ("ligands", "poses", "logs"):
        os.makedirs(os.path.join(d, sub), exist_ok=True)
    return d

def load_panel():
    return {r["target_set"]: r for r in csv.DictReader(open(os.path.join(LIT, "panel.csv"), encoding="utf-8-sig"))}

def load_receptors():
    return {r["target_set"]: r for r in csv.DictReader(open(os.path.join(LIT, "receptors.csv"), encoding="utf-8-sig"))}

def build_manifest(s, subseed, ratio=20, cap=1000):
    d = set_dir(s)
    acts = list(csv.DictReader(open(os.path.join(LIT, "splits", f"{s}_actives_val.csv"), encoding="utf-8-sig")))
    inacts = list(csv.DictReader(open(os.path.join(LIT, "subsamples", f"{s}_inactives_s{subseed}.csv"), encoding="utf-8-sig")))
    n_want = max(0, min(ratio * len(acts), cap, len(inacts)))
    inacts = inacts[:n_want]  # subsample file is already stratified-random; A2: 20:1 cap 1000
    rows = []
    for r in acts:
        sid = r.get("pubchem_sid") or r.get("rank") or f"a{len(rows)}"
        rows.append({"ligand_id": "act_" + str(sid), "label": 1, "smiles": r["smiles"].strip(), "source": "official_val"})
    for r in inacts:
        sid = r.get("pubchem_sid") or r.get("rank") or f"i{len(rows)}"
        rows.append({"ligand_id": "inact_" + str(sid), "label": 0, "smiles": r["inactive_smiles"].strip(), "source": "official_inactive_sub"})
    for r in rows:
        r["pdbqt"] = os.path.join(d, "ligands", r["ligand_id"] + ".pdbqt")
    man = os.path.join(d, f"manifest_s{subseed}.csv")
    with open(man, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "smiles", "pdbqt", "source"]); w.writeheader(); w.writerows(rows)
    print(json.dumps({"set": s, "subseed": subseed, "n_act": len(acts), "n_inact": len(inacts), "manifest": man}), flush=True)
    return man

def prep_child(man, chunk, nchunks):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    rows = list(csv.DictReader(open(man, encoding="utf-8-sig")))
    mine = [r for i, r in enumerate(rows) if i % nchunks == chunk]
    ok = skip = fail = 0
    for r in mine:
        if os.path.exists(r["pdbqt"]) and os.path.getsize(r["pdbqt"]) > 100:
            skip += 1; continue
        try:
            m = Chem.MolFromSmiles(r["smiles"])
            if m is None:
                fail += 1; continue
            m = Chem.AddHs(m)
            p = AllChem.ETKDGv3(); p.randomSeed = 42
            if AllChem.EmbedMolecule(m, p) != 0:
                AllChem.EmbedMolecule(m, useRandomCoords=True, randomSeed=42)
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
    print(json.dumps({"chunk": chunk, "ok": ok, "skip": skip, "fail": fail}), flush=True)

def launch(args_list, logs):
    procs = []
    for args, log in zip(args_list, logs):
        os.makedirs(os.path.dirname(log), exist_ok=True)
        lf = open(log, "w", encoding="utf-8", errors="replace")
        procs.append((subprocess.Popen([sys.executable, os.path.abspath(__file__)] + args,
                                       stdout=lf, stderr=subprocess.STDOUT, cwd=ROOT), lf))
    for p, lf in procs:
        p.wait(); lf.close()

def cmd_prep(s, subseed, nproc, ratio=20, cap=1000):
    man = build_manifest(s, subseed, ratio, cap)
    d = set_dir(s)
    args = [["prep-child", "--manifest", man, "--chunk", str(c), "--nchunks", str(nproc)] for c in range(nproc)]
    logs = [os.path.join(d, "logs", f"prep_{c}.log") for c in range(nproc)]
    t0 = time.time(); launch(args, logs)
    rows = list(csv.DictReader(open(man, encoding="utf-8-sig")))
    n_ok = sum(1 for r in rows if os.path.exists(r["pdbqt"]) and os.path.getsize(r["pdbqt"]) > 100)
    print(json.dumps({"prep_ok": n_ok, "n": len(rows), "elapsed_s": round(time.time()-t0,1)}), flush=True)

def dock_child(s, man, chunk, nchunks, cpu, vseed, exh=4, subseed=42):
    d = set_dir(s)
    panel = load_panel()[s]; rec = load_receptors()[s]
    rec_path = os.path.join(ROOT, rec["pdbqt_path"])
    box = [float(panel["grid_center"].strip("[]").split(",")[i]) for i in range(3)]
    rows = list(csv.DictReader(open(man, encoding="utf-8-sig")))
    mine = [r for i, r in enumerate(rows) if i % nchunks == chunk]
    outp = os.path.join(d, f"scores_vseed{vseed}_part{chunk}.csv")
    with open(outp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f); w.writerow(["ligand_id", "label", "score", "status", "vina_seed"])
        for r in mine:
            lid = r["ligand_id"]
            if not (os.path.exists(r["pdbqt"]) and os.path.getsize(r["pdbqt"]) > 100):
                w.writerow([lid, r["label"], "", "prep_missing", vseed]); f.flush(); continue
            pose = os.path.join(d, "poses", f"{lid}_sub{subseed}_v{vseed}_out.pdbqt")
            logf = os.path.join(d, "logs", f"vina_{lid}_v{vseed}.log")
            cmd = [VINA, "--receptor", rec_path, "--ligand", r["pdbqt"],
                   "--center_x", str(box[0]), "--center_y", str(box[1]), "--center_z", str(box[2]),
                   "--size_x", "22", "--size_y", "22", "--size_z", "22",
                   "--exhaustiveness", str(exh), "--num_modes", "9", "--seed", str(vseed), "--cpu", str(cpu), "--out", pose]
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
            w.writerow([lid, r["label"], "" if score is None else score,
                        "ok" if score is not None else "fail", vseed]); f.flush()
    print(json.dumps({"chunk": chunk, "done": len(mine)}), flush=True)

def cmd_dock(s, subseed, vseed, nproc, cpu, exh=4):
    d = set_dir(s)
    man = os.path.join(d, f"manifest_s{subseed}.csv")
    args = [["dock-child", "--set", s, "--manifest", man, "--chunk", str(c), "--nchunks", str(nproc),
             "--cpu", str(cpu), "--vseed", str(vseed), "--exhaustiveness", str(exh), "--subseed", str(subseed)] for c in range(nproc)]
    logs = [os.path.join(d, "logs", f"dock_v{vseed}_{c}.log") for c in range(nproc)]
    t0 = time.time(); launch(args, logs)
    rows = []
    for c in range(nproc):
        p = os.path.join(d, f"scores_s{subseed}_vseed{vseed}_part{c}.csv")
        if os.path.exists(p):
            rows += list(csv.DictReader(open(p, encoding="utf-8-sig")))
    sp = os.path.join(d, f"scores_s{subseed}_v{vseed}.csv")
    with open(sp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "score", "status", "vina_seed"]); w.writeheader(); w.writerows(rows)
    n_ok = sum(1 for r in rows if r["status"] == "ok")
    print(json.dumps({"set": s, "subseed": subseed, "vseed": vseed, "n": len(rows), "n_ok": n_ok,
                      "elapsed_s": round(time.time()-t0,1)}, indent=1), flush=True)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("prep"); p1.add_argument("--set", required=True); p1.add_argument("--subseed", type=int, default=42); p1.add_argument("--nproc", type=int, default=12); p1.add_argument("--ratio", type=int, default=20); p1.add_argument("--cap", type=int, default=1000)
    p2 = sub.add_parser("dock"); p2.add_argument("--set", required=True); p2.add_argument("--subseed", type=int, default=42)
    p2.add_argument("--vseed", type=int, default=42); p2.add_argument("--nproc", type=int, default=16); p2.add_argument("--cpu", type=int, default=1); p2.add_argument("--exhaustiveness", type=int, default=4)
    p3 = sub.add_parser("prep-child"); p3.add_argument("--manifest", required=True); p3.add_argument("--chunk", type=int, required=True); p3.add_argument("--nchunks", type=int, required=True)
    p4 = sub.add_parser("dock-child"); p4.add_argument("--set", required=True); p4.add_argument("--manifest", required=True)
    p4.add_argument("--chunk", type=int, required=True); p4.add_argument("--nchunks", type=int, required=True)
    p4.add_argument("--cpu", type=int, default=1); p4.add_argument("--vseed", type=int, default=42); p4.add_argument("--exhaustiveness", type=int, default=4); p4.add_argument("--subseed", type=int, default=42)
    a = ap.parse_args()
    if a.cmd == "prep": cmd_prep(a.set, a.subseed, a.nproc, a.ratio, a.cap)
    elif a.cmd == "dock": cmd_dock(a.set, a.subseed, a.vseed, a.nproc, a.cpu, a.exhaustiveness)
    elif a.cmd == "prep-child": prep_child(a.manifest, a.chunk, a.nchunks)
    elif a.cmd == "dock-child": dock_child(a.set, a.manifest, a.chunk, a.nchunks, a.cpu, a.vseed, a.exhaustiveness, a.subseed)
