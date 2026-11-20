#!/usr/bin/env python
"""Retry ligand 3D prep for prep_missing rows (Lead request 2026-11-17).

For each manifest row without a valid PDBQT:
  fallback chain: ETKDGv3(useRandomCoords, maxAttempts=200) -> ETKDGv3(ignoreSmoothingFailures)
                  -> EmbedMolecule(useRandomCoords, ignoreSmoothingFailures)
  then meeko mk_prepare_ligand; record final reason on failure.
Outputs: results/trackA/full/<target>/prep_failures.csv (ligand_id,label,reason,rescued)
         + prints actives/decoys counts (before/after rescue).
Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\_b_prep_retry.py --target LpxC_Ec --nproc 4
"""
import argparse, csv, json, os, subprocess, sys
ROOT = r"D:\harness\projects\d1_gramneg_vs"
OUTROOT = os.path.join(ROOT, "results", "trackA", "full")
MEEPREP = os.path.join(ROOT, "tools", "envs", "boltz", "Scripts", "mk_prepare_ligand.exe")
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

def embed_and_write(smiles, sdf):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, "smiles_invalid"
    m = Chem.AddHs(m)
    ok = False
    try:
        p = AllChem.ETKDGv3(); p.randomSeed = SEED; p.useRandomCoords = True
        p.maxAttempts = 200; p.ignoreSmoothingFailures = True
        ok = AllChem.EmbedMolecule(m, p) == 0
    except Exception:
        ok = False
    if not ok:
        try:
            p = AllChem.ETKDGv3(); p.randomSeed = SEED; p.ignoreSmoothingFailures = True
            ok = AllChem.EmbedMolecule(m, p) == 0
        except Exception:
            ok = False
    if not ok:
        try:
            ok = AllChem.EmbedMolecule(m, useRandomCoords=True, ignoreSmoothingFailures=True, randomSeed=SEED) == 0
        except Exception:
            ok = False
    if not ok:
        return None, "embed_fail"
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=200)
    except Exception:
        pass
    Chem.MolToMolFile(m, sdf)
    return sdf, "ok"

def child(target, chunk, nchunks):
    d = os.path.join(OUTROOT, target)
    man = os.path.join(d, "manifest.csv")
    rows = list(csv.DictReader(open(man, encoding="utf-8-sig")))
    todo = [r for r in rows if not (os.path.exists(r["pdbqt"]) and os.path.getsize(r["pdbqt"]) > 100)]
    mine = [r for i, r in enumerate(todo) if i % nchunks == chunk]
    res = []
    for r in mine:
        sdf_fail_reason = ""
        try:
            sdf = r["pdbqt"].replace(".pdbqt", ".sdf")
            got, reason = embed_and_write(r["smiles"], sdf)
            if got is None:
                pdbqt_fail = reason
            else:
                pr = subprocess.run([MEEPREP, "-i", sdf, "-o", r["pdbqt"]], capture_output=True,
                                    text=True, encoding="utf-8", errors="replace", timeout=300)
                if pr.returncode == 0 and os.path.exists(r["pdbqt"]) and os.path.getsize(r["pdbqt"]) > 100:
                    pdbqt_fail = ""
                else:
                    pdbqt_fail = "meeko_fail"
        except Exception as e:
            pdbqt_fail = f"exception_{type(e).__name__}"
        res.append({"ligand_id": r["ligand_id"], "label": r["label"],
                    "reason": pdbqt_fail, "rescued": "True" if not pdbqt_fail else "False"})
    with open(os.path.join(d, f"prep_retry_part{chunk}.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "reason", "rescued"]); w.writeheader(); w.writerows(res)
    print(json.dumps({"chunk": chunk, "n": len(res), "rescued": sum(1 for x in res if x["rescued"] == "True")}), flush=True)

def main_cmd(target, nproc):
    d = os.path.join(OUTROOT, target)
    man = os.path.join(d, "manifest.csv")
    rows = list(csv.DictReader(open(man, encoding="utf-8-sig")))
    missing = [r for r in rows if not (os.path.exists(r["pdbqt"]) and os.path.getsize(r["pdbqt"]) > 100)]
    n_act = sum(1 for r in missing if r["label"] == "1")
    print(json.dumps({"target": target, "n_rows": len(rows), "n_missing_before": len(missing),
                      "missing_actives": n_act, "missing_decoys": len(missing) - n_act}), flush=True)
    args = [["retry-child", "--target", target, "--chunk", str(c), "--nchunks", str(nproc)] for c in range(nproc)]
    logs = [os.path.join(d, "logs", f"retry_{c}.log") for c in range(nproc)]
    launch(args, logs)
    parts = []
    for c in range(nproc):
        p = os.path.join(d, f"prep_retry_part{c}.csv")
        if os.path.exists(p):
            parts += list(csv.DictReader(open(p, encoding="utf-8-sig")))
    still = [x for x in parts if x["rescued"] != "True"]
    # final failures with reason
    failed = [{"ligand_id": x["ligand_id"], "label": x["label"], "reason": x["reason"] or "unknown"} for x in still]
    with open(os.path.join(d, "prep_failures.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "reason"]); w.writeheader(); w.writerows(failed)
    fa = sum(1 for x in failed if x["label"] == "1"); fd = len(failed) - fa
    print(json.dumps({"rescued": len(parts) - len(failed), "final_fail": len(failed),
                      "final_fail_actives": fa, "final_fail_decoys": fd,
                      "fail_reasons": {r: sum(1 for x in failed if x["reason"] == r) for r in set(x["reason"] for x in failed)},
                      "prep_failures_csv": os.path.join(d, "prep_failures.csv")}, indent=1), flush=True)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("retry"); p1.add_argument("--target", required=True); p1.add_argument("--nproc", type=int, default=4)
    p2 = sub.add_parser("retry-child"); p2.add_argument("--target", required=True); p2.add_argument("--chunk", type=int, required=True); p2.add_argument("--nchunks", type=int, required=True)
    a = ap.parse_args()
    if a.cmd == "retry":
        main_cmd(a.target, a.nproc)
    else:
        child(a.target, a.chunk, a.nchunks)
