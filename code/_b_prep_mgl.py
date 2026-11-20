#!/usr/bin/env python
"""MGLTools fallback for meeko-failed decoy ligands (Lead request).

For rows in <target>/prep_failures.csv with reason=meeko_fail and an existing .sdf:
  -> RDKit SDF->PDB, then MGLTools 1.5.7 prepare_ligand4.py -A hydrogens
     (run with cwd = ligand dir and RELATIVE filenames, per A-line note)
Updates prep_failures.csv (rescued rows removed) and prints final counts.
Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\_b_prep_mgl.py --target LpxC_Ec
"""
import argparse, csv, json, os, shutil, subprocess, sys

ROOT = r"D:\harness\projects\d1_gramneg_vs"
OUTROOT = os.path.join(ROOT, "results", "trackA", "full")
MGL_PY = r"C:\Program Files (x86)\MGLTools-1.5.7\python.exe"
MGL_LIG = r"C:\Program Files (x86)\MGLTools-1.5.7\Lib\site-packages\AutoDockTools\Utilities24\prepare_ligand4.py"

def main(target):
    d = os.path.join(OUTROOT, target)
    pf = os.path.join(d, "prep_failures.csv")
    if not os.path.exists(pf):
        sys.exit("no prep_failures.csv")
    rows = list(csv.DictReader(open(pf, encoding="utf-8-sig")))
    man = {r["ligand_id"]: r for r in csv.DictReader(open(os.path.join(d, "manifest.csv"), encoding="utf-8-sig"))}
    from rdkit import Chem
    rescued, still = [], []
    for r in rows:
        lid = r["ligand_id"]
        m = man.get(lid, {})
        pdbqt = m.get("pdbqt", "")
        sdf = pdbqt.replace(".pdbqt", ".sdf") if pdbqt else ""
        ok = False
        try:
            if sdf and os.path.exists(sdf):
                mol = Chem.MolFromMolFile(sdf, removeHs=False)
                if mol is not None:
                    pdb = sdf.replace(".sdf", ".pdb")
                    Chem.MolToPDBFile(mol, pdb)
                    lig_dir = os.path.dirname(pdb) or "."
                    rel_pdb, rel_out = os.path.basename(pdb), os.path.basename(pdbqt)
                    res = subprocess.run([MGL_PY, MGL_LIG, "-l", rel_pdb, "-o", rel_out, "-A", "hydrogens"],
                                         capture_output=True, text=True, encoding="utf-8",
                                         errors="replace", timeout=300, cwd=lig_dir)
                    if res.returncode == 0 and os.path.exists(pdbqt) and os.path.getsize(pdbqt) > 100:
                        ok = True
                    else:
                        r["reason"] = "mgl_fail:" + ((res.stderr or res.stdout or "")[-80:].replace("\n", " "))
            else:
                r["reason"] = "no_sdf"
        except Exception as e:
            r["reason"] = f"mgl_exception:{type(e).__name__}"
        if ok:
            rescued.append(lid)
        else:
            still.append(r)
    with open(pf, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "reason"]); w.writeheader(); w.writerows(still)
    fa = sum(1 for x in still if x["label"] == "1"); fd = len(still) - fa
    print(json.dumps({"target": target, "rescue_mgl": len(rescued), "final_fail": len(still),
                      "final_fail_actives": fa, "final_fail_decoys": fd,
                      "fail_reasons": {rr: sum(1 for x in still if x["reason"].startswith(rr)) for rr in set(x["reason"].split(":")[0] for x in still)},
                      "rescued_examples": rescued[:10]}, indent=1, ensure_ascii=False), flush=True)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--target", required=True)
    main(ap.parse_args().target)
