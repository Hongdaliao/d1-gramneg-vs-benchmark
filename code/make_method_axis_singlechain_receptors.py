#!/usr/bin/env python
"""Extract single-chain sensitivity receptors for Amendment #12 multichain targets.

PREP-ONLY (no docking / no predictions). For GyrA_Ec (2 chains), FabI_Ec (8 chains),
LepB_Ec (4 chains) DiffDock main analysis uses the full Vina receptor (all chains).
The pre-registered single-chain sensitivity arm uses only the chain whose CA is
nearest to the Vina grid centre (chain recorded in receptor_prep.json).

Inputs  (frozen): results/trackA/method_axis/receptor_prep.json + receptors/<target>.pdb
Outputs          : results/trackA/method_axis/receptors_singlechain/<target>.pdb
                   results/trackA/method_axis/receptor_singlechain_prep.json  (sha256 + checks)
"""
import json
import os
import hashlib
import subprocess
import tempfile

ROOT = r"D:\harness\projects\d1_gramneg_vs"
BASE = os.path.join(ROOT, "results", "trackA", "method_axis")
PREP = os.path.join(BASE, "receptor_prep.json")
OUTDIR = os.path.join(BASE, "receptors_singlechain")
MULTICHAIN = ["GyrA_Ec", "FabI_Ec", "LepB_Ec"]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def prody_check(prody_python, pdb_path):
    code = """
import json, sys
import prody
p = prody.parsePDB(sys.argv[1], altloc='all')
ca = p.ca
print(json.dumps({"prody_ok": ca is not None,
                  "ca_atoms": int(ca.numAtoms()) if ca is not None else 0,
                  "chains": list(dict.fromkeys(ca.getChids())) if ca is not None else []}))
"""
    fd, tmp = tempfile.mkstemp(suffix=".py", dir=OUTDIR)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(code)
        r = subprocess.run([prody_python, tmp, pdb_path], capture_output=True, text=True, timeout=120)
        lines = [ln for ln in r.stdout.strip().splitlines() if ln.strip().startswith("{")]
        return json.loads(lines[-1]) if lines else {"prody_ok": False, "stderr": r.stderr[-150:]}
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass


def main(prody_python=None):
    prep = json.load(open(PREP, encoding="utf-8"))
    os.makedirs(OUTDIR, exist_ok=True)
    report = {
        "amendment": "#12 v3 multichain single-chain sensitivity arm (exploratory)",
        "policy": "extract only the chain with CA nearest to Vina grid centre from the frozen full receptor PDB; "
                  "protein ATOM records only (DiffDock parser models protein CA only); no docking scoring",
        "input_receptor_prep_sha256": sha256_file(PREP),
        "targets": {},
    }
    for t in MULTICHAIN:
        d = prep["targets"][t]
        full_pdb = os.path.join(ROOT, d["pdb"])
        chain = d["boltz2_chain"]
        assert chain, f"no chosen chain for {t}"
        assert len(d["chains"]) > 1, f"{t} is not multichain"
        lines, n_atom = [], 0
        for ln in open(full_pdb, encoding="ascii", errors="replace"):
            if ln.startswith("ATOM") and ln[21] == chain:
                lines.append(ln)
                n_atom += 1
        if not lines:
            raise RuntimeError(f"no ATOM records for chain {chain} in {full_pdb}")
        lines.append("END\n")
        out = os.path.join(OUTDIR, f"{t}_singlechain.pdb")
        with open(out, "w", encoding="ascii", errors="replace") as f:
            f.writelines(lines)
        heavy = sum(1 for ln in lines if ln.startswith("ATOM") and ln[76:78].strip() != "H")
        info = d["chain_info"][chain]
        rec = {
            "target": t,
            "full_receptor_pdb": d["pdb"],
            "full_receptor_chains": d["chains"],
            "selected_chain": chain,
            "selection_rule": "chain with CA nearest to Vina grid centre (frozen in receptor_prep.json)",
            "min_ca_dist_to_box_A": info["min_ca_dist_to_box_A"],
            "singlechain_pdb": os.path.relpath(out, ROOT),
            "singlechain_sha256": sha256_file(out),
            "n_atoms": n_atom,
            "n_heavy": heavy,
            "n_residues_expected": info["n_residues"],
            "sequence_len_expected": len(info["sequence"]),
        }
        if prody_python:
            rec["prody_parser"] = prody_check(prody_python, out)
            rec["ca_equals_residues"] = (rec["prody_parser"].get("ca_atoms") == info["n_residues"])
        report["targets"][t] = rec
        print(f"{t:10s} chain={chain} atoms={n_atom} heavy={heavy} "
              f"residues={info['n_residues']} prody={rec.get('prody_parser')}")
    with open(os.path.join(BASE, "receptor_singlechain_prep.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1, ensure_ascii=False)
    print("wrote", os.path.join(BASE, "receptor_singlechain_prep.json"))
    print("singlechain prep sha256:", sha256_file(os.path.join(BASE, "receptor_singlechain_prep.json")))


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--prody-python", default=None)
    args = ap.parse_args()
    main(args.prody_python)
