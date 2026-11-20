#!/usr/bin/env python
"""Prepare method-axis receptor PDBs for Amendment #12 (DiffDock / Boltz-2).

READ-ONLY over existing receptors; produces ONLY prep/validation artifacts
(no docking scores, no predictions).

Policy (Amendment #12 v2, 2026-11-19):
  * DiffDock receptor = SAME PDBQT used by Vina, converted to PDB, preserving
    existing HETATM records (Zn / Mg / NDP / NAD).  Sanity check on LpxC_Ec + DHFR_Ec:
    heavy-atom count identical, coordinates identical (<=1e-3 A), metals/cofactors present.
    Note: DiffDock parses receptors with ProDy `parsePDB(...).ca`, i.e. protein C-alpha
    atoms only -> HETATM cofactors/metals are NOT modelled even though present in input.
  * Boltz-2 input = SINGLE-CHAIN protein sequence of the chain nearest to the Vina grid
    centre (docking-relevant chain); multi-chain receptors record all chains and the
    selected chain; `msa: empty` (single sequence), no ion entity (documented limitation).

Outputs (only when executed):
  results/trackA/method_axis/receptors/<target>.pdb
  results/trackA/method_axis/receptor_prep.json

Run (dry-run allowed; no scoring):
  tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\make_method_axis_receptors.py \
      --prody-python tools\\envs\\diffdock2\\Scripts\\python.exe
"""
import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile

ROOT = r"D:\harness\projects\d1_gramneg_vs"
OUTDIR = os.path.join(ROOT, "results", "trackA", "method_axis")
RECEPTOR_DIR = os.path.join(OUTDIR, "receptors")
MANIFEST = os.path.join(ROOT, "data", "receptors", "receptor_manifest.csv")

TARGETS = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec",
           "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]
RECEPTOR_FILE = {
    "LpxC_Ec": "LpxC_Ec_holo.pdbqt", "LpxC_Pa": "LpxC_Pa_holo.pdbqt",
    "DHFR_Ec": "DHFR_Ec_holo_3DAU.pdbqt", "GyrB_Ec": "GyrB_Ec_holo.pdbqt",
    "Gyrase_Ec": "GyrB_Ec_holo.pdbqt", "GyrA_Ec": "GyrA_Ec_holo.pdbqt",
    "FabI_Ec": "FabI_Ec_holo.pdbqt", "TEM1_Ec": "TEM1_Ec_holo.pdbqt",
    "LepB_Ec": "LepB_Ec_holo.pdbqt",
}
# Gyrase_Ec has no own grid centre in the manifest (surrogate = GyrB_Ec subunit receptor)
GRID_ROW = {"Gyrase_Ec": "GyrB_Ec"}
SANITY_TARGETS = ["LpxC_Ec", "DHFR_Ec"]
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
# PDBQT atom type -> element (case-sensitive; PDBQT uses e.g. "NA"=acceptor N, "Na"=sodium)
ATOM_TYPE_ELEMENT = {
    "H": "H", "HD": "H", "HS": "H", "HO": "H",
    "C": "C", "A": "C", "CG0": "C", "CG1": "C", "CG2": "C", "CG3": "C",
    "N": "N", "NA": "N", "NS": "N", "NP": "N",
    "O": "O", "OA": "O", "OS": "O", "OXT": "O",
    "S": "S", "SA": "S", "P": "P",
    "F": "F", "Cl": "Cl", "CL": "Cl", "Br": "Br", "BR": "Br", "I": "I",
    "Zn": "Zn", "ZN": "Zn", "Mg": "Mg", "MG": "Mg", "Mn": "Mn", "MN": "Mn",
    "Fe": "Fe", "FE": "Fe", "Ca": "Ca", "CA": "Ca", "K": "K", "Na": "Na",
}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_grid_centres():
    rows = list(csv.DictReader(open(MANIFEST, encoding="utf-8-sig")))
    by = {}
    for r in rows:
        if str(r.get("role")) != "holo":
            continue
        key = r["target_key"]
        try:
            c = (float(r["grid_center_x"]), float(r["grid_center_y"]), float(r["grid_center_z"]))
        except (KeyError, TypeError, ValueError):
            continue
        by[key] = {"center": c, "box_size": float(r.get("box_size") or 22.0),
                   "pdbqt": r.get("pdbqt", ""), "metal_handling": r.get("metal_handling", "")}
    out = {}
    for t in TARGETS:
        src = GRID_ROW.get(t, t)
        if src not in by:
            raise RuntimeError(f"no grid centre for {t} (looked up {src})")
        out[t] = dict(by[src])
        out[t]["grid_source"] = src
    return out


def parse_pdbqt(path):
    """Parse ATOM/HETATM records; return atoms with element + type + residue metadata."""
    atoms = []
    for ln in open(path, encoding="utf-8", errors="replace"):
        rec = ln[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        atype = ln[77:79].strip() if len(ln) > 77 else ""
        # robust fallback: atom type is last whitespace token on PDBQT lines
        toks = ln.split()
        if not atype and toks:
            atype = toks[-1]
        atoms.append({
            "record": rec,
            "serial": ln[6:11].strip(),
            "name": ln[12:16],
            "altloc": ln[16],
            "resname": ln[17:20].strip(),
            "chain": ln[21] if len(ln) > 21 else " ",
            "resseq": ln[22:26].strip(),
            "icode": ln[26] if len(ln) > 26 else " ",
            "xyz": (float(ln[30:38]), float(ln[38:46]), float(ln[46:54])),
            "atom_type": atype,
            "element": ATOM_TYPE_ELEMENT.get(atype, atype[:1].upper() if atype else "C"),
            "is_h": ATOM_TYPE_ELEMENT.get(atype, "") == "H",
        })
    return atoms


def pdb_line(a, serial):
    el = a["element"]
    # atom-name alignment: 1-char element starts at col 14, 2-char at col 13
    name = a["name"]
    if not name or name.strip() == "":
        name = " %-3s" % (el.upper() if len(el) == 1 else el)
    if len(el) == 1:
        name = (" " + name.strip())[:4].ljust(4)
    else:
        name = name.strip()[:4].ljust(4)
    resn = (a["resname"] or "UNK")[:3].rjust(3)
    ch = (a["chain"] or " ")[:1]
    resseq = str(a["resseq"] or "1")[:4].rjust(4)
    icode = (a["icode"] or " ")[:1]
    x, y, z = a["xyz"]
    line = (f"{a['record']:<6}{serial:>5} {name}{' ':1}{resn} {ch}"
            f"{resseq}{icode}   {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
            f"          {el:>2}")
    return line.rstrip() + "\n"


def write_pdb(atoms, path):
    lines = []
    serial = 0
    for a in atoms:
        serial += 1
        if serial > 99999:
            serial = 99999
        lines.append(pdb_line(a, serial))
    lines.append("END\n")
    with open(path, "w", encoding="ascii", errors="replace") as f:
        f.writelines(lines)


def chain_info_from_atoms(atoms, center):
    """Per-chain: protein residue list (ATOM only, standard AA), sequence, min CA distance to center."""
    import math
    chains = {}
    for a in atoms:
        if a["record"] != "ATOM":
            continue
        if a["resname"] not in AA3:
            continue
        ch = a["chain"]
        key = (ch, a["resseq"], a["icode"])
        d = chains.setdefault(ch, {"res": {}, "ca": []})
        d["res"][key] = AA3[a["resname"]]
        if a["name"].strip() == "CA":
            d["ca"].append(a["xyz"])
    out = {}
    for ch, d in chains.items():
        order = sorted(d["res"], key=lambda k: (str(k[0]), int(k[1]) if str(k[1]).isdigit() else 0, k[2]))
        seq = "".join(d["res"][k] for k in order)
        dists = [math.dist(p, center) for p in d["ca"]] if d["ca"] else []
        out[ch] = {"n_residues": len(order), "sequence": seq,
                   "min_ca_dist_to_box_A": round(min(dists), 2) if dists else None}
    return out


def prody_check(prody_python, pdb_path):
    """Use DiffDock's parser stack (ProDy) to confirm PDB parses and determine HETATM modelling."""
    code = """
import json, sys
try:
    import prody
    p = prody.parsePDB(sys.argv[1], altloc='all')
    ca = p.ca
    if ca is None:
        print(json.dumps({"prody_ok": False, "reason": "no CA atoms"})); sys.exit(0)
    chains = list(dict.fromkeys(ca.getChids()))
    resnames = set()
    for a in ca:
        resnames.add(a.getResname())
    print(json.dumps({"prody_ok": True, "ca_atoms": int(ca.numAtoms()),
                      "chains": chains, "n_chains": len(chains),
                      "ca_resnames": sorted(resnames)[:8]}))
except Exception as e:
    print(json.dumps({"prody_ok": False, "reason": repr(e)[:200]}))
"""
    fd, tmp = tempfile.mkstemp(suffix=".py", dir=os.path.dirname(pdb_path) or None)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(code)
        r = subprocess.run([prody_python, tmp, pdb_path], capture_output=True, text=True, timeout=120)
        lines = [ln for ln in r.stdout.strip().splitlines() if ln.strip().startswith("{")]
        return json.loads(lines[-1]) if lines else {"prody_ok": False, "reason": r.stderr[-200:]}
    except Exception as e:
        return {"prody_ok": False, "reason": repr(e)[:200]}
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass


def build(prody_python=None):
    centres = read_grid_centres()
    os.makedirs(RECEPTOR_DIR, exist_ok=True)
    report = {
        "amendment": "#12 v2 method-class complementarity axis",
        "policy": "DiffDock: convert the SAME PDBQT used by Vina to PDB, keep HETATM (Zn/Mg/NDP/NAD). "
                  "Boltz-2: single-chain sequence of the chain nearest to Vina grid centre; msa=empty; no ion entity.",
        "conversion": "PDBQT -> PDB, coordinates verbatim (8.3f), occupancy/tempFactor 1.00/0.00, element from PDBQT atom type",
        "receptor_file": {t: os.path.relpath(os.path.join(ROOT, "data", "receptors", RECEPTOR_FILE[t]), ROOT) for t in TARGETS},
        "targets": {},
        "sanity_checks": {},
        "note": "Prep-only artifacts; no docking scores / predictions generated.",
    }
    for t in TARGETS:
        pdbqt = os.path.join(ROOT, "data", "receptors", RECEPTOR_FILE[t])
        atoms = parse_pdbqt(pdbqt)
        out_pdb = os.path.join(RECEPTOR_DIR, f"{t}.pdb")
        write_pdb(atoms, out_pdb)
        heavy = [a for a in atoms if not a["is_h"]]
        het = [a for a in atoms if a["record"] == "HETATM"]
        het_res = sorted({a["resname"] for a in het})
        metal_types = sorted({a["atom_type"] for a in het})
        chains_pdb = sorted({a["chain"] for a in atoms})
        info = chain_info_from_atoms(atoms, centres[t]["center"])
        chosen = None
        if info:
            valid = [(ch, d["min_ca_dist_to_box_A"]) for ch, d in info.items() if d["min_ca_dist_to_box_A"] is not None]
            if valid:
                chosen = min(valid, key=lambda x: x[1])[0]
        # coordinate fidelity: same order, first min(len) atoms
        max_dev = 0.0
        for i in range(min(len(atoms), len(atoms))):
            pass
        max_dev = 0.0  # conversion writes coordinates verbatim; validate by round-trip read below
        # validate PDB heavy-atom count by re-reading written file
        n_pdb_atoms = sum(1 for ln in open(out_pdb, encoding="ascii", errors="replace")
                          if ln.startswith(("ATOM", "HETATM")))
        n_pdb_heavy = sum(1 for ln in open(out_pdb, encoding="ascii", errors="replace")
                          if ln.startswith(("ATOM", "HETATM")) and ln[76:78].strip() != "H")
        report["targets"][t] = {
            "target": t,
            "pdbqt": os.path.relpath(pdbqt, ROOT),
            "pdbqt_sha256": sha256_file(pdbqt),
            "pdb": os.path.relpath(out_pdb, ROOT),
            "pdb_sha256": sha256_file(out_pdb),
            "grid_center": list(centres[t]["center"]),
            "grid_source": centres[t]["grid_source"],
            "box_size": centres[t]["box_size"],
            "n_atoms_pdbqt": len(atoms),
            "n_heavy_pdbqt": len(heavy),
            "n_hydrogen_pdbqt": len(atoms) - len(heavy),
            "n_atoms_pdb": n_pdb_atoms,
            "n_heavy_pdb": n_pdb_heavy,
            "chains": chains_pdb,
            "chain_info": info,
            "boltz2_chain": chosen,
            "boltz2_sequence": info[chosen]["sequence"] if chosen else "",
            "boltz2_sequence_len": len(info[chosen]["sequence"]) if chosen else 0,
            "hetatm_resnames": het_res,
            "hetatm_atom_types": metal_types,
            "hetatm_n_atoms": len(het),
            "cofactors_retained": [r for r in het_res if r in ("NDP", "NAD")],
            "metals_retained": [r for r in het_res if r in ("ZN", "MG", "MN", "FE", "CA", "K", "NA")],
            "chain_selection_rule": "chain with CA nearest to Vina grid centre",
            "multichain": len(chains_pdb) > 1,
            "diffdock_receptor_chain_note": "DiffDock reads the full PDB with ProDy (protein C-alpha only); "
                                            "all chains in this file are used, incl. multimeric copies",
        }
    # post-hoc coordinate round-trip fidelity check (write -> parse)
    for t in TARGETS:
        pdbqt = os.path.join(ROOT, "data", "receptors", RECEPTOR_FILE[t])
        atoms = parse_pdbqt(pdbqt)
        pdb = os.path.join(RECEPTOR_DIR, f"{t}.pdb")
        coords = []
        for ln in open(pdb, encoding="ascii", errors="replace"):
            if ln.startswith(("ATOM", "HETATM")):
                coords.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
        dev = max((abs(a["xyz"][0] - c[0]) + abs(a["xyz"][1] - c[1]) + abs(a["xyz"][2] - c[2]))
                  for a, c in zip(atoms, coords)) if atoms and len(coords) == len(atoms) else None
        report["targets"][t]["max_coord_deviation_A"] = round(dev, 6) if dev is not None else None
        report["targets"][t]["atom_count_match"] = (len(coords) == len(atoms))
        # HETATM rows preserved exactly by resname
        het_pdb = set()
        for ln in open(pdb, encoding="ascii", errors="replace"):
            if ln.startswith("HETATM"):
                het_pdb.add(ln[17:20].strip())
        report["targets"][t]["hetatm_resnames_in_pdb"] = sorted(het_pdb)
        report["targets"][t]["hetatm_preserved"] = (set(report["targets"][t]["hetatm_resnames"]) == het_pdb)
    # prody / DiffDock parser check
    if prody_python:
        for t in TARGETS:
            r = prody_check(prody_python, os.path.join(RECEPTOR_DIR, f"{t}.pdb"))
            report["targets"][t]["prody_parser"] = r
    # parser flags: ProDy CA atoms should equal number of protein residues; HETATM never modelled
    for t in TARGETS:
        d = report["targets"][t]
        pr = d.get("prody_parser") or {}
        n_prot_res = sum(c["n_residues"] for c in d["chain_info"].values())
        d["diffdock_parser"] = {
            "prody_ok": bool(pr.get("prody_ok")),
            "ca_atoms": pr.get("ca_atoms"),
            "ca_equals_protein_residues": (pr.get("ca_atoms") == n_prot_res) if pr.get("prody_ok") else None,
            "parser_models_hetatm": False if (pr.get("prody_ok") and d["hetatm_n_atoms"] > 0) else None,
            "note": "DiffDock receptor parser uses ProDy parsePDB(...).ca (protein C-alpha only); "
                    "HETATM records (metals/cofactors) are kept in the input PDB but not modelled.",
        }
    # strict sanity checks for 2 pre-registered targets
    for t in SANITY_TARGETS:
        d = report["targets"][t]
        ok = (d["atom_count_match"] and d["max_coord_deviation_A"] is not None
              and d["max_coord_deviation_A"] <= 1e-3
              and d["n_heavy_pdb"] == d["n_heavy_pdbqt"]
              and d["hetatm_preserved"]
              and bool(d["metals_retained"] or d["cofactors_retained"]))
        report["sanity_checks"][t] = {
            "ok": bool(ok),
            "heavy_atom_count_pdbqt": d["n_heavy_pdbqt"],
            "heavy_atom_count_pdb": d["n_heavy_pdb"],
            "heavy_atom_count_match": d["n_heavy_pdb"] == d["n_heavy_pdbqt"],
            "atom_count_match": d["atom_count_match"],
            "max_coord_deviation_A": d["max_coord_deviation_A"],
            "metals_retained": d["metals_retained"],
            "cofactors_retained": d["cofactors_retained"],
            "hetatm_preserved": d["hetatm_preserved"],
            "diffdock_parser_ok": d["diffdock_parser"]["prody_ok"],
            "diffdock_ca_equals_protein_residues": d["diffdock_parser"]["ca_equals_protein_residues"],
            "parser_models_hetatm": d["diffdock_parser"]["parser_models_hetatm"],
            "note": "DiffDock/ProDy receptor parser uses protein C-alpha atoms only; "
                    "retained HETATM (Zn/Mg/NDP/NAD) are present in the input PDB but not modelled by DiffDock.",
        }
    with open(os.path.join(OUTDIR, "receptor_prep.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1, ensure_ascii=False)
    print("wrote", os.path.join(OUTDIR, "receptor_prep.json"))
    for t in TARGETS:
        d = report["targets"][t]
        print(f"{t:10s} pdb={os.path.basename(d['pdb']):16s} atoms={d['n_atoms_pdbqt']:6d} "
              f"heavy={d['n_heavy_pdbqt']:6d} het={d['hetatm_n_atoms']:4d} {d['hetatm_resnames']} "
              f"chains={d['chains']} boltz_chain={d['boltz2_chain']} seq_len={d['boltz2_sequence_len']} "
              f"dev={d['max_coord_deviation_A']}")
    print("sanity:", json.dumps(report["sanity_checks"], ensure_ascii=False)[:500])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prody-python", default=None,
                    help="python with prody installed (DiffDock env) to test parser; e.g. tools\\envs\\diffdock2\\Scripts\\python.exe")
    args = ap.parse_args()
    build(args.prody_python)


if __name__ == "__main__":
    main()
