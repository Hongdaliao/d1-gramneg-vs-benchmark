#!/usr/bin/env python
"""Round-2 feasibility: build DrugCLIP pocket inputs for Track A (do not touch frozen axis).

For each Track A holo target:
  * locate the original co-crystal ligand copy (resname == manifest ligand_comp_id) in the raw PDB
  * choose the copy whose heavy-atom centroid is closest to the Vina grid centre
  * write a trimmed PDB (all protein ATOM records + only that ligand copy as HETATM)
    named "<Target-with-dashes>_<LIG>.pdb" (DrugCLIP encodes pocket residues within 6 A)
  * record ligand centroid and offset vs grid centre (sanity)

Also writes ligand CSVs (5 actives + 5 decoys from the frozen subsets) for the two
sanity-check targets (LpxC_Ec, DHFR_Ec).
"""
import csv
import json
import os
import sys

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

ROOT = r"D:\harness\projects\d1_gramneg_vs"
OUT = os.path.join(ROOT, "results", "trackA", "method_axis")  # read-only source
R2 = os.path.join(ROOT, "results", "tool_feasibility_round2")
POCKETS = os.path.join(R2, "drugclip_pockets")
CSV_OUT = os.path.join(R2, "ligand_sets")
MANIFEST = os.path.join(ROOT, "data", "receptors", "receptor_manifest.csv")
SUBSETS = os.path.join(OUT, "subsets_diffdock.csv")
TARGETS = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec",
           "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]


def dash(t):
    return t.replace("_", "-")


class ProteinPlusLigand(Select):
    def __init__(self, keep_residue):
        self.keep = keep_residue

    def accept_residue(self, res):
        if res.id[0] == " ":
            return True
        return res is self.keep


def main():
    os.makedirs(POCKETS, exist_ok=True)
    os.makedirs(CSV_OUT, exist_ok=True)
    man = list(csv.DictReader(open(MANIFEST, encoding="utf-8-sig")))
    holo = {}
    for r in man:
        if str(r.get("role")) != "holo":
            continue
        holo[r["target_key"]] = r
    parser = PDBParser(QUIET=True)
    report = {}
    for t in TARGETS:
        src = holo.get("GyrB_Ec") if t == "Gyrase_Ec" else holo[t]
        pdb_path = os.path.join(ROOT, src["pdb_raw"])
        lig_code = src["ligand_comp_id"]
        center = np.array([float(src["grid_center_x"]), float(src["grid_center_y"]),
                           float(src["grid_center_z"])])
        stem = "".join(
            ch for ch in os.path.basename(pdb_path) if ch.isalnum() or ch in "._-"
        ).replace(".pdb", "")
        structure = parser.get_structure(stem, pdb_path)
        model = structure[0]
        # collect ligand copies
        copies = []
        for chain in model:
            for res in chain:
                if res.id[0] != " " and res.resname == lig_code:
                    coords = np.array([a.coord for a in res.get_atoms() if a.element != "H"])
                    if len(coords) == 0:
                        continue
                    copies.append({"chain": chain.id, "resid": res.id, "res": res,
                                   "center": coords.mean(axis=0), "n_atoms": int(len(coords))})
        if not copies:
            report[t] = {"status": "no_ligand_copy", "ligand_code": lig_code, "pdb": src["pdb_raw"]}
            print(t, "NO LIGAND COPY", lig_code, src["pdb_raw"])
            continue
        dists = [float(np.linalg.norm(c["center"] - center)) for c in copies]
        best = copies[int(np.argmin(dists))]
        out_pdb = os.path.join(POCKETS, f"{dash(t)}_{lig_code}.pdb")
        io = PDBIO()
        io.set_structure(structure)
        io.save(out_pdb, ProteinPlusLigand(best["res"]))
        report[t] = {
            "status": "ok", "ligand_code": lig_code, "pdb_raw": src["pdb_raw"],
            "pocket_pdb": os.path.relpath(out_pdb, ROOT),
            "pdb_id": src.get("pdb_id", ""),
            "grid_center": center.tolist(),
            "ligand_centroid": best["center"].tolist(),
            "ligand_copy": f"chain {best['chain']} res {best['resid'][1]}",
            "offset_A": round(float(np.linalg.norm(best["center"] - center)), 3),
            "n_ligand_heavy_atoms": best["n_atoms"],
            "all_ligand_copies": [{"chain": c["chain"], "resid": str(c["resid"][1]),
                                   "center": c["center"].tolist(),
                                   "offset_A": round(float(np.linalg.norm(c["center"] - center)), 3)}
                                  for c in copies],
        }
        print(f"{t:10s} {lig_code} pocket={os.path.basename(out_pdb):22s} "
              f"offset={report[t]['offset_A']:.3f} A copies={len(copies)}")
    json.dump(report, open(os.path.join(R2, "drugclip_pocket_prep.json"), "w"),
              indent=1, ensure_ascii=False)
    # ligand CSVs: 5 actives + 5 decoys per sanity target (frozen subsets)
    subs = list(csv.DictReader(open(SUBSETS, encoding="utf-8-sig")))
    for t in ["LpxC_Ec", "DHFR_Ec"]:
        rows = [r for r in subs if r["target"] == t]
        act = [r for r in rows if r["label"] == "1"][:5]
        dec = [r for r in rows if r["label"] == "0"][:5]
        sel = act + dec
        out_csv = os.path.join(CSV_OUT, f"{t}_10lig.csv")
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["mol_id", "smiles"])
            for r in sel:
                w.writerow([r["ligand_id"], r["smiles"]])
        print("wrote", out_csv, len(sel), "ligands")
    print("done; report:", os.path.join(R2, "drugclip_pocket_prep.json"))


if __name__ == "__main__":
    main()
