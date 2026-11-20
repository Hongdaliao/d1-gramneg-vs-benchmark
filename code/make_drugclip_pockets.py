#!/usr/bin/env python
"""Build DrugCLIP pocket inputs for the 9 Track A targets (Amendment #13, prep only).

Scheme 1 (pre-registered): use each target's native co-crystal ligand copy whose
heavy-atom centroid matches the Vina grid centre (Track A grid centres were derived
from those ligand centres; observed offsets 0.000-0.001 A).

For each target:
  * copy from data/receptors/raw/<pdb_id>.pdb the protein ATOM records + ONLY the chosen
    native ligand copy (resname == manifest ligand_comp_id) as HETATM
  * name it "<Target-with-dashes>_<LIG>.pdb" (DrugCLIP rule: text between first '_' and
    '.pdb' is the ligand residue name; target part must not contain underscores)
  * record ligand code, heavy atom count, 6 A pocket atom/residue counts, centroid,
    grid-centre offset, sha256

Outputs (prep only, NO screening):
  results/trackA/method_axis/drugclip/pockets/<Target-dash>_<LIG>.pdb
  results/trackA/method_axis/drugclip/pocket_manifest.json
  results/trackA/method_axis/drugclip/pocket_manifest.csv
"""
import csv
import hashlib
import json
import os

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

ROOT = r"D:\harness\projects\d1_gramneg_vs"
OUT = os.path.join(ROOT, "results", "trackA", "method_axis", "drugclip")
POCKETS = os.path.join(OUT, "pockets")
MANIFEST = os.path.join(ROOT, "data", "receptors", "receptor_manifest.csv")
TARGETS = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec",
           "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]


def dash(t):
    return t.replace("_", "-")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class ProteinPlusLigand(Select):
    def __init__(self, keep_residue):
        self.keep = keep_residue

    def accept_residue(self, res):
        if res.id[0] == " ":
            return True
        return res is self.keep


def sixA_counts(prot_coords, lig_coords):
    atoms = 0
    residues = set()
    for a in prot_coords:
        d = np.min(np.linalg.norm(lig_coords - np.array(a["coord"]), axis=1))
        if d <= 6.0:
            atoms += 1
            residues.add((a["chain"], a["resid"]))
    return atoms, len(residues)


def main():
    os.makedirs(POCKETS, exist_ok=True)
    man = list(csv.DictReader(open(MANIFEST, encoding="utf-8-sig")))
    holo = {r["target_key"]: r for r in man if str(r.get("role")) == "holo"}
    parser = PDBParser(QUIET=True)
    report = {}
    rows = []
    for t in TARGETS:
        src = holo.get("GyrB_Ec") if t == "Gyrase_Ec" else holo[t]
        pdb_path = os.path.join(ROOT, src["pdb_raw"])
        lig_code = src["ligand_comp_id"]
        center = np.array([float(src["grid_center_x"]), float(src["grid_center_y"]),
                           float(src["grid_center_z"])])
        stem = "".join(ch for ch in os.path.basename(pdb_path)
                       if ch.isalnum() or ch in "._-").replace(".pdb", "")
        structure = parser.get_structure(stem, pdb_path)
        model = structure[0]
        copies = []
        for chain in model:
            for res in chain:
                if res.id[0] != " " and res.resname == lig_code:
                    coords = np.array([a.coord for a in res.get_atoms() if a.element != "H"])
                    if len(coords):
                        copies.append({"chain": chain.id, "resid": res.id, "res": res,
                                       "center": coords.mean(axis=0), "n_atoms": int(len(coords))})
        if not copies:
            report[t] = {"status": "no_ligand_copy", "ligand_code": lig_code}
            continue
        pick = min(copies, key=lambda c: float(np.linalg.norm(c["center"] - center)))
        out_pdb = os.path.join(POCKETS, f"{dash(t)}_{lig_code}.pdb")
        io = PDBIO()
        io.set_structure(structure)
        io.save(out_pdb, ProteinPlusLigand(pick["res"]))
        # pocket stats from written file (6 A around chosen ligand)
        s2 = parser.get_structure("p", out_pdb)
        prot = [{"coord": a.coord, "chain": a.get_parent().get_parent().id,
                 "resid": a.get_parent().id[1]} for a in s2.get_atoms()
                if a.get_parent().id[0] == " " and a.element != "H"]
        lig = np.array([a.coord for a in s2.get_atoms()
                        if a.get_parent().id[0] != " " and a.element != "H"])
        n_atoms, n_res = sixA_counts(prot, lig)
        rec = {
            "status": "ok", "target": t, "pdb_id": src.get("pdb_id", ""),
            "ligand_code": lig_code, "pocket_pdb": os.path.relpath(out_pdb, ROOT),
            "pocket_sha256": sha256_file(out_pdb),
            "pdb_raw": src["pdb_raw"],
            "grid_center": center.tolist(),
            "ligand_centroid": pick["center"].tolist(),
            "offset_A": round(float(np.linalg.norm(pick["center"] - center)), 4),
            "ligand_copy": f"chain {pick['chain']} res {pick['resid'][1]}",
            "n_ligand_heavy_atoms": pick["n_atoms"],
            "pocket_6A_atoms": n_atoms, "pocket_6A_residues": n_res,
            "n_protein_heavy_atoms_in_pdb": len(prot),
            "n_ligand_copies": len(copies),
        }
        report[t] = rec
        rows.append({k: (json.dumps(v) if isinstance(v, (list, dict)) else v)
                     for k, v in rec.items()})
        print(f"{t:10s} {lig_code:4s} pocket6A atoms={n_atoms:4d} res={n_res:3d} "
              f"offset={rec['offset_A']:.4f} A")
    meta = {
        "amendment": "#13 DrugCLIP pocket scheme 1 (native co-crystal ligand)",
        "rule": "DrugCLIP encodes standard-AA residues within 6 A of the ligand residue named "
                "in the PDB filename (text between first '_' and '.pdb'); target part uses dashes",
        "targets": report,
    }
    json.dump(meta, open(os.path.join(OUT, "pocket_manifest.json"), "w"), indent=1,
              ensure_ascii=False)
    with open(os.path.join(OUT, "pocket_manifest.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote:", os.path.join(OUT, "pocket_manifest.json"),
          "| json sha256:", sha256_file(os.path.join(OUT, "pocket_manifest.json")))


if __name__ == "__main__":
    main()
