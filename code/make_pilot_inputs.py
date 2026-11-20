#!/usr/bin/env python
"""G1 LpxC pilot input builder (B-line).

Samples 20 actives + 200 decoys for LpxC_Ec, converts SMILES -> PDBQT via meeko
(mk_prepare_ligand), and derives the docking box from the holo ligand (GVR in 4ISA).

Run with the rdkit-env python (has meeko 0.7.1):
  D:\\anconda\\envs\\rdkit-env\\python.exe code\\trackA\\make_pilot_inputs.py

Outputs (results/trackA/pilot_lpxc/):
  pilot_manifest.csv   ligand_id,label,pdbqt,source
  box.json             center+size (from 4ISA GVR + 4 A margin)
  protocol.md          sampling seed, source files, hashes, box rule
"""
import csv, hashlib, json, os, random, subprocess, sys, datetime

ROOT = r"D:\harness\projects\d1_gramneg_vs"
ACT = os.path.join(ROOT, "data", "trackA", "LpxC_Ec_actives.csv")
DEC = os.path.join(ROOT, "data", "decoys", "LpxC_Ec_decoys_s1.csv")
RAWPDB = os.path.join(ROOT, "data", "receptors", "raw", "4ISA.pdb")
OUT = os.path.join(ROOT, "results", "trackA", "pilot_lpxc")
LIG = os.path.join(OUT, "ligands")
MKPREP = r"D:\harness\projects\d1_gramneg_vs\tools\envs\boltz\Scripts\mk_prepare_ligand.exe"
SEED = 20261116
N_ACT, N_DEC = 20, 200

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def box_from_ligand(pdb, resname="GVR", margin=4.0, lo=18.0, hi=26.0):
    xs, ys, zs = [], [], []
    with open(pdb, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(("HETATM", "ATOM")) and line[17:20].strip() == resname:
                xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
    if not xs:
        raise SystemExit(f"ligand {resname} not found in {pdb}")
    center = [round((min(v) + max(v)) / 2, 2) for v in (xs, ys, zs)]
    size = [round(min(hi, max(lo, (max(v) - min(v)) + 2 * margin)), 1) for v in (xs, ys, zs)]
    return {"center": center, "size": size,
            "rule": f"{resname} heavy-atom bbox center; size = bbox + {margin} A margin, clamped [{lo},{hi}]",
            "n_atoms": len(xs)}

def main():
    os.makedirs(LIG, exist_ok=True)
    rng = random.Random(SEED)
    with open(ACT, newline="", encoding="utf-8-sig") as f:
        act = [r for r in csv.DictReader(f) if (r.get("canonical_smiles") or "").strip()]
    with open(DEC, newline="", encoding="utf-8-sig") as f:
        dec = [r for r in csv.DictReader(f) if (r.get("decoy_canonical_smiles") or "").strip()]
    print(f"actives available={len(act)} decoys available={len(dec)}", flush=True)
    act_s = rng.sample(act, min(N_ACT, len(act)))
    dec_s = rng.sample(dec, min(N_DEC, len(dec)))
    box = box_from_ligand(RAWPDB)
    rows = []
    def convert(lid, smiles):
        sdf = os.path.join(LIG, lid + ".sdf")
        pdbqt = os.path.join(LIG, lid + ".pdbqt")
        if os.path.exists(pdbqt) and os.path.getsize(pdbqt) > 100:
            return pdbqt, "cached"
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            m = Chem.MolFromSmiles(smiles)
            if m is None:
                return None, "bad smiles"
            m = Chem.AddHs(m)
            params = AllChem.ETKDGv3(); params.randomSeed = 0xF00D
            if AllChem.EmbedMolecule(m, params) != 0:
                AllChem.EmbedMolecule(m, useRandomCoords=True, randomSeed=0xF00D)
            try:
                AllChem.MMFFOptimizeMolecule(m, maxIters=200)
            except Exception:
                pass
            Chem.MolToMolFile(m, sdf)
        except Exception as e:
            return None, f"rdkit {type(e).__name__}: {e}"
        res = subprocess.run([MKPREP, "-i", sdf, "-o", pdbqt], capture_output=True, text=True,
                             timeout=300, cwd=LIG)
        if res.returncode != 0 or not os.path.exists(pdbqt):
            return None, (res.stderr or res.stdout or "")[-200:]
        return pdbqt, "built"
    ok = {"act": 0, "dec": 0}; fail = []
    for i, r in enumerate(act_s, 1):
        lid = f"act_{i:03d}_{r['compound_chembl_id']}"
        p, st = convert(lid, r["canonical_smiles"].strip())
        if p: ok["act"] += 1; rows.append({"ligand_id": lid, "label": 1, "pdbqt": p,
                                          "source": "ChEMBL actives LpxC_Ec " + r["compound_chembl_id"]})
        else: fail.append((lid, st))
    for i, r in enumerate(dec_s, 1):
        lid = f"dec_{i:03d}_{r['decoy_chembl_id']}"
        p, st = convert(lid, r["decoy_canonical_smiles"].strip())
        if p: ok["dec"] += 1; rows.append({"ligand_id": lid, "label": 0, "pdbqt": p,
                                          "source": "property-matched decoys s1 " + r["decoy_chembl_id"]})
        else: fail.append((lid, st))
    man = os.path.join(OUT, "pilot_manifest.csv")
    with open(man, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "pdbqt", "source"])
        w.writeheader(); w.writerows(rows)
    with open(os.path.join(OUT, "box.json"), "w", encoding="utf-8") as f:
        json.dump(box, f, indent=1, ensure_ascii=False)
    proto = f"""# G1 LpxC_Ec pilot protocol (B-line, generated {datetime.datetime.now().isoformat(timespec='seconds')})

- sampling seed: {SEED}; actives n={len(act_s)} (of {len(act)}), decoys n={len(dec_s)} (of {len(dec)})
- actives source: data/trackA/LpxC_Ec_actives.csv (sha256={sha256(ACT)[:16]}...)
- decoys source : data/decoys/LpxC_Ec_decoys_s1.csv (sha256={sha256(DEC)[:16]}...)
- receptor: data/receptors/LpxC_Ec_holo.pdbqt (4ISA holo, GVR ligand; A-line prepared)
- box: {json.dumps(box)}
- ligand prep: meeko 0.7.1 mk_prepare_ligand (rdkit-env), 1 pose per ligand
- docking: Vina 1.2.5, exhaustiveness 8, seed 42, num_modes 9, cpu 16
- metrics: EF1%, EF5%, AUC, BEDROC20 + bootstrap CI (1000x, seed 42)
- ligand conversion: ok actives={ok['act']} decoys={ok['dec']}, failures={len(fail)}
"""
    with open(os.path.join(OUT, "protocol.md"), "w", encoding="utf-8") as f:
        f.write(proto)
    print(json.dumps({"ok": ok, "fail": len(fail), "fail_examples": fail[:3],
                      "manifest": man, "box": box}, indent=1, ensure_ascii=False), flush=True)

if __name__ == "__main__":
    main()
