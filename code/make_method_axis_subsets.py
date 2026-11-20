#!/usr/bin/env python
"""Build frozen method-axis subsets for Amendment #12 (DiffDock / Boltz-2).

READ-ONLY over existing Track A outputs. Does NOT run any docking / prediction.
Outputs (only when executed; Lead instructs):
  results/trackA/method_axis/subsets_diffdock.csv
  results/trackA/method_axis/subsets_boltz2.csv
  results/trackA/method_axis/subsets_sha256.json   (CSV sha256 + input-file hashes + script sha)

Sampling protocol (pre-registered; see notes/preregistration_amendment12_method_axis_draft.md):
  * eligible actives = label==1, status=="ok", box_clipped != True,
    has >=1 matched decoy that is itself status=="ok"
  * per target: rng = random.Random(42) (fresh per target, iterated in fixed TARGETS order)
    DiffDock: k = min(50, n_eligible); sample k actives WITHOUT replacement from the
    ligand_id-sorted eligible list (sampling order retained).
    For each sampled active in sampling order: matched_pool = sorted ok matched decoys
    (dec_<inchikey>) not yet used; pick rng.choice(matched_pool);
    if empty -> fallback: rng.choice(sorted unused eligible decoys matched to any active
    of that target), decoy_match="fallback"; if still empty -> decoy_match="none".
    Nested Boltz-2 subset = first 10 actives (sampling order) of the DiffDock subset and
    their assigned decoys (first 10 decoy assignments); if fewer, take all.
  * target order fixed TARGETS; seed=42; no post-hoc change.
Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\make_method_axis_subsets.py
"""
import csv
import hashlib
import json
import os
import random
import sys

ROOT = r"D:\harness\projects\d1_gramneg_vs"
OUTDIR = os.path.join(ROOT, "results", "trackA", "method_axis")
SCORES = os.path.join(ROOT, "results", "trackA", "full")
DECOYS = os.path.join(ROOT, "data", "decoys")
ACTIVES = os.path.join(ROOT, "data", "trackA")
RECEPTORS = os.path.join(ROOT, "data", "receptors")

TARGETS = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec",
           "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]
SEED = 42
N_DIFFDOCK = 50          # actives per target (target total = 2x)
N_BOLTZ = 10             # actives per target (target total = 2x)
RECEPTOR_FILE = {
    "LpxC_Ec": "LpxC_Ec_holo.pdbqt", "LpxC_Pa": "LpxC_Pa_holo.pdbqt",
    "DHFR_Ec": "DHFR_Ec_holo_3DAU.pdbqt", "GyrB_Ec": "GyrB_Ec_holo.pdbqt",
    "Gyrase_Ec": "GyrB_Ec_holo.pdbqt", "GyrA_Ec": "GyrA_Ec_holo.pdbqt",
    "FabI_Ec": "FabI_Ec_holo.pdbqt", "TEM1_Ec": "TEM1_Ec_holo.pdbqt",
    "LepB_Ec": "LepB_Ec_holo.pdbqt",
}
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()

def read_csv(path):
    if not os.path.exists(path):
        return []
    return list(csv.DictReader(open(path, encoding="utf-8-sig")))

def read_actives(target):
    """compound_chembl_id -> inchikey."""
    out = {}
    for a in read_csv(os.path.join(ACTIVES, f"{target}_actives.csv")):
        cid = (a.get("compound_chembl_id") or "").strip()
        ik = (a.get("inchikey") or "").strip()
        if cid and ik:
            out[cid] = ik
    return out

def read_decoy_pool(target):
    """matched_active_inchikey -> sorted set of decoy inchikeys (from s1 pool)."""
    pool = {}
    for r in read_csv(os.path.join(DECOYS, f"{target}_decoys_s1.csv")):
        ik = (r.get("decoy_inchikey") or "").strip()
        ma = (r.get("matched_active_inchikey") or "").strip()
        if ik and ma:
            pool.setdefault(ma, set()).add(ik)
    return {k: sorted(v) for k, v in pool.items()}

def protein_sequence_guess(pdbqt_path):
    """Derive a single-chain protein sequence from a prepared receptor PDBQT.

    MGLTools PDBQT: ATOM records carry residue names; chain id may be blank.
    Returns (sequence, chain_id, n_residues, note).
    """
    if not os.path.exists(pdbqt_path):
        return "", "", 0, "receptor_missing"
    chains = {}
    for ln in open(pdbqt_path, encoding="utf-8", errors="replace"):
        if not ln.startswith("ATOM"):
            continue
        resn = ln[17:20].strip()
        if resn not in AA3:
            continue
        ch = ln[21] if len(ln) > 21 else " "
        try:
            resseq = int(ln[22:26])
        except Exception:
            continue
        chains.setdefault(ch, {})[(resseq, ln[26])] = AA3[resn]
    if not chains:
        return "", "", 0, "no_standard_residues"
    # choose the chain with most residues (document choice in column)
    ch = max(chains, key=lambda c: len(chains[c]))
    seq = "".join(chains[ch][k] for k in sorted(chains[ch]))
    note = f"single chain '{ch}' from prepared PDBQT; multichain={len(chains) > 1}; " \
           f"other_chains=" + ",".join(sorted(c for c in chains if c != ch))
    return seq, ch, len(seq), note

def build():
    os.makedirs(OUTDIR, exist_ok=True)
    inputs = {}
    diffdock, boltz2 = [], []
    summary = {"seed": SEED, "n_diffdock_actives": N_DIFFDOCK, "n_boltz_actives": N_BOLTZ,
               "targets": {}, "algorithm": __doc__.split("Sampling protocol", 1)[-1][:1200]}
    for target in TARGETS:
        d = os.path.join(SCORES, target)
        scores_path = os.path.join(d, f"{target}_scores.csv")
        manifest_path = os.path.join(d, "manifest.csv")
        actives_path = os.path.join(ACTIVES, f"{target}_actives.csv")
        decoys_path = os.path.join(DECOYS, f"{target}_decoys_s1.csv")
        for p in [scores_path, manifest_path, actives_path, decoys_path]:
            inputs.setdefault(target, {})[os.path.relpath(p, ROOT)] = sha256_file(p)
        scores = {r["ligand_id"]: r for r in read_csv(scores_path)}
        man = {r["ligand_id"]: r for r in read_csv(manifest_path)}
        act_ik = read_actives(target)
        pool = read_decoy_pool(target)
        ok_dec = {lid for lid, r in scores.items()
                  if str(r.get("label")) != "1" and r.get("status") == "ok"}
        # eligible actives (sorted ligand_id; deterministic universe)
        eligible = []
        for lid in sorted(scores):
            r = scores[lid]
            if str(r.get("label")) != "1" or r.get("status") != "ok":
                continue
            if str(r.get("box_clipped", "False")).lower() == "true":
                continue
            chembl = lid[4:] if lid.startswith("act_") else lid
            ik = act_ik.get(chembl, "")
            matched = [x for x in pool.get(ik, []) if ("dec_" + x) in ok_dec]
            if matched:
                eligible.append({"ligand_id": lid, "inchikey": ik, "matched": matched})
        # deterministic per-target RNG (fresh seed per target; target order fixed)
        rng = random.Random(SEED)
        k = min(N_DIFFDOCK, len(eligible))
        sampled = rng.sample(eligible, k)  # sampling order retained
        used_dec = set()
        all_ok_matched = sorted({x for e in eligible for x in e["matched"]})
        assignments = []
        for rank, e in enumerate(sampled):
            candidates = [x for x in e["matched"] if x not in used_dec]
            match_kind = "matched"
            if not candidates:
                fallback = [x for x in all_ok_matched if x not in used_dec]
                if fallback:
                    candidates, match_kind = fallback, "fallback"
            if candidates:
                pick = rng.choice(candidates)
                used_dec.add(pick)
                assignments.append((e, "dec_" + pick, match_kind))
            else:
                assignments.append((e, "", "none"))
        # write DiffDock rows (active then decoy, deterministic block order)
        for rank, (e, decoy_id, kind) in enumerate(assignments):
            lid = e["ligand_id"]
            m = man.get(lid, {})
            sc = scores.get(lid, {})
            diffdock.append({
                "subset_id": f"{target}_dd_{rank:03d}",
                "method": "diffdock",
                "target": target,
                "ligand_id": lid,
                "label": "1",
                "smiles": m.get("smiles", ""),
                "ligand_source": "trackA_manifest_pdbqt",
                "ligand_pdbqt": m.get("pdbqt", ""),
                "receptor_pdbqt": os.path.join(RECEPTORS, RECEPTOR_FILE[target]),
                "receptor_pdb": os.path.join(OUTDIR, "receptors", f"{target}.pdb"),
                "vina_score": sc.get("score", ""),
                "vina_status": sc.get("status", ""),
                "box_clipped": sc.get("box_clipped", "False"),
                "temporal_split": sc.get("temporal_split", ""),
                "matched_active_inchikey": e["inchikey"],
                "decoy_match": "active",
                "note": "Amendment #12 exploratory subset; actives block",
            })
            if decoy_id:
                dm = man.get(decoy_id, {})
                dsc = scores.get(decoy_id, {})
                diffdock.append({
                    "subset_id": f"{target}_dd_{rank:03d}",
                    "method": "diffdock",
                    "target": target,
                    "ligand_id": decoy_id,
                    "label": "0",
                    "smiles": dm.get("smiles", ""),
                    "ligand_source": "trackA_manifest_pdbqt",
                    "ligand_pdbqt": dm.get("pdbqt", ""),
                    "receptor_pdbqt": os.path.join(RECEPTORS, RECEPTOR_FILE[target]),
                    "receptor_pdb": os.path.join(OUTDIR, "receptors", f"{target}.pdb"),
                    "vina_score": dsc.get("score", ""),
                    "vina_status": dsc.get("status", ""),
                    "box_clipped": dsc.get("box_clipped", "False"),
                    "temporal_split": "decoy",
                    "matched_active_inchikey": e["inchikey"],
                    "decoy_match": kind,
                    "note": "Amendment #12 exploratory subset; decoy block",
                })
        # nested Boltz-2 subset: first N_BOLTZ assignments (actives + assigned decoys)
        nb = min(N_BOLTZ, len(assignments))
        seq, chain, nres, seq_note = protein_sequence_guess(os.path.join(RECEPTORS, RECEPTOR_FILE[target]))
        for rank, (e, decoy_id, kind) in enumerate(assignments[:nb]):
            lid = e["ligand_id"]
            m = man.get(lid, {})
            sc = scores.get(lid, {})
            boltz2.append({
                "subset_id": f"{target}_bz_{rank:03d}", "method": "boltz2", "target": target,
                "ligand_id": lid, "label": "1", "smiles": m.get("smiles", ""),
                "receptor_pdbqt": os.path.join(RECEPTORS, RECEPTOR_FILE[target]),
                "protein_sequence": seq, "receptor_chain": chain, "receptor_n_residues": nres,
                "sequence_note": seq_note,
                "vina_score": sc.get("score", ""), "vina_status": sc.get("status", ""),
                "box_clipped": sc.get("box_clipped", "False"),
                "temporal_split": sc.get("temporal_split", ""),
                "merged_from_diffdock_matched_active": e["inchikey"],
                "decoy_match": "active", "note": "Amendment #12 exploratory subset; Boltz-2 actives block",
            })
            if decoy_id:
                dm = man.get(decoy_id, {})
                dsc = scores.get(decoy_id, {})
                boltz2.append({
                    "subset_id": f"{target}_bz_{rank:03d}", "method": "boltz2", "target": target,
                    "ligand_id": decoy_id, "label": "0", "smiles": dm.get("smiles", ""),
                    "receptor_pdbqt": os.path.join(RECEPTORS, RECEPTOR_FILE[target]),
                    "protein_sequence": seq, "receptor_chain": chain, "receptor_n_residues": nres,
                    "sequence_note": seq_note,
                    "vina_score": dsc.get("score", ""), "vina_status": dsc.get("status", ""),
                    "box_clipped": dsc.get("box_clipped", "False"), "temporal_split": "decoy",
                    "merged_from_diffdock_matched_active": e["inchikey"],
                    "decoy_match": kind, "note": "Amendment #12 exploratory subset; Boltz-2 decoy block",
                })
        summary["targets"][target] = {
            "n_eligible_actives": len(eligible),
            "diffdock_actives": len(sampled),
            "diffdock_decoys": sum(1 for _, dec, _ in assignments if dec),
            "diffdock_decoy_match": {"matched": sum(1 for _, _, kd in assignments if kd == "matched"),
                                     "fallback": sum(1 for _, _, kd in assignments if kd == "fallback"),
                                     "none": sum(1 for _, _, kd in assignments if kd == "none")},
            "boltz2_actives": nb,
            "boltz2_decoys": sum(1 for _, dec, _ in assignments[:nb] if dec),
            "receptor_pdbqt": os.path.join("data", "receptors", RECEPTOR_FILE[target]),
            "receptor_sequence_len_guess": nres,
            "receptor_sequence_note": seq_note,
        }
    # write CSVs
    out = {}
    for name, rows in [("subsets_diffdock.csv", diffdock), ("subsets_boltz2.csv", boltz2)]:
        path = os.path.join(OUTDIR, name)
        if not rows:
            raise RuntimeError(f"empty subset: {name}")
        fields = list(rows[0].keys())
        with open(path, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        out[name] = {"path": path, "rows": len(rows), "sha256": sha256_file(path)}
        print(name, "rows:", len(rows), "sha256:", out[name]["sha256"])
    meta = {
        "amendment": "#12 method-class complementarity axis (post-hoc exploratory, draft)",
        "generated_at": __import__("time").strftime("%Y-%m-%dT%H:%M:%S"),
        "script": os.path.relpath(os.path.abspath(__file__), ROOT),
        "script_sha256": sha256_file(os.path.abspath(__file__)),
        "seed": SEED,
        "targets_order": TARGETS,
        "outputs": out,
        "inputs": inputs,
        "target_summary": summary["targets"],
        "note": "READ-ONLY subset generation; no docking/prediction performed; DiffDock/Boltz-2 runs require Lead go-ahead.",
    }
    with open(os.path.join(OUTDIR, "subsets_sha256.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=1, ensure_ascii=False)
    print("wrote", os.path.join(OUTDIR, "subsets_sha256.json"))

if __name__ == "__main__":
    build()
