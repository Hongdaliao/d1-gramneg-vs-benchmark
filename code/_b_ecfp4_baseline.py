#!/usr/bin/env python
"""ECFP4+Tanimoto 1-NN baseline for Track A (task-30, CPU-only).

Score each ligand by max Tanimoto similarity (ECFP4, 2048 bits) to the TRAIN actives
(official temporal split; actives use leave-one-out), then compute AUC/EF5/EF1 vs labels.
Outputs: results/trackA/ecfp4_baseline/<target>.csv (ligand_id,label,score,status)
         results/trackA/ecfp4_baseline/baseline_metrics.csv
Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\_b_ecfp4_baseline.py
"""
import csv, json, math, os, sys

ROOT = r"D:\harness\projects\d1_gramneg_vs"
TRACKA = os.path.join(ROOT, "results", "trackA", "full")
SPLITS = os.path.join(ROOT, "data", "trackA", "splits")
OUT = os.path.join(ROOT, "results", "trackA", "ecfp4_baseline")
os.makedirs(OUT, exist_ok=True)
TARGETS = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec", "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

def fp(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)

def metrics(data):
    """data: list of (score,label); higher score = better (NN similarity). Pair-count AUC."""
    data = [(s, l) for s, l in data if s is not None]
    data.sort(key=lambda x: -x[0])
    N = len(data); nact = sum(1 for s, l in data if l == 1); nneg = N - nact
    if nact == 0 or nneg == 0:
        return {"N": N, "n_act": nact, "auc": None, "ef1": None, "ef5": None, "bedroc20": None}
    pos = [s for s, l in data if l == 1]
    neg = [s for s, l in data if l == 0]
    wins = 0.0
    for ps in pos:
        for ns in neg:
            if ps > ns:
                wins += 1.0
            elif ps == ns:
                wins += 0.5
    auc = wins / (nact * nneg)
    def ef(f):
        k = max(1, int(f * N))
        hits = sum(1 for x in data[:k] if x[1] == 1)
        return round((hits / nact) / (k / N), 4)
    alpha = 20.0; ra = nact / N
    s = sum(math.exp(-alpha * i / N) for i, x in enumerate(data) if x[1] == 1)
    denom = (1 - math.exp(-alpha)) / (math.exp(alpha / N) - 1) + ra * alpha / (1 - math.exp(-alpha))
    return {"N": N, "n_act": nact, "n_neg": nneg, "auc": round(auc, 4), "ef5": ef(0.05), "ef1": ef(0.01),
            "bedroc20": round(s / (nact * denom), 4)}

def main():
    rows_metrics = []
    for t in TARGETS:
        d = os.path.join(TRACKA, t)
        man = os.path.join(d, "manifest.csv")
        if not (os.path.exists(man) and os.path.exists(os.path.join(d, f"{t}_scores.csv"))):
            print(t, "skip (no manifest/scores)"); continue
        manrows = list(csv.DictReader(open(man, encoding="utf-8-sig")))
        # train actives from official split file
        sp = os.path.join(SPLITS, f"{t}_temporal_latest20.csv")
        train_smis = []
        if os.path.exists(sp):
            for r in csv.DictReader(open(sp, encoding="utf-8-sig")):
                if (r.get("split_latest20") or "") == "train" and r.get("canonical_smiles"):
                    train_smis.append(r["canonical_smiles"])
        train_fps = [f for f in (fp(s) for s in train_smis) if f is not None]
        if not train_fps:
            print(t, "no train actives -> skip"); continue
        # score all manifest ligands
        out = []
        for r in manrows:
            f = fp(r.get("smiles", ""))
            if f is None:
                out.append({"ligand_id": r["ligand_id"], "label": r["label"], "score": "", "status": "smiles_invalid"})
                continue
            sims = DataStructs.BulkTanimotoSimilarity(f, train_fps)
            # actives: exclude self-similarity if identical fingerprint (leave-one-out approx)
            s = max(sims) if sims else 0.0
            if r["label"] == "1":
                # remove one maximal value if it corresponds to self (exact fp match)
                same = sum(1 for tf in train_fps if tf == f)
                if same > 0:
                    srt = sorted(sims, reverse=True)
                    srt = srt[same:] or sorted(sims, reverse=True)
                    s = srt[0] if srt else 0.0
            out.append({"ligand_id": r["ligand_id"], "label": r["label"], "score": round(float(s), 4), "status": "ok"})
        p = os.path.join(OUT, f"{t}.csv")
        with open(p, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=["ligand_id", "label", "score", "status"]); w.writeheader(); w.writerows(out)
        m = metrics([(float(r["score"]), 1 if str(r["label"]) == "1" else 0) for r in out if r["status"] == "ok"])
        m.update({"target": t})
        rows_metrics.append(m)
        print(t, json.dumps(m), flush=True)
    pm = os.path.join(OUT, "baseline_metrics.csv")
    cols = ["target", "N", "n_act", "n_neg", "auc", "ef5", "ef1", "bedroc20"]
    with open(pm, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows_metrics)
    print("wrote", pm)

if __name__ == "__main__":
    main()
