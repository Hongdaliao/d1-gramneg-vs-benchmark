#!/usr/bin/env python
"""Post-process DrugCLIP retrieval output (mol_id,score txt) into ranked CSV with SMILES."""
import csv
import pickle
import sys

import lmdb


def main(mol_lmdb, txt_path, out_csv):
    smiles = {}
    env = lmdb.open(mol_lmdb, subdir=False, readonly=True, lock=False, readahead=False, meminit=False)
    with env.begin() as txn:
        for key, value in txn.cursor():
            try:
                d = pickle.loads(value)
                mid = d.get("smi", "")
                smi = d.get("_smiles", d.get("smi", ""))
                if mid and smi:
                    smiles[mid] = smi
            except Exception:
                pass
    env.close()
    results = []
    with open(txt_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.rsplit(",", 1)
            if len(parts) == 2:
                try:
                    results.append((float(parts[1]), parts[0]))
                except ValueError:
                    continue
    results.sort(key=lambda x: x[0], reverse=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["rank", "mol_id", "drugclip_score", "smiles"])
        for i, (score, mid) in enumerate(results):
            w.writerow([i + 1, mid, f"{score:.4f}", smiles.get(mid, "N/A")])
    print(f"{len(results)} results -> {out_csv}")
    for score, mid in results[:10]:
        print(f"  #{results.index((score, mid))+1:2d} {mid:>16s} {score:7.4f}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: post.py MOL_LMDB RETRIEVAL_TXT OUT_CSV")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
