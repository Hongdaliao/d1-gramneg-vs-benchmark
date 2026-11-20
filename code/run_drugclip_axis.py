#!/usr/bin/env python
"""DrugCLIP 9-target screening runner (Amendment #13; frozen inputs only).

Per target: copy pocket PDB -> per-target dir; encode pocket (6 folds, DrugCLIP
encode_pockets.py); retrieval over the frozen 806-molecule library; post-process.
Outputs scores in schema: target,ligand_id,label,score,status (+ NA rows), plus a
combined CSV, NA coverage, progress JSON and per-target logs.

Run AFTER Boltz-2 finishes and after Lead freezes Amendments #13/#14:
  <conda-python> code\\trackA\\run_drugclip_axis.py [--targets T1,T2] [--fold-version 6_folds]
"""
import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import time
from datetime import datetime

ROOT = r"D:\harness\projects\d1_gramneg_vs"
B = os.path.join(ROOT, "results", "trackA", "method_axis", "drugclip")
POCKETS = os.path.join(B, "pockets")
LIB_CSV = os.path.join(B, "inputs", "library_806.csv")
LIB_LMDB = os.path.join(B, "inputs", "library_806.lmdb")
SCORES = os.path.join(B, "scores")
SCREEN = os.path.join(B, "screening")
LOGS = os.path.join(B, "logs")
PROGRESS = os.path.join(B, "progress.json")
DRUGCLIP_ROOT = r"F:\Drugclip\Drugclip_run"
ENV_PY = os.path.join(DRUGCLIP_ROOT, "env", "python.exe")
POST = os.path.join(ROOT, "code", "trackA", "feasibility_drugclip_post.py")
TARGETS = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec",
           "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]
POCKET_FILE = {
    "LpxC_Ec": "LpxC-Ec_GVR.pdb", "LpxC_Pa": "LpxC-Pa_EUY.pdb",
    "GyrB_Ec": "GyrB-Ec_4QR.pdb", "Gyrase_Ec": "Gyrase-Ec_4QR.pdb",
    "DHFR_Ec": "DHFR-Ec_MTX.pdb", "GyrA_Ec": "GyrA-Ec_SM8.pdb",
    "FabI_Ec": "FabI-Ec_TCL.pdb", "TEM1_Ec": "TEM1-Ec_FTA.pdb",
    "LepB_Ec": "LepB-Ec_1PN.pdb",
}
SCHEMA = ["target", "ligand_id", "label", "score", "status"]


def now():
    return datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(p):
    if not os.path.exists(p):
        return []
    return list(csv.DictReader(open(p, encoding="utf-8-sig")))


def write_csv(p, rows, fields=SCHEMA):
    tmp = p + ".tmp"
    with open(tmp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, p)


def load_progress():
    if os.path.exists(PROGRESS):
        try:
            return json.load(open(PROGRESS, encoding="utf-8"))
        except Exception:
            pass
    return {"amendment": "#13 frozen", "started_at": now(), "updated_at": now(),
            "targets": {}, "fold_version": "6_folds"}


def save_progress(p):
    p["updated_at"] = now()
    tmp = PROGRESS + ".tmp"
    json.dump(p, open(tmp, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    os.replace(tmp, PROGRESS)


def run(cmd, log_path, cwd):
    t0 = time.time()
    with open(log_path, "a", encoding="utf-8", errors="replace") as lf:
        lf.write(f"\n===== {now()} {' '.join(cmd)[:300]} =====\n")
        r = subprocess.run(cmd, cwd=cwd, stdout=lf, stderr=subprocess.STDOUT, timeout=3 * 3600)
    return r.returncode, round(time.time() - t0, 1)


def post_process(target, lmdb, txt, out_csv):
    smiles = {}
    import lmdb as _lmdb
    import pickle
    env = _lmdb.open(lmdb, subdir=False, readonly=True, lock=False, readahead=False, meminit=False)
    with env.begin() as txn:
        for key, value in txn.cursor():
            try:
                d = pickle.loads(value)
                mid = d.get("smi", "")
                if mid:
                    smiles[mid] = d.get("_smiles", mid)
            except Exception:
                pass
    env.close()
    results = []
    with open(txt, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.strip().rsplit(",", 1)
            if len(parts) == 2:
                try:
                    results.append((float(parts[1]), parts[0]))
                except ValueError:
                    continue
    results.sort(key=lambda x: x[0], reverse=True)
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["rank", "mol_id", "drugclip_score", "smiles"])
        for i, (s, m) in enumerate(results):
            w.writerow([i + 1, m, f"{s:.4f}", smiles.get(m, "N/A")])
    return len(results)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default=None)
    ap.add_argument("--fold-version", default="6_folds")
    args = ap.parse_args()
    targets = [x.strip() for x in args.targets.split(",")] if args.targets else TARGETS
    os.makedirs(SCORES, exist_ok=True)
    os.makedirs(SCREEN, exist_ok=True)
    os.makedirs(LOGS, exist_ok=True)
    lib = read_csv(LIB_CSV)
    by_mol = {r["mol_id"]: r for r in lib}
    prog = load_progress()
    for t in targets:
        t0 = time.time()
        tdir = os.path.join(SCREEN, t)
        os.makedirs(tdir, exist_ok=True)
        pocket = os.path.join(tdir, POCKET_FILE[t])
        if not os.path.exists(pocket):
            shutil.copy2(os.path.join(POCKETS, POCKET_FILE[t]), pocket)
        log = os.path.join(LOGS, f"{t}.log")
        # encode pocket
        enc_cmd = [ENV_PY, os.path.join(DRUGCLIP_ROOT, "unimol", "encode_pockets.py"),
                   "--user-dir", "./unimol", "./dict", "--valid-subset", "test",
                   "--num-workers", "0", "--ddp-backend=c10d", "--batch-size", "8",
                   "--task", "drugclip", "--loss", "in_batch_softmax", "--arch", "drugclip",
                   "--max-pocket-atoms", "511", "--seed", "1", "--log-interval", "100",
                   "--log-format", "simple", "--pocket-dir", tdir, "--fp16"]
        rc, sec = run(enc_cmd, log, DRUGCLIP_ROOT)
        if rc != 0 or not os.path.exists(os.path.join(tdir, "pocket.lmdb")):
            prog["targets"][t] = {"status": "encode_failed", "rc": rc, "sec": sec, "at": now()}
            save_progress(prog)
            continue
        # retrieval
        txt = os.path.join(tdir, "retrieval.txt")
        raw_csv = os.path.join(tdir, "raw_scores.csv")
        ret_cmd = [ENV_PY, os.path.join(DRUGCLIP_ROOT, "unimol", "retrieval.py"),
                   "--user-dir", "./unimol", "./dict", "--valid-subset", "test",
                   "--num-workers", "0", "--ddp-backend=c10d", "--batch-size", "8",
                   "--task", "drugclip", "--loss", "in_batch_softmax", "--arch", "drugclip",
                   "--max-pocket-atoms", "511", "--fp16", "--fp16-init-scale", "4",
                   "--fp16-scale-window", "256", "--seed", "1", "--log-interval", "100",
                   "--log-format", "simple", "--mol-path", LIB_LMDB,
                   "--pocket-path", os.path.join(tdir, "pocket.lmdb"),
                   "--fold-version", args.fold_version, "--use-cache", "False",
                   "--save-path", txt]
        rc2, sec2 = run(ret_cmd, log, DRUGCLIP_ROOT)
        n = post_process(t, LIB_LMDB, txt, raw_csv) if rc2 == 0 else 0
        raw = {r["mol_id"]: r for r in read_csv(raw_csv)}
        rows = []
        for r in lib:
            if r["target"] != t:
                continue
            mid = r["mol_id"]
            if mid in raw:
                rows.append({"target": t, "ligand_id": r["ligand_id"], "label": r["label"],
                             "score": raw[mid]["drugclip_score"], "status": "ok"})
            else:
                rows.append({"target": t, "ligand_id": r["ligand_id"], "label": r["label"],
                             "score": "", "status": "NA"})
        n_ok = sum(1 for x in rows if x["status"] == "ok")
        out = os.path.join(SCORES, f"{t}_{args.fold_version}.csv")
        write_csv(out, rows)
        prog["targets"][t] = {"status": "ok" if rc2 == 0 else "retrieval_failed",
                              "n": len(rows), "n_ok": n_ok, "n_na": len(rows) - n_ok,
                              "encode_sec": sec, "retrieval_sec": sec2,
                              "scores_csv": os.path.relpath(out, ROOT),
                              "scores_sha256": sha256(out), "at": now(),
                              "total_sec": round(time.time() - t0, 1)}
        save_progress(prog)
    # aggregate
    all_rows = []
    for t in targets:
        p = os.path.join(SCORES, f"{t}_{args.fold_version}.csv")
        for i, r in enumerate(read_csv(p)):
            r["rank_within_target"] = i + 1
            all_rows.append(r)
    agg = os.path.join(B, f"drugclip_scores_{args.fold_version}.csv")
    if all_rows:
        write_csv(agg, all_rows, SCHEMA + ["rank_within_target"])
    na = [r for r in all_rows if r["status"] != "ok"]
    na_csv = os.path.join(B, f"drugclip_na_coverage_{args.fold_version}.csv")
    write_csv(na_csv, na, ["target", "ligand_id", "label", "status"])
    prog["aggregate"] = {"csv": os.path.relpath(agg, ROOT), "sha256": sha256(agg) if all_rows else None,
                         "n": len(all_rows), "n_ok": len(all_rows) - len(na), "n_na": len(na),
                         "na_csv": os.path.relpath(na_csv, ROOT), "at": now()}
    save_progress(prog)
    print(json.dumps(prog.get("aggregate"), indent=1))


if __name__ == "__main__":
    main()
