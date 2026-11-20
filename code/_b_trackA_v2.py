#!/usr/bin/env python
"""Track A full docking driver v2 (amendment #7 protocol).

Protocol (frozen-by-lead 2026-11-17):
  * actives: all
  * decoys: s1 rows de-duplicated by decoy_inchikey, then per-active sample of `--ratio`
    property-matched decoys (matched_active_inchikey), deterministic seed 42
    -> decoy_subsample_sha256 recorded
  * Vina: exhaustiveness (main panel 4; non-EF1 targets may use 2), 16 workers x cpu1,
    num_modes=3, seed=42, box 22^3 at A-line center
  * outputs: results/trackA/full/<target>_scores.csv, <target>/summary.json,
    results/trackA/full/target_metadata.csv

Subcommands: manifest | prep | dock | metadata | child
Parallelism = subprocess.Popen children (no multiprocessing pipes in sandbox).
"""
import argparse, csv, hashlib, json, os, random, re, subprocess, sys, time

ROOT = r"D:\harness\projects\d1_gramneg_vs"
OUTROOT = os.path.join(ROOT, "results", "trackA", "full")
VINA = os.path.join(ROOT, "tools", "vina", "vina.exe")
MEEPREP = os.path.join(ROOT, "tools", "envs", "boltz", "Scripts", "mk_prepare_ligand.exe")
REC_DIR = os.path.join(ROOT, "data", "receptors")
SEED = 42

# target -> (holo receptor pdbqt, box center, panel role, ef1_eligible)
TARGETS = {
    "LpxC_Ec":  ("LpxC_Ec_holo.pdbqt",  [8.621, -32.447, -15.695], "main", True),
    "LpxC_Pa":  ("LpxC_Pa_holo.pdbqt",  [-4.066, 10.272, 12.367], "main", True),
    "DHFR_Ec":  ("DHFR_Ec_holo_3DAU.pdbqt", [27.158, 9.763, -17.934], "main", True),
    "GyrB_Ec":  ("GyrB_Ec_holo.pdbqt",  [-28.103, 7.970, -1.874], "main", True),
    "GyrA_Ec":  ("GyrA_Ec_holo.pdbqt",  [-7.561, 50.828, 29.374], "supplementary", False),
    "FabI_Ec":  ("FabI_Ec_holo.pdbqt",  [12.555, 15.073, 11.276], "main", False),
    "TEM1_Ec":  ("TEM1_Ec_holo.pdbqt",  [3.088, 13.906, 18.215], "main", False),
    "LepB_Ec":  ("LepB_Ec_holo.pdbqt",  [64.921, 22.645, 10.289], "supplementary", False),
    "Gyrase_Ec": ("GyrB_Ec_holo.pdbqt", [-28.103, 7.970, -1.874], "main", True),  # GyrB main arm
}
BOX_SIZES = {"LpxC_Ec": 26, "GyrB_Ec": 26, "Gyrase_Ec": 26, "LepB_Ec": 26}  # amendments #8/#10

ACTIVE_FILES = {t: f"{t}_actives.csv" for t in TARGETS}
DECOY_FILES = {t: f"{t}_decoys_s1.csv" for t in TARGETS}
DECOY_FILES["Gyrase_Ec"] = "Gyrase_Ec_decoys_s1.csv"

def now():
    return __import__("datetime").datetime.now().isoformat(timespec="seconds")

def tdir(target):
    d = os.path.join(OUTROOT, target)
    for sub in ("ligands", "poses", "logs"):
        os.makedirs(os.path.join(d, sub), exist_ok=True)
    return d

def sha256_text(s):
    return hashlib.sha256(s.encode()).hexdigest()

def read_actives(target):
    p = os.path.join(ROOT, "data", "trackA", ACTIVE_FILES[target])
    out = []
    for r in csv.DictReader(open(p, encoding="utf-8-sig")):
        smi = (r.get("canonical_smiles") or "").strip()
        if not smi:
            continue
        cid = (r.get("compound_chembl_id") or r.get("inchikey") or f"a{len(out)}").strip()
        out.append({"ligand_id": "act_" + cid, "label": 1, "smiles": smi,
                    "inchikey": (r.get("inchikey") or "").strip()})
    return out

def read_decoys_unique(target):
    """Deduplicate s1 decoys by decoy_inchikey; keep matched_active_inchikey."""
    p = os.path.join(ROOT, "data", "decoys", DECOY_FILES[target])
    seen, out = set(), []
    for r in csv.DictReader(open(p, encoding="utf-8-sig")):
        smi = (r.get("decoy_canonical_smiles") or "").strip()
        ik = (r.get("decoy_inchikey") or "").strip()
        if not smi:
            continue
        if ik and ik in seen:
            continue
        if ik:
            seen.add(ik)
        out.append({"ligand_id": "dec_" + (ik or r.get("decoy_chembl_id", f"d{len(out)}")),
                    "label": 0, "smiles": smi, "inchikey": ik,
                    "matched_active_inchikey": (r.get("matched_active_inchikey") or "").strip()})
    return out

def build_manifest(target, ratio=20):
    d = tdir(target)
    acts = read_actives(target)
    decs = read_decoys_unique(target)
    rng = random.Random(SEED)
    # group decoys by matched active inchikey
    groups = {}
    for x in decs:
        groups.setdefault(x["matched_active_inchikey"], []).append(x)
    for k in groups:
        groups[k].sort(key=lambda z: z["ligand_id"])
    chosen, used = [], set()
    for a in acts:
        pool = [z for z in groups.get(a["inchikey"], []) if z["ligand_id"] not in used]
        k = min(ratio, len(pool))
        if k > 0:
            pick = rng.sample(pool, k)
            for z in pick:
                used.add(z["ligand_id"]); chosen.append(z)
    rows = acts + chosen
    # amendment #10: mark A-line box-clipped actives (skip docking but keep rows)
    clipped_path = os.path.join(ROOT, "data", "trackA", "box_clipped_ligands.csv")
    clipped_ik = set()
    if os.path.exists(clipped_path):
        for r in csv.DictReader(open(clipped_path, encoding="utf-8-sig")):
            if r.get("target", "").strip() == target and r.get("inchikey"):
                clipped_ik.add(r["inchikey"].strip())
    for r in rows:
        r["box_clipped"] = "True" if (r.get("inchikey") in clipped_ik) else "False"
    for r in rows:
        r["pdbqt"] = os.path.join(d, "ligands", r["ligand_id"] + ".pdbqt")
        r["status"] = "pending"
    subsample = ",".join(sorted(z["inchikey"] for z in chosen if z["inchikey"]))
    subsample_hash = sha256_text(subsample)
    man = os.path.join(d, "manifest.csv")
    with open(man, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "smiles", "pdbqt", "status", "box_clipped"], extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    meta = {"target": target, "n_act": len(acts), "n_dec_pool_unique": len(decs),
            "n_dec_sampled": len(chosen), "decoy_ratio": ratio, "sampling_seed": SEED,
            "decoy_subsample_sha256": subsample_hash, "n_rows": len(rows),
            "manifest": os.path.relpath(man, ROOT), "generated_at": now()}
    json.dump(meta, open(os.path.join(d, "sampling_meta.json"), "w"), indent=1)
    print(json.dumps(meta, indent=1), flush=True)
    return man, meta

def prep_child(manifest, chunk, nchunks):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    rows = list(csv.DictReader(open(manifest, encoding="utf-8-sig")))
    mine = [r for i, r in enumerate(rows) if i % nchunks == chunk]
    ok = fail = skip = 0
    for r in mine:
        if os.path.exists(r["pdbqt"]) and os.path.getsize(r["pdbqt"]) > 100:
            skip += 1; continue
        try:
            m = Chem.MolFromSmiles(r["smiles"])
            if m is None:
                fail += 1; continue
            m = Chem.AddHs(m)
            p = AllChem.ETKDGv3(); p.randomSeed = SEED
            if AllChem.EmbedMolecule(m, p) != 0:
                AllChem.EmbedMolecule(m, useRandomCoords=True, randomSeed=SEED)
            try:
                AllChem.MMFFOptimizeMolecule(m, maxIters=200)
            except Exception:
                pass
            sdf = r["pdbqt"].replace(".pdbqt", ".sdf")
            Chem.MolToMolFile(m, sdf)
            res = subprocess.run([MEEPREP, "-i", sdf, "-o", r["pdbqt"]], capture_output=True,
                                 text=True, encoding="utf-8", errors="replace", timeout=300)
            if res.returncode == 0 and os.path.exists(r["pdbqt"]):
                ok += 1
            else:
                fail += 1
        except Exception:
            fail += 1
    print(json.dumps({"chunk": chunk, "ok": ok, "skip": skip, "fail": fail}), flush=True)

def vina_one(rec, lig, box, cpu, exh, out_pdbqt, logf, size=22):
    if os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 500 and os.path.exists(logf):
        txt = open(logf, encoding="utf-8", errors="replace").read()
        for ln in txt.splitlines():
            m = re.match(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$", ln)
            if m:
                return float(m.group(2)), "ok(reuse)"
    cmd = [VINA, "--receptor", rec, "--ligand", lig, "--center_x", str(box[0]), "--center_y", str(box[1]),
           "--center_z", str(box[2]), "--size_x", str(size), "--size_y", str(size), "--size_z", str(size),
           "--exhaustiveness", str(exh), "--seed", str(SEED), "--num_modes", "3", "--cpu", str(cpu),
           "--out", out_pdbqt]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=3600)
        txt = (res.stdout or "") + (res.stderr or "")
    except subprocess.TimeoutExpired:
        return None, "timeout"
    open(logf, "w", encoding="utf-8", errors="replace").write(txt)
    for ln in txt.splitlines():
        m = re.match(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$", ln)
        if m:
            return float(m.group(2)), "ok"
    return None, "fail"

def dock_child(target, manifest, chunk, nchunks, cpu, exh):
    d = tdir(target)
    rec = os.path.join(REC_DIR, TARGETS[target][0])
    box = TARGETS[target][1]
    rows = list(csv.DictReader(open(manifest, encoding="utf-8-sig")))
    mine = [r for i, r in enumerate(rows) if i % nchunks == chunk]
    outp = os.path.join(d, f"scores_part{chunk}.csv")
    done = 0
    with open(outp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f); w.writerow(["ligand_id", "label", "score", "status", "seed", "exhaustiveness"])
        partial = os.path.join(d, f"scores.partial{chunk}.csv")
        if os.path.exists(partial):
            os.remove(partial)
        pf = open(partial, "a", newline="", encoding="utf-8-sig"); pw = csv.writer(pf); nflush = 0
        for r in mine:
            if str(r.get("box_clipped", "")).lower() == "true":
                row = [r["ligand_id"], r["label"], "", "box_clipped_not_docked", SEED, exh]
                w.writerow(row); f.flush()
                pw.writerow(row); pf.flush(); nflush += 1
                continue
            if not (os.path.exists(r["pdbqt"]) and os.path.getsize(r["pdbqt"]) > 100):
                w.writerow([r["ligand_id"], r["label"], "", "prep_missing", SEED, exh]); f.flush(); continue
            score, status = vina_one(rec, r["pdbqt"], box, cpu, exh,
                                     os.path.join(d, "poses", r["ligand_id"] + "_out.pdbqt"),
                                     os.path.join(d, "logs", "vina_" + r["ligand_id"] + ".log"), size=BOX_SIZES.get(target, 22))
            row = [r["ligand_id"], r["label"], "" if score is None else score, status, SEED, exh]
            w.writerow(row); f.flush()
            pw.writerow(row); pf.flush(); nflush += 1
            done += 1
    try:
        pf.close()
    except Exception:
        pass
    print(json.dumps({"chunk": chunk, "docked": done, "flushed": nflush}), flush=True)

def launch(args_list, logs):
    procs = []
    for args, log in zip(args_list, logs):
        os.makedirs(os.path.dirname(log), exist_ok=True)
        lf = open(log, "w", encoding="utf-8", errors="replace")
        procs.append((subprocess.Popen([sys.executable, os.path.abspath(__file__)] + args,
                                       stdout=lf, stderr=subprocess.STDOUT, cwd=ROOT), lf))
    for p, lf in procs:
        p.wait(); lf.close()

def cmd_prep(target, manifest, nproc):
    d = tdir(target)
    args = [["prep-child", "--manifest", manifest, "--chunk", str(c), "--nchunks", str(nproc)] for c in range(nproc)]
    logs = [os.path.join(d, "logs", f"prep_{c}.log") for c in range(nproc)]
    t0 = time.time(); launch(args, logs)
    rows = list(csv.DictReader(open(manifest, encoding="utf-8-sig")))
    n_ok = sum(1 for r in rows if os.path.exists(r["pdbqt"]) and os.path.getsize(r["pdbqt"]) > 100)
    print(json.dumps({"prep_ok": n_ok, "n": len(rows), "elapsed_s": round(time.time()-t0, 1)}), flush=True)

def cmd_dock(target, manifest, nproc, cpu, exh):
    d = tdir(target)
    args = [["dock-child", "--target", target, "--manifest", manifest, "--chunk", str(c),
             "--nchunks", str(nproc), "--cpu", str(cpu), "--exhaustiveness", str(exh)] for c in range(nproc)]
    logs = [os.path.join(d, "logs", f"dock_{c}.log") for c in range(nproc)]
    t0 = time.time(); launch(args, logs)
    rows = []
    for c in range(nproc):
        p = os.path.join(d, f"scores_part{c}.csv")
        if os.path.exists(p):
            rows += list(csv.DictReader(open(p, encoding="utf-8-sig")))
    out = os.path.join(d, f"{target}_scores.csv")
    with open(out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["ligand_id", "label", "score", "status", "seed", "exhaustiveness"])
        w.writeheader(); w.writerows(rows)
    n_ok = sum(1 for r in rows if r["status"].startswith("ok"))
    n_act = sum(1 for r in rows if str(r["label"]) == "1"); n_dec = len(rows) - n_act
    el = round(time.time() - t0, 1)
    summary = {"target": target, "date": "2026-11-17", "n_rows": len(rows), "n_act": n_act, "n_dec": n_dec,
               "n_ok": n_ok, "n_fail": len(rows) - n_ok,
               "fail_status": {s: sum(1 for r in rows if r["status"] == s) for s in set(r["status"] for r in rows) if not s.startswith("ok")},
               "receptor": TARGETS[target][0], "box": TARGETS[target][1], "box_size": BOX_SIZES.get(target, 22),
               "exhaustiveness": exh, "seed": SEED, "num_modes": 3, "workers": nproc, "cpu_per_worker": cpu,
               "elapsed_s": el, "ligands_per_min": round(len(rows) / max(el/60, 1e-9), 2),
               "panel_role": TARGETS[target][2],
               "source": "B-line _b_trackA_v2.py (amendment #7)"}
    json.dump(summary, open(os.path.join(d, "summary.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1), flush=True)

def cmd_metadata():
    rows = []
    for t in TARGETS:
        d = os.path.join(OUTROOT, t)
        sm = os.path.join(d, "sampling_meta.json"); su = os.path.join(d, "summary.json")
        meta = json.load(open(sm, encoding="utf-8")) if os.path.exists(sm) else {}
        summ = json.load(open(su, encoding="utf-8")) if os.path.exists(su) else {}
        rows.append({"target": t,
                     "decoy_ratio": meta.get("decoy_ratio", ""),
                     "sampling_seed": meta.get("sampling_seed", SEED),
                     "decoy_subsample_sha256": meta.get("decoy_subsample_sha256", ""),
                     "exhaustiveness": summ.get("exhaustiveness", ""),
                     "panel_role": TARGETS[t][2],
                     "box_size": BOX_SIZES.get(t, 22),
                     "n_act": meta.get("n_act", ""), "n_dec_sampled": meta.get("n_dec_sampled", ""),
                     "n_rows": meta.get("n_rows", ""), "status": "done" if summ else "pending"})
    p = os.path.join(OUTROOT, "target_metadata.csv")
    with open(p, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(json.dumps(rows, indent=1), flush=True)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("manifest"); p1.add_argument("--target", required=True); p1.add_argument("--ratio", type=int, default=20)
    p2 = sub.add_parser("prep"); p2.add_argument("--target", required=True); p2.add_argument("--manifest", required=True); p2.add_argument("--nproc", type=int, default=12)
    p3 = sub.add_parser("dock"); p3.add_argument("--target", required=True); p3.add_argument("--manifest", required=True)
    p3.add_argument("--nproc", type=int, default=16); p3.add_argument("--cpu", type=int, default=1)
    p3.add_argument("--exhaustiveness", type=int, default=4)
    sub.add_parser("metadata")
    p4 = sub.add_parser("prep-child"); p4.add_argument("--manifest", required=True); p4.add_argument("--chunk", type=int, required=True); p4.add_argument("--nchunks", type=int, required=True)
    p5 = sub.add_parser("dock-child"); p5.add_argument("--target", required=True); p5.add_argument("--manifest", required=True); p5.add_argument("--chunk", type=int, required=True); p5.add_argument("--nchunks", type=int, required=True); p5.add_argument("--cpu", type=int, default=1); p5.add_argument("--exhaustiveness", type=int, default=4)
    a = ap.parse_args()
    if a.cmd == "manifest":
        build_manifest(a.target, a.ratio)
    elif a.cmd == "prep":
        cmd_prep(a.target, a.manifest, a.nproc)
    elif a.cmd == "dock":
        cmd_dock(a.target, a.manifest, a.nproc, a.cpu, a.exhaustiveness)
    elif a.cmd == "metadata":
        cmd_metadata()
    elif a.cmd == "prep-child":
        prep_child(a.manifest, a.chunk, a.nchunks)
    elif a.cmd == "dock-child":
        dock_child(a.target, a.manifest, a.chunk, a.nchunks, a.cpu, a.exhaustiveness)
