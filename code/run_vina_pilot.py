#!/usr/bin/env python
"""Track A / LIT-PCBA docking pilot driver (B-line, task-11 skeleton).

Purpose
-------
Parameterised AutoDock Vina runner + enrichment metrics. Works for:
  - G1 LpxC pilot   : 20 actives + 200 decoys, single receptor, single box
  - Track A full    : per-target receptors (data/receptors/*.pdbqt) + actives/decoys CSV
  - LIT-PCBA anchor : 15 targets x 3 seeds (W5-W6), pure CPU, no GPU lock needed

What it does
------------
1. Reads a receptor PDBQT + a ligand manifest CSV (columns: ligand_id,label[,pdbqt|smiles])
2. Docks every ligand with the official Vina 1.2.5 binary (subprocess, reproducible seed)
3. Writes per-ligand scores and computes EF1%, EF5%, AUC (with bootstrap CI) + hit lists
4. Emits a fully reproducible run manifest (file sha256, vina --version, params, date)

Usage
-----
  python code/trackA/run_vina_pilot.py \
    --receptor data/receptors/lpxc_ec_4is9.pdbqt \
    --ligands data/trackA/lpxc_pilot.csv \
    --box 10.2 -5.3 28.1 22 22 22 \
    --outdir results/trackA/pilot_lpxc \
    --exhaustiveness 8 --seed 42 --cpu 16

Ligand manifest CSV (one row per ligand):
  ligand_id,label,pdbqt            # preferred: prebuilt pdbqt
  act_001,1,/path/act_001.pdbqt
  dec_017,0,/path/dec_017.pdbqt
or with SMILES (requires meeko / mk_prepare_ligand.py on PATH or --mk-prepare):
  ligand_id,label,smiles
  act_001,1,CC(=O)Nc1ccc(O)cc1

Box
---
Either pass --box CX CY CZ SX SY SZ, or --box-json {"center":[x,y,z],"size":[sx,sy,sz]}.
Box size guidance (Vina units = Angstrom): 18-25 A per side for a typical drug-like
pocket; keep <= 30 A (search volume explodes and EF degrades).

Metrics
-------
EF1%  = (actives in top 1% of ranked ligands / total actives) / 0.01
EF5%  = (actives in top 5% / total actives) / 0.05
AUC   = ROC-AUC of -Vina score as a binder score (ties handled by rank)
BEDROC (alpha=20) is also reported (early-recognition metric, useful for LIT-PCBA).
Bootstrap CI: resample ligands with replacement, 1000x, report 2.5-97.5 pct.

This script never writes outside --outdir. It does not need GPU; Vina is CPU-only.
"""
from __future__ import annotations

import argparse, csv, hashlib, json, math, os, random, re, shutil, subprocess, sys, tempfile, time
import datetime

# ---------------------------------------------------------------- utilities
def sha256(path: str, n: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(n)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def now() -> str:
    return datetime.datetime.now().isoformat(timespec="seconds")

def vina_version(vina: str) -> str:
    try:
        out = subprocess.run([vina, "--version"], capture_output=True, text=True, timeout=30)
        return (out.stdout or out.stderr).strip().splitlines()[0]
    except Exception as e:
        return f"UNKNOWN({e})"

# ---------------------------------------------------------------- box parsing
def parse_box(args):
    if args.box:
        if len(args.box) != 6:
            sys.exit("--box needs 6 numbers: CX CY CZ SX SY SZ")
        return {"center": [float(x) for x in args.box[:3]], "size": [float(x) for x in args.box[3:]]}
    if args.box_json:
        j = json.load(open(args.box_json, encoding="utf-8"))
        return {"center": list(map(float, j["center"])), "size": list(map(float, j["size"]))}
    sys.exit("need --box or --box-json")

# ---------------------------------------------------------------- ligand prep
def prepare_ligands(rows, outdir, mk_prepare, args):
    """Return rows with 'pdbqt' resolved; convert SMILES -> pdbqt when needed."""
    ligdir = os.path.join(outdir, "ligands")
    os.makedirs(ligdir, exist_ok=True)
    n_conv = 0
    for r in rows:
        if r.get("pdbqt"):
            p = r["pdbqt"]
            if not os.path.isabs(p):
                p = os.path.join(args.base_dir, p)
            r["pdbqt"] = os.path.abspath(p)
            continue
        smi = (r.get("smiles") or "").strip()
        if not smi:
            sys.exit(f"ligand {r['ligand_id']} has neither pdbqt nor smiles")
        p = os.path.join(ligdir, f"{r['ligand_id']}.pdbqt")
        if not os.path.exists(p):
            if not mk_prepare:
                sys.exit("SMILES given but --mk-prepare not set; install meeko (mk_prepare_ligand.py)")
            smifile = os.path.join(ligdir, f"{r['ligand_id']}.smi")
            open(smifile, "w").write(smi + "\n")
            cmd = [sys.executable, mk_prepare, "-i", smifile, "-o", p]
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            if res.returncode != 0 or not os.path.exists(p):
                raise RuntimeError(f"meeko failed for {r['ligand_id']}: {res.stderr[-500:]}")
            n_conv += 1
        r["pdbqt"] = p
    return rows, n_conv

# ---------------------------------------------------------------- docking
SCORE_RE = re.compile(r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$")

def dock_one(vina, receptor, ligand, box, outdir, exhaustiveness, seed, cpu, timeout=1800):
    """Run Vina on one ligand; return (best_affinity, all_affinities, log)."""
    lig_out = os.path.join(outdir, "poses", os.path.basename(ligand).replace(".pdbqt", "_out.pdbqt"))
    os.makedirs(os.path.dirname(lig_out), exist_ok=True)
    cmd = [vina, "--receptor", receptor, "--ligand", ligand,
           "--center_x", str(box["center"][0]), "--center_y", str(box["center"][1]),
           "--center_z", str(box["center"][2]),
           "--size_x", str(box["size"][0]), "--size_y", str(box["size"][1]),
           "--size_z", str(box["size"][2]),
           "--exhaustiveness", str(exhaustiveness), "--seed", str(seed),
           "--cpu", str(cpu), "--num_modes", "3",
           "--out", lig_out]
    t0 = time.time()
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, [], "TIMEOUT"
    log = (res.stdout or "") + (res.stderr or "")
    scores = [float(m.group(2)) for m in SCORE_RE.finditer(log)]
    return (scores[0] if scores else None), scores, log + f"\n[elapsed {time.time()-t0:.1f}s]"

def run_docking(args, rows, box, outdir):
    rec = os.path.abspath(os.path.join(args.base_dir, args.receptor) if not os.path.isabs(args.receptor) else args.receptor)
    if not os.path.exists(rec):
        sys.exit(f"receptor not found: {rec}")
    vina = args.vina
    ver = vina_version(vina)
    scores_path = os.path.join(outdir, "docking_scores.csv")
    pose_dir = os.path.join(outdir, "poses")
    os.makedirs(pose_dir, exist_ok=True)
    done = 0
    with open(scores_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["ligand_id", "label", "best_affinity", "n_modes", "status", "elapsed_s", "date"])
        for r in rows:
            best, allsc, log = dock_one(vina, rec, r["pdbqt"], box, outdir,
                                        args.exhaustiveness, args.seed, args.cpu)
            outfile = os.path.join(outdir, "poses", os.path.basename(r["pdbqt"]).replace(".pdbqt", "_out.pdbqt"))
            out_exists = os.path.exists(outfile) and os.path.getsize(outfile) > 100
            if best is None and log != "TIMEOUT" and out_exists:
                tail = log.split("-----+")[-1] if "-----+" in log else log
                nums = [float(x) for x in re.findall(r"-?\d+\.\d+", tail)]
                if nums:
                    best = nums[0]
            status = "ok" if best is not None else ("timeout" if log == "TIMEOUT" else "fail")
            r["best_affinity"] = best
            r["status"] = status
            w.writerow([r["ligand_id"], r["label"], "" if best is None else best,
                        len(allsc), status, "", now()])
            try:
                os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)
                open(os.path.join(outdir, "logs", r["ligand_id"] + ".log"), "w", encoding="utf-8", errors="replace").write(log)
            except Exception:
                pass
            done += 1
            if args.verbose or done % 25 == 0:
                print(f"[dock] {done}/{len(rows)} {r['ligand_id']} best={best} ({status})", flush=True)
    manifest = {"date": now(), "vina": ver, "vina_path": vina,
                "receptor": rec, "receptor_sha256": sha256(rec),
                "box": box, "exhaustiveness": args.exhaustiveness, "seed": args.seed,
                "cpu": args.cpu, "n_ligands": len(rows),
                "ligand_manifest": os.path.abspath(args.ligands),
                "ligand_manifest_sha256": sha256(os.path.abspath(args.ligands)) if os.path.exists(args.ligands) else "",
                "source": "B-line run_vina_pilot.py v0.1"}
    json.dump(manifest, open(os.path.join(outdir, "run_manifest.json"), "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)
    return rec, ver, manifest

# ---------------------------------------------------------------- metrics
def ef_at(ranked_labels, frac, n_total):
    k = max(1, int(round(frac * n_total)))
    hits = sum(ranked_labels[:k])
    n_act = sum(ranked_labels)
    if n_act == 0:
        return float("nan")
    return (hits / n_act) / frac

def bedroc(ranked_labels, alpha=20.0):
    n = len(ranked_labels); n_act = sum(ranked_labels)
    if n_act == 0:
        return float("nan")
    ra = n_act / n
    s = 0.0
    for i, y in enumerate(ranked_labels):
        if y:
            s += math.exp(-alpha * i / n)
    denom = (1 - math.exp(-alpha)) / (math.exp(alpha / n) - 1 + 1e-12) + ra * alpha / (1 - math.exp(-alpha))
    return s / (n_act * denom) if denom else float("nan")

def auc_score(labels, scores):
    """ROC-AUC via rank statistic (higher score = binder). scores: -affinity."""
    pairs = sorted(zip(scores, labels), reverse=True)
    n_pos = sum(labels); n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    tp = fp = 0; auc = 0.0
    for _, y in pairs:
        if y:
            tp += 1
        else:
            fp += 1
            auc += tp
    return auc / (n_pos * n_neg)

def bootstrap_ci(rows, fn, n_boot=1000, seed=42):
    rng = random.Random(seed)
    n = len(rows)
    vals = []
    for _ in range(n_boot):
        sample = [rows[rng.randrange(n)] for _ in range(n)]
        v = fn(sample)
        if not math.isnan(v):
            vals.append(v)
    if not vals:
        return (float("nan"), float("nan"))
    vals.sort()
    return (vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals)) - 1])

def compute_metrics(rows):
    ok = [r for r in rows if r.get("status") == "ok" and r.get("best_affinity") is not None]
    # rank: more negative affinity = better -> score = -affinity
    ranked = sorted(ok, key=lambda r: r["best_affinity"])
    labels = [int(r["label"]) for r in ranked]
    scores = [-r["best_affinity"] for r in ranked]
    out = {"n_ligands": len(rows), "n_docked": len(ok), "n_actives": sum(labels),
           "n_decoys": len(labels) - sum(labels),
           "EF1": ef_at(labels, 0.01, len(labels)),
           "EF5": ef_at(labels, 0.05, len(labels)),
           "AUC": auc_score(labels, scores),
           "BEDROC20": bedroc(labels, 20.0)}
    out["EF1_CI95"] = bootstrap_ci(ok, lambda s: ef_at(
        [int(r["label"]) for r in sorted(s, key=lambda r: r["best_affinity"])], 0.01, len(s)))
    out["EF5_CI95"] = bootstrap_ci(ok, lambda s: ef_at(
        [int(r["label"]) for r in sorted(s, key=lambda r: r["best_affinity"])], 0.05, len(s)))
    out["AUC_CI95"] = bootstrap_ci(ok, lambda s: auc_score(
        [int(r["label"]) for r in sorted(s, key=lambda r: r["best_affinity"])],
        [-r["best_affinity"] for r in sorted(s, key=lambda r: r["best_affinity"])]))
    out["top10"] = [{"ligand_id": r["ligand_id"], "label": int(r["label"]),
                     "affinity": r["best_affinity"]} for r in ranked[:10]]
    # rank-normalised score for cross-target comparability
    if ranked:
        for i, r in enumerate(ranked):
            r["rank"] = i + 1
            r["rank_score"] = 1 - i / max(1, len(ranked) - 1)
    return out, ranked

# ---------------------------------------------------------------- main
def read_manifest(path):
    if not os.path.exists(path):
        sys.exit(f"ligand manifest not found: {path}")
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        rd = csv.DictReader(f)
        need = {"ligand_id", "label"}
        if not need.issubset(set(rd.fieldnames or [])):
            sys.exit(f"manifest needs columns {need}; got {rd.fieldnames}")
        for r in rd:
            r["label"] = int(float(r["label"]))
            rows.append(r)
    return rows

def main():
    ap = argparse.ArgumentParser(description="Vina docking pilot + enrichment metrics")
    ap.add_argument("--receptor", required=True, help="receptor PDBQT path")
    ap.add_argument("--ligands", required=True, help="ligand manifest CSV (ligand_id,label[,pdbqt|smiles])")
    ap.add_argument("--box", nargs=6, type=float, help="CX CY CZ SX SY SZ")
    ap.add_argument("--box-json", help='{"center":[x,y,z],"size":[sx,sy,sz]}')
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--vina", default=os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "tools", "vina", "vina.exe"))
    ap.add_argument("--mk-prepare", default=None, help="path to meeko mk_prepare_ligand.py")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--cpu", type=int, default=16)
    ap.add_argument("--base-dir", default=os.getcwd())
    ap.add_argument("--score-only", action="store_true",
                    help="skip docking, read existing docking_scores.csv in outdir and compute metrics")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    box = parse_box(args)
    rows = read_manifest(os.path.abspath(args.ligands))
    if args.score_only:
        sp = os.path.join(outdir, "docking_scores.csv")
        with open(sp, newline="", encoding="utf-8-sig") as f:
            rows = list(csv.DictReader(f))
        for r in rows:
            r["label"] = int(float(r["label"]))
            r["best_affinity"] = float(r["best_affinity"]) if r.get("best_affinity") not in (None, "") else None
    else:
        rows, n_conv = prepare_ligands(rows, outdir, args.mk_prepare, args)
        print(f"[prep] {len(rows)} ligands ({n_conv} converted from SMILES)", flush=True)
        rec, ver, manifest = run_docking(args, rows, box, outdir)
        print(f"[docking] vina={ver} receptor={rec}", flush=True)

    metrics, ranked = compute_metrics(rows)
    with open(os.path.join(outdir, "metrics.json"), "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=1, ensure_ascii=False, default=str)
    with open(os.path.join(outdir, "ranked_ligands.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["rank", "ligand_id", "label", "best_affinity", "rank_score"])
        for r in ranked:
            w.writerow([r.get("rank", ""), r["ligand_id"], r["label"],
                        r["best_affinity"], round(r.get("rank_score", float("nan")), 6)])
    print(json.dumps(metrics, indent=1, ensure_ascii=False, default=str), flush=True)
    print(f"\n[done] results in {outdir}", flush=True)

if __name__ == "__main__":
    main()
