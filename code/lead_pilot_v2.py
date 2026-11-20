# -*- coding: utf-8 -*-
"""Lead-owned pilot v2: Zn-retained vs Zn-stripped LpxC, and DHFR non-metal control.

Runs three docking configurations with the same Vina settings and reports
enrichment metrics + cross-run comparison. Writes to results/trackA/pilot_v2/.
"""
import csv, json, math, os, subprocess, sys, time
import concurrent.futures as cf

ROOT = r"D:\harness\projects\d1_gramneg_vs"
VINA = os.path.join(ROOT, r"tools\vina\vina.exe")
MKPREP = os.path.join(ROOT, r"tools\envs\boltz\Scripts\mk_prepare_ligand.exe")
OUT = os.path.join(ROOT, r"results\trackA\pilot_v2")
LIGDIR = os.path.join(OUT, "dhfr_ligands")
os.makedirs(OUT, exist_ok=True)
os.makedirs(LIGDIR, exist_ok=True)

BOX = dict(cx=8.621, cy=-32.447, cz=-15.695, sx=22.0, sy=22.0, sz=22.0)
BOX_DHFR = dict(cx=27.158, cy=9.763, cz=-17.934, sx=22.0, sy=22.0, sz=22.0)
EXHAUST, SEED, CPU, WORKERS, NUM_MODES = 8, 42, 6, 2, 3


def dock_one(vina, receptor, ligand, box, outdir):
    outp = os.path.join(outdir, os.path.basename(ligand).replace(".pdbqt", "_out.pdbqt"))
    cmd = [vina, "--receptor", receptor, "--ligand", ligand,
           "--center_x", str(box["cx"]), "--center_y", str(box["cy"]), "--center_z", str(box["cz"]),
           "--size_x", str(box["sx"]), "--size_y", str(box["sy"]), "--size_z", str(box["sz"]),
           "--exhaustiveness", str(EXHAUST), "--seed", str(SEED), "--cpu", str(CPU),
           "--num_modes", str(NUM_MODES), "--out", outp]
    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=600)
        rc = r.returncode
    except subprocess.TimeoutExpired:
        return None, "timeout", round(time.time() - t0, 1)
    affs = []
    if os.path.exists(outp):
        for line in open(outp, encoding="utf-8", errors="replace"):
            if line.startswith("REMARK VINA RESULT:"):
                try:
                    affs.append(float(line.split()[3]))
                except (ValueError, IndexError):
                    pass
    best = min(affs) if affs else None
    status = "ok" if best is not None else ("fail" if rc != 0 else "nodata")
    return best, status, round(time.time() - t0, 1)


def run_lpxc(tag):
    man = os.path.join(ROOT, r"results\trackA\pilot_lpxc\pilot_manifest.csv")
    rows = list(csv.DictReader(open(man, encoding="utf-8-sig")))
    rec = os.path.join(ROOT, r"data\receptors", "LpxC_Ec_holo.pdbqt" if tag == "lpxc_zn" else "LpxC_Ec_holo_noZn.pdbqt")
    return rows, rec, BOX


def prep_dhfr():
    man = os.path.join(ROOT, r"data\receptors\DHFR_pilot_manifest.csv")
    rows = []
    for r in csv.DictReader(open(man, encoding="utf-8-sig")):
        lid = r["compound_id"] + "_" + r["label"]
        p = os.path.join(LIGDIR, lid + ".pdbqt")
        if not os.path.exists(p):
            smi = os.path.join(LIGDIR, lid + ".smi")
            open(smi, "w", encoding="utf-8").write(r["smiles"] + "\n")
            res = subprocess.run([MKPREP, "-i", smi, "-o", p], capture_output=True, text=True)
            if res.returncode != 0 or not os.path.exists(p):
                continue
        rows.append(dict(ligand_id=lid, label="1" if r["label"] == "active" else "0", pdbqt=p))
    return rows, os.path.join(ROOT, r"data\receptors\DHFR_Ec_holo_3DAU.pdbqt"), BOX_DHFR


def metrics(rows):
    ok = [r for r in rows if r["status"] == "ok" and r["vina_score"] is not None]
    labels = [1 if str(r["label"]) == "1" else 0 for r in ok]
    scores = [r["vina_score"] for r in ok]
    n, na = len(labels), sum(labels)
    rk = sorted(zip(scores, labels))
    lab = [y for _, y in rk]
    def ef(f):
        k = max(1, int(round(f * n)))
        return (sum(lab[:k]) / na) / f if na else float("nan")
    pos = [i for i, y in enumerate(lab) if y]
    neg = [i for i, y in enumerate(lab) if not y]
    auc = (sum(1 for p in pos for q in neg if p < q) / (len(pos) * len(neg))) if pos and neg else float("nan")
    rie = sum(math.exp(-20.0 * (i + 1) / n) for i in pos)
    bedroc = rie * (1 - math.exp(-20.0)) / (na * (1 - math.exp(-20.0 * na / n))) if na else float("nan")
    return dict(n=n, n_active=na, ef1=ef(0.01), ef5=ef(0.05), ef10=ef(0.10), auc=auc, bedroc20=bedroc)


def main():
    runs = {}
    for tag, (rows, rec, box) in {
        "lpxc_zn": run_lpxc("lpxc_zn"),
        "lpxc_nozn": run_lpxc("lpxc_nozn"),
    }.items():
        print(f"[{tag}] receptor={rec} n={len(rows)}", flush=True)
        outdir = os.path.join(OUT, tag, "poses")
        os.makedirs(outdir, exist_ok=True)
        res = []
        t0 = time.time()
        with cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
            futs = {ex.submit(dock_one, VINA, rec, r["pdbqt"], box, outdir): r for r in rows}
            for i, f in enumerate(cf.as_completed(futs), 1):
                r = futs[f]
                best, status, el = f.result()
                res.append(dict(ligand_id=r["ligand_id"], label=r["label"], vina_score=best, status=status, elapsed=el))
                if i % 25 == 0:
                    print(f"[{tag}] {i}/{len(rows)} ok={sum(1 for x in res if x['status']=='ok')}", flush=True)
        with open(os.path.join(OUT, f"{tag}_scores.csv"), "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.writer(fh)
            w.writerow(["ligand_id", "label", "vina_score", "status", "elapsed_s"])
            for r in res:
                w.writerow([r["ligand_id"], r["label"], "" if r["vina_score"] is None else r["vina_score"], r["status"], r["elapsed"]])
        runs[tag] = dict(rows=res, metrics=metrics(res), elapsed=round(time.time() - t0, 1),
                         receptor=os.path.basename(rec), box=box)
        print(f"[{tag}] done {runs[tag]['elapsed']}s metrics={ {k: (round(v,3) if isinstance(v,float) else v) for k,v in runs[tag]['metrics'].items()} }", flush=True)

    # DHFR control
    rows, rec, box = prep_dhfr()
    tag = "dhfr"
    print(f"[dhfr] receptor={rec} n={len(rows)}", flush=True)
    outdir = os.path.join(OUT, tag, "poses")
    os.makedirs(outdir, exist_ok=True)
    res = []
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(dock_one, VINA, rec, r["pdbqt"], box, outdir): r for r in rows}
        for i, f in enumerate(cf.as_completed(futs), 1):
            r = futs[f]
            best, status, el = f.result()
            res.append(dict(ligand_id=r["ligand_id"], label=r["label"], vina_score=best, status=status, elapsed=el))
            if i % 25 == 0:
                print(f"[dhfr] {i}/{len(rows)} ok={sum(1 for x in res if x['status']=='ok')}", flush=True)
    with open(os.path.join(OUT, f"{tag}_scores.csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["ligand_id", "label", "vina_score", "status", "elapsed_s"])
        for r in res:
            w.writerow([r["ligand_id"], r["label"], "" if r["vina_score"] is None else r["vina_score"], r["status"], r["elapsed"]])
    runs[tag] = dict(rows=res, metrics=metrics(res), elapsed=round(time.time() - t0, 1),
                     receptor=os.path.basename(rec), box=box)
    print(f"[dhfr] done {runs[tag]['elapsed']}s metrics={ {k: (round(v,3) if isinstance(v,float) else v) for k,v in runs[tag]['metrics'].items()} }", flush=True)

    # summary table + delta
    summary = {tag: dict(n=v["metrics"]["n"], n_active=v["metrics"]["n_active"],
                         ef1=round(v["metrics"]["ef1"], 4), ef5=round(v["metrics"]["ef5"], 4),
                         ef10=round(v["metrics"]["ef10"], 4), auc=round(v["metrics"]["auc"], 4),
                         bedroc20=round(v["metrics"]["bedroc20"], 4),
                         receptor=v["receptor"], elapsed_s=v["elapsed"]) for tag, v in runs.items()}
    summary["delta_zn_minus_nozn"] = dict(ef1=round(runs["lpxc_zn"]["metrics"]["ef1"] - runs["lpxc_nozn"]["metrics"]["ef1"], 4),
                                          ef5=round(runs["lpxc_zn"]["metrics"]["ef5"] - runs["lpxc_nozn"]["metrics"]["ef5"], 4),
                                          auc=round(runs["lpxc_zn"]["metrics"]["auc"] - runs["lpxc_nozn"]["metrics"]["auc"], 4))
    json.dump(summary, open(os.path.join(OUT, "pilot_v2_metrics.json"), "w", encoding="utf-8"), indent=2, ensure_ascii=False)
    with open(os.path.join(OUT, "pilot_v2_metrics.csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["run", "n", "n_active", "ef1", "ef5", "ef10", "auc", "bedroc20", "receptor", "elapsed_s"])
        for tag in ("lpxc_zn", "lpxc_nozn", "dhfr"):
            s = summary[tag]
            w.writerow([tag, s["n"], s["n_active"], s["ef1"], s["ef5"], s["ef10"], s["auc"], s["bedroc20"], s["receptor"], s["elapsed_s"]])
    # figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(13, 4))
        for ax, tag, title in zip(axes, ("lpxc_zn", "lpxc_nozn", "dhfr"),
                                  ("LpxC_Ec Zn-retained", "LpxC_Ec Zn-stripped", "DHFR (non-metal control)")):
            ok = [r for r in runs[tag]["rows"] if r["status"] == "ok"]
            ok.sort(key=lambda r: r["vina_score"])
            lab = [1 if str(r["label"]) == "1" else 0 for r in ok]
            n = len(lab); tot = sum(lab); c = 0; cum = []
            for y in lab:
                c += y
                cum.append(c / tot if tot else 0)
            xs = [i / n for i in range(1, n + 1)]
            m = runs[tag]["metrics"]
            ax.plot(xs, cum, lw=2)
            ax.plot([0, 1], [0, 1], ":", color="k", lw=0.8)
            ax.set_title(f"{title}\nAUC={m['auc']:.3f}, EF5={m['ef5']:.2f}", fontsize=10)
            ax.set_xlabel("fraction of ranked database")
            ax.set_ylabel("fraction of actives recovered")
        fig.tight_layout()
        fig.savefig(os.path.join(ROOT, r"figures\pilot_v2_zn_nozn_dhfr.png"), dpi=200)
        fig.savefig(os.path.join(ROOT, r"figures\pilot_v2_zn_nozn_dhfr.pdf"))
        print("[figure] saved", flush=True)
    except Exception as exc:  # noqa: BLE001
        print(f"[figure] failed: {exc}", flush=True)
    print(json.dumps(summary, ensure_ascii=False, indent=1), flush=True)


if __name__ == "__main__":
    main()
