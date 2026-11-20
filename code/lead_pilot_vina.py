# -*- coding: utf-8 -*-
"""Lead-owned G1 LpxC pilot: Vina 1.2.5 docking + enrichment metrics.

Self-contained fallback because the B-line runner had a score-parsing bug.
Box is the A-line authoritative one (4ISA GVR active-site copy nearest catalytic Zn).
Writes to results/trackA/pilot_lpxc_lead/ (does not touch B-line outputs).
"""
import csv, json, math, os, subprocess, sys, time
import concurrent.futures as cf

ROOT = r"D:\harness\projects\d1_gramneg_vs"
VINA = os.path.join(ROOT, r"tools\vina\vina.exe")
REC = os.path.join(ROOT, r"data\receptors\LpxC_Ec_holo.pdbqt")
MAN = os.path.join(ROOT, r"results\trackA\pilot_lpxc\pilot_manifest.csv")
OUT = os.path.join(ROOT, r"results\trackA\pilot_lpxc_lead")
POSE = os.path.join(OUT, "poses")
os.makedirs(POSE, exist_ok=True)

BOX = dict(cx=8.621, cy=-32.447, cz=-15.695, sx=22.0, sy=22.0, sz=22.0)
EXHAUST = 8
SEED = 42
CPU = 8           # two workers x 8 = 16 threads
WORKERS = 2
NUM_MODES = 3


def dock(row):
    lid = row["ligand_id"]
    lig = row["pdbqt"]
    label = int(row["label"])
    outp = os.path.join(POSE, lid.replace("/", "_") + "_out.pdbqt")
    cmd = [VINA, "--receptor", REC, "--ligand", lig,
           "--center_x", str(BOX["cx"]), "--center_y", str(BOX["cy"]), "--center_z", str(BOX["cz"]),
           "--size_x", str(BOX["sx"]), "--size_y", str(BOX["sy"]), "--size_z", str(BOX["sz"]),
           "--exhaustiveness", str(EXHAUST), "--seed", str(SEED), "--cpu", str(CPU),
           "--num_modes", str(NUM_MODES), "--out", outp]
    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=600)
        rc = r.returncode
        err = (r.stderr or "")[-300:]
    except subprocess.TimeoutExpired:
        return dict(ligand_id=lid, label=label, best_affinity="", n_modes=0,
                    status="timeout", elapsed=round(time.time() - t0, 1), rc=-1, err="timeout")
    affs = []
    if os.path.exists(outp):
        with open(outp, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    try:
                        affs.append(float(parts[3]))
                    except (ValueError, IndexError):
                        pass
    best = min(affs) if affs else None
    if best is not None:
        status = "ok"
    elif rc != 0:
        status = "fail"
    else:
        status = "nodata"
    return dict(ligand_id=lid, label=label,
                best_affinity="" if best is None else best,
                n_modes=len(affs), status=status,
                elapsed=round(time.time() - t0, 1), rc=rc, err=err)


def ef_at(ranked_labels, frac):
    n = len(ranked_labels)
    k = max(1, int(round(frac * n)))
    n_act = sum(ranked_labels)
    if n_act == 0:
        return float("nan")
    hits = sum(ranked_labels[:k])
    return (hits / n_act) / frac


def bedroc(ranked_labels, alpha=20.0):
    n = len(ranked_labels)
    n_act = sum(ranked_labels)
    if n_act == 0:
        return float("nan")
    ra = n_act / n
    s = sum(math.exp(-alpha * i / n) for i, y in enumerate(ranked_labels) if y)
    denom = (1 - math.exp(-alpha)) / (math.exp(alpha / n) - 1 + 1e-12) + ra * alpha / (1 - math.exp(-alpha))
    return s / (n_act * denom) if denom else float("nan")


def auc_score(labels, scores):
    # higher is better: use -affinity
    pairs = sorted(zip(scores, labels), reverse=True)
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    rank_sum = 0.0
    for i, (s, y) in enumerate(pairs, 1):
        if y:
            rank_sum += i
    return (rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def main():
    rows = list(csv.DictReader(open(MAN, encoding="utf-8-sig")))
    print(f"[pilot] {len(rows)} ligands, receptor={REC}", flush=True)
    t0 = time.time()
    results = []
    with cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for i, res in enumerate(ex.map(dock, rows), 1):
            results.append(res)
            if i % 10 == 0 or i == len(rows):
                nok = sum(1 for r in results if r["status"] == "ok")
                print(f"[pilot] {i}/{len(rows)} ok={nok} last={res['ligand_id']} {res['status']} "
                      f"best={res['best_affinity']}", flush=True)
    # write scores
    sc_path = os.path.join(OUT, "pilot_scores.csv")
    with open(sc_path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["ligand_id", "label", "vina_score", "n_modes", "status", "elapsed_s", "returncode"])
        for r in results:
            w.writerow([r["ligand_id"], r["label"], r["best_affinity"], r["n_modes"],
                        r["status"], r["elapsed"], r["rc"]])
    ok = [r for r in results if r["status"] == "ok"]
    metrics = {"n_total": len(rows), "n_ok": len(ok),
               "n_failed": len(rows) - len(ok),
               "n_active": sum(r["label"] for r in ok),
               "n_decoy": sum(1 - r["label"] for r in ok),
               "box": BOX, "exhaustiveness": EXHAUST, "seed": SEED,
               "vina": "1.2.5", "box_rule": "A-line 4ISA GVR active-site copy nearest catalytic Zn",
               "elapsed_total_s": round(time.time() - t0, 1)}
    if ok:
        ranked = sorted(ok, key=lambda r: r["best_affinity"])  # most negative first
        labels = [r["label"] for r in ranked]
        neg = [-r["best_affinity"] for r in ranked]
        metrics.update({
            "ef1": ef_at(labels, 0.01), "ef5": ef_at(labels, 0.05), "ef10": ef_at(labels, 0.10),
            "auc": auc_score(labels, neg), "bedroc20": bedroc(labels, 20.0),
            "top10_actives": sum(labels[:10]),
        })
    json.dump(metrics, open(os.path.join(OUT, "pilot_metrics.json"), "w", encoding="utf-8"),
              indent=2, ensure_ascii=False)
    # enrichment curve + ROC
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig_dir = os.path.join(ROOT, "figures")
        os.makedirs(fig_dir, exist_ok=True)
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
        if ok:
            n = len(labels)
            xs = [i / n for i in range(1, n + 1)]
            cum = []
            total = sum(labels)
            c = 0
            for y in labels:
                c += y
                cum.append(c / total if total else 0)
            ideal = [min(1.0, x / (total / n)) for x in xs]
            axes[0].plot(xs, cum, lw=2, label="Vina (pilot)")
            axes[0].plot(xs, ideal, "--", color="grey", lw=1, label="ideal")
            axes[0].plot([0, 1], [0, 1], ":", color="k", lw=0.8, label="random")
            axes[0].set_xlabel("fraction of ranked database")
            axes[0].set_ylabel("fraction of actives recovered")
            axes[0].set_title(f"Enrichment (n_act={total}, n_dec={n - total})")
            axes[0].legend(frameon=False)
            # ROC
            pairs = sorted(zip(neg, labels), reverse=True)
            tp = fp = 0
            tpr = [0.0]
            fpr = [0.0]
            nt = sum(labels)
            nd = n - nt
            for s, y in pairs:
                if y:
                    tp += 1
                else:
                    fp += 1
                tpr.append(tp / nt)
                fpr.append(fp / nd)
            axes[1].plot(fpr, tpr, lw=2, label=f"AUC={metrics.get('auc', float('nan')):.3f}")
            axes[1].plot([0, 1], [0, 1], ":", color="k", lw=0.8)
            axes[1].set_xlabel("false positive rate")
            axes[1].set_ylabel("true positive rate")
            axes[1].set_title("ROC (Vina pilot)")
            axes[1].legend(frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "pilot_lpxc_lead_enrichment.png"), dpi=200)
        fig.savefig(os.path.join(fig_dir, "pilot_lpxc_lead_enrichment.pdf"))
        print("[pilot] figure saved", flush=True)
    except Exception as exc:  # noqa: BLE001
        print(f"[pilot] figure failed: {exc}", flush=True)
    print(json.dumps(metrics, ensure_ascii=False, indent=1), flush=True)


if __name__ == "__main__":
    main()
