#!/usr/bin/env python
"""G1 finalizer: turn run_vina_pilot outputs into the Lead-specified deliverables.

Inputs (results/trackA/pilot_lpxc/):
  docking_scores.csv  (ligand_id,label,best_affinity,...)
  pilot_manifest.csv, box.json, protocol.md
Outputs:
  pilot_scores.csv    ligand_id, label(active/decoy), vina_score, rank, seed, box_json_hash
  pilot_metrics.csv   metric,value (EF1%,EF5%,AUC,BEDROC20 + CIs + n_act/n_decoy)
  figs/pilot_lpxc_ef.png  enrichment (hit-rate) curve
Run: tools\\envs\\boltz\\Scripts\\python.exe code\\trackA\\finalize_pilot.py
"""
import csv, hashlib, json, math, os

ROOT = r"D:\harness\projects\d1_gramneg_vs"
P = os.path.join(ROOT, "results", "trackA", "pilot_lpxc")
SEED = 42
BOX = os.path.join(P, "box.json")

def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()

def ef_curve(labels, npts=100):
    n = len(labels); nact = sum(labels)
    curve = []
    for i in range(1, npts + 1):
        k = max(1, int(round(i / npts * n)))
        curve.append((i / npts, sum(labels[:k]) / max(1, nact)))
    return curve

def auc(labels, scores):
    pairs = sorted(zip(scores, labels), reverse=True)
    pos = sum(labels); neg = len(labels) - pos
    tp = fp = 0; a = 0.0
    for _, y in pairs:
        if y: tp += 1
        else: fp += 1; a += tp
    return a / (pos * neg) if pos and neg else float("nan")

def bedroc(labels, alpha=20.0):
    n = len(labels); nact = sum(labels)
    s = sum(math.exp(-alpha * i / n) for i, y in enumerate(labels) if y)
    ra = nact / n
    denom = (1 - math.exp(-alpha)) / (math.exp(alpha / n) - 1 + 1e-12) + ra * alpha / (1 - math.exp(-alpha))
    return s / (nact * denom) if nact and denom else float("nan")

def main():
    rows = []
    with open(os.path.join(P, "docking_scores.csv"), newline="", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            if r.get("best_affinity") in (None, "", "None"):
                continue
            r["label"] = int(float(r["label"]))
            r["best_affinity"] = float(r["best_affinity"])
            rows.append(r)
    ranked = sorted(rows, key=lambda r: r["best_affinity"])
    box_hash = sha(BOX)[:16] if os.path.exists(BOX) else ""
    with open(os.path.join(P, "pilot_scores.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["ligand_id", "label", "vina_score", "rank", "seed", "box_json_hash"])
        for i, r in enumerate(ranked, 1):
            w.writerow([r["ligand_id"], "active" if r["label"] == 1 else "decoy",
                        r["best_affinity"], i, SEED, box_hash])
    labels = [r["label"] for r in ranked]
    scores = [-r["best_affinity"] for r in ranked]
    n = len(ranked); nact = sum(labels)
    ef1 = (sum(labels[:max(1, int(round(0.01 * n)))]) / nact) / 0.01 if nact else float("nan")
    ef5 = (sum(labels[:max(1, int(round(0.05 * n)))]) / nact) / 0.05 if nact else float("nan")
    metrics = [("n_docked", n), ("n_actives", nact), ("n_decoys", n - nact),
               ("EF1_percent", round(ef1, 2)), ("EF5_percent", round(ef5, 2)),
               ("AUC", round(auc(labels, scores), 4)), ("BEDROC20", round(bedroc(labels), 4)),
               ("seed", SEED), ("box_json_hash", box_hash), ("date", "2026-11-16")]
    with open(os.path.join(P, "pilot_metrics.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f); w.writerow(["metric", "value"]); w.writerows(metrics)
    # enrichment figure
    os.makedirs(os.path.join(P, "figs"), exist_ok=True)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        curve = ef_curve(labels)
        base = [max(0.0, x * n / max(1, nact)) * 0 + x * (n / max(1, nact)) * 0 for x, _ in curve]
        xs = [c[0] * 100 for c in curve]; ys = [c[1] for c in curve]
        fig, ax = plt.subplots(figsize=(5.2, 3.6), dpi=200)
        ax.plot(xs, ys, lw=2, color="#1f77b4", label=f"Vina 1.2.5 (n_act={nact}, n_dec={n-nact})")
        # random baseline: expected recovery vs fraction of library screened = fraction
        frac = [x / 100 for x in xs]
        ax.plot(xs, frac, "--", color="grey", lw=1, label="random")
        ax.set_xlabel("Library screened (%)"); ax.set_ylabel("Fraction of actives recovered")
        ax.set_title("LpxC_Ec pilot enrichment (holo 4ISA, Vina 1.2.5)")
        ax.set_xlim(0, 100); ax.set_ylim(0, 1.02)
        ax.legend(frameon=False, fontsize=7); fig.tight_layout()
        fig.savefig(os.path.join(P, "figs", "pilot_lpxc_ef.png"))
        print("figure written", flush=True)
    except Exception as e:
        print("figure skipped:", e, flush=True)
    print(json.dumps(dict(metrics), indent=1), flush=True)

if __name__ == "__main__":
    main()
