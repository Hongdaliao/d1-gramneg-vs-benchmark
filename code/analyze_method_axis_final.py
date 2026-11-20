"""Final six-method-class analysis from the frozen combined_method_axis_scores.csv.

Formal statistics only (not the preview columns in the manifest):
  - main 8-target pooled: Vina, DiffDock, Boltz-2, DrugCLIP
  - Uni-Dock pooled: 7 targets (exclude FabI)
  - common 7-target sensitivity: all methods excluding FabI
  - paired bootstrap vs Vina; target-level pooled bootstrap.
No v3 edits; no GPU.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(r"D:\harness\projects\d1_gramneg_vs")
BASE = ROOT / "results/trackA/method_axis"
OUT = BASE / "analysis"
LONG = BASE / "combined_method_axis_scores.csv"
SUBSET = BASE / "subsets_diffdock.csv"
TARGETS = ["LpxC_Ec", "LpxC_Pa", "DHFR_Ec", "Gyrase_Ec", "GyrA_Ec", "GyrB_Ec", "FabI_Ec", "TEM1_Ec"]
ALL9 = TARGETS + ["LepB_Ec"]
COMMON7 = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec", "GyrA_Ec", "TEM1_Ec"]
MAIN8 = TARGETS
N_BOOT = 2000
SEED = 42


def auc(y, s):
    y = np.asarray(y, int); s = np.asarray(s, float)
    if y.sum() == 0 or y.sum() == len(y):
        return float("nan")
    try:
        return float(roc_auc_score(y, s))
    except ValueError:
        return float("nan")


def boot_ci(y, s, seed):
    y = np.asarray(y, int); s = np.asarray(s, float)
    pos = np.where(y == 1)[0]; neg = np.where(y == 0)[0]
    rng = np.random.default_rng(seed); vals = []
    for _ in range(N_BOOT):
        idx = np.r_[rng.choice(pos, len(pos), replace=True), rng.choice(neg, len(neg), replace=True)]
        v = auc(y[idx], s[idx])
        if np.isfinite(v): vals.append(v)
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))) if vals else (float("nan"), float("nan"))


def top5(y, s):
    y = np.asarray(y, int); s = np.asarray(s, float); k = min(5, len(y))
    return float(y[np.argsort(-s)[:k]].sum() / min(5, max(1, y.sum())))


def rho(x, y):
    x = pd.Series(np.asarray(x, float)); y = pd.Series(np.asarray(y, float)); ok = x.notna() & y.notna()
    return float(x[ok].rank().corr(y[ok].rank())) if ok.sum() >= 3 else float("nan")


def load():
    long = pd.read_csv(LONG)
    sub = pd.read_csv(SUBSET)[["target", "ligand_id", "label", "subset_id", "vina_score"]]
    sub["vina_score"] = pd.to_numeric(sub["vina_score"], errors="coerce")
    sub = sub.rename(columns={"label": "sub_label"})
    long = long.merge(sub, on=["target", "ligand_id"], how="left")
    return long, sub


def analyze_method(long, method, targets):
    d = long[(long["method"] == method) & (long["target"].isin(targets))].copy()
    d = d[d["status"].astype(str).str.lower().eq("ok") & d["score"].notna() & d["vina_score"].notna()].copy()
    # orient so higher score = better
    d["oriented"] = np.where(d["higher_is_better"].astype(int).eq(1), d["score"], -d["score"])
    rows, paired = [], []
    for t in targets:
        s = d[d["target"] == t].copy()
        if s.empty: continue
        y = s["label"].to_numpy(int); sc = s["oriented"].to_numpy(float)
        if t == "LepB_Ec":
            continue
        lo, hi = boot_ci(y, sc, SEED + ALL9.index(t) + method_offset(method))
        rows.append({"target": t, "method": method, "n_ligands": int((long[(long["method"] == method) & (long["target"] == t)]).shape[0]),
                     "n_analysis": len(s), "n_active": int(y.sum()), "n_decoy": int((y == 0).sum()),
                     "auc": auc(y, sc), "auc_ci_lo": lo, "auc_ci_hi": hi,
                     "top5_recall": top5(y, sc), "rho_vina": rho(s["score"], s["vina_score"]),
                     "status": "ok"})
        # paired bootstrap vs vina on same subset pairs
        pairs = []
        for _, g in s.groupby("subset_id"):
            if len(g) == 2 and set(g["label"].astype(int)) == {0, 1}:
                pairs.append(g)
        if pairs:
            rng = np.random.default_rng(SEED + 600 + ALL9.index(t) + method_offset(method)); deltas = []
            for _ in range(N_BOOT):
                sample = pd.concat([pairs[i] for i in rng.integers(0, len(pairs), len(pairs))], ignore_index=True)
                yy = sample["label"].to_numpy(int)
                a = auc(yy, sample["oriented"].to_numpy(float))
                b = auc(yy, -sample["vina_score"].to_numpy(float))
                if np.isfinite(a) and np.isfinite(b): deltas.append(a - b)
            paired.append({"target": t, "method": method, "n_pairs": len(pairs),
                           "delta_vs_vina": auc(y, sc) - auc(y, -s["vina_score"].to_numpy(float)),
                           "delta_ci_lo": float(np.percentile(deltas, 2.5)) if deltas else float("nan"),
                           "delta_ci_hi": float(np.percentile(deltas, 97.5)) if deltas else float("nan")})
    return pd.DataFrame(rows), pd.DataFrame(paired)


def method_offset(m):
    return {"vina": 0, "unidock": 100, "diffdock_main": 200, "boltz2": 300, "drugclip": 400, "diffdock_singlechain": 500}.get(m, 0)


def pooled(frame, targets):
    vals = frame[frame["target"].isin(targets)].set_index("target")["auc"].dropna()
    if vals.empty:
        return {"n_targets": 0, "mean_auc": None, "ci95": [None, None], "targets": {}}
    arr = vals.to_numpy(float); rng = np.random.default_rng(SEED)
    draws = [float(np.mean(arr[rng.integers(0, len(arr), len(arr))])) for _ in range(N_BOOT)]
    return {"n_targets": len(vals), "mean_auc": float(np.mean(arr)),
            "ci95": [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))], "targets": vals.to_dict()}


def coverage(long):
    rows = []
    for (m, t), g in long.groupby(["method", "target"]):
        if t not in ALL9: continue
        ok = int(g["status"].astype(str).str.lower().eq("ok").sum())
        st = g["status"].fillna("NA").astype(str).value_counts().to_dict()
        rows.append({"method": m, "target": t, "attempted": len(g), "ok": ok, "na": len(g) - ok,
                     "status_counts": json.dumps(st)})
    return pd.DataFrame(rows)


def make_figure(metrics, pooled_main, paired, cov):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4))
    methods = ["vina", "unidock", "diffdock_main", "boltz2", "drugclip"]
    labels = {"vina": "Vina", "unidock": "Uni-Dock", "diffdock_main": "DiffDock", "boltz2": "Boltz-2", "drugclip": "DrugCLIP"}
    colors = {"vina": "#225ea8", "unidock": "#41b6c4", "diffdock_main": "#cb181d", "boltz2": "#6a51a3", "drugclip": "#e6550d"}
    # A per-target AUC
    ax = axes[0, 0]
    x = np.arange(len(TARGETS)); w = 0.16
    for i, m in enumerate(methods):
        d = metrics[metrics["method"] == m].set_index("target").reindex(TARGETS)
        ax.bar(x + (i - 2) * w, d["auc"], width=w, color=colors[m], label=labels[m])
    ax.set_xticks(x); ax.set_xticklabels(TARGETS, rotation=40, ha="right", fontsize=6)
    ax.set_ylim(0, 1); ax.set_ylabel("AUC"); ax.set_title("A. Per-target AUC (5 methods)")
    ax.legend(fontsize=5.5, frameon=False, ncol=2)
    # B pooled
    ax = axes[0, 1]
    for i, m in enumerate(["vina", "diffdock_main", "boltz2", "drugclip", "unidock"]):
        p = pooled_main.get(m, {})
        mean = p.get("mean_auc", np.nan); lo, hi = p.get("ci95", [np.nan, np.nan])
        ax.bar(i, mean, color=colors[m], width=0.6)
        ax.errorbar(i, mean, yerr=[[mean - lo], [hi - mean]], fmt="none", color="0.3", capsize=3, lw=0.8)
    ax.set_xticks(range(5)); ax.set_xticklabels(["Vina", "DiffDock", "Boltz-2", "DrugCLIP", "Uni-Dock"], rotation=25, ha="right", fontsize=6)
    ax.set_ylim(0, 1); ax.set_ylabel("Pooled AUC"); ax.set_title("B. Pooled AUC (8; Uni-Dock 7)")
    # C paired delta vs Vina
    ax = axes[0, 2]
    for i, m in enumerate(["diffdock_main", "boltz2", "drugclip", "unidock"]):
        d = paired[paired["method"] == m].set_index("target").reindex(TARGETS)
        ax.scatter(d["delta_vs_vina"], np.arange(len(TARGETS)) + (i - 1.5) * 0.08, s=8, color=colors[m], label=labels[m])
    ax.axvline(0, color="0.5", ls="--", lw=0.7); ax.set_yticks(range(len(TARGETS))); ax.set_yticklabels(TARGETS, fontsize=6)
    ax.set_xlabel("ΔAUC vs Vina"); ax.set_title("C. Paired per-target ΔAUC")
    # D singlechain
    ax = axes[1, 0]
    sc = metrics[metrics["method"] == "diffdock_singlechain"].set_index("target").reindex(["GyrA_Ec", "FabI_Ec"])
    ax.bar([0, 1], sc["auc"].values, color="#cb181d")
    ax.set_xticks([0, 1], ["GyrA_Ec", "FabI_Ec"]); ax.set_ylim(0, 1); ax.set_ylabel("AUC"); ax.set_title("D. DiffDock single-chain (descriptive)")
    # E common-7 sensitivity
    ax = axes[1, 1]
    for i, m in enumerate(["vina", "diffdock_main", "boltz2", "drugclip", "unidock"]):
        # placeholder: use main pooled targets excluding FabI from metrics
        d = metrics[(metrics["method"] == m) & metrics["target"].isin(COMMON7)].set_index("target")["auc"].dropna()
        if len(d):
            ax.bar(i, d.mean(), color=colors[m], width=0.6)
    ax.set_xticks(range(5)); ax.set_xticklabels(["Vina", "DiffDock", "Boltz-2", "DrugCLIP", "Uni-Dock"], rotation=25, ha="right", fontsize=6)
    ax.set_ylim(0, 1); ax.set_ylabel("Mean AUC"); ax.set_title("E. Common 7-target sensitivity")
    # F coverage
    ax = axes[1, 2]
    for i, m in enumerate(methods):
        d = cov[cov["method"] == m]
        ax.bar(i, d["ok"].sum(), color=colors[m], width=0.6)
        ax.set_xticks(range(5)); ax.set_xticklabels(["Vina", "Uni-Dock", "DiffDock", "Boltz-2", "DrugCLIP"], rotation=25, ha="right", fontsize=6)
    ax.set_ylabel("ok rows"); ax.set_title("F. Coverage (ok rows)")
    fig.tight_layout(); fig.savefig(ROOT / "figures/v3/fig7_method_axis.png", dpi=300, bbox_inches="tight")
    fig.savefig(ROOT / "figures/v3/fig7_method_axis.pdf", bbox_inches="tight"); plt.close(fig)


def main():
    long, sub = load()
    all_metrics, all_paired = [], []
    for m in ["vina", "diffdock_main", "boltz2", "drugclip", "unidock"]:
        mf, pf = analyze_method(long, m, TARGETS)
        all_metrics.append(mf); all_paired.append(pf)
    sc = long[(long["method"] == "diffdock_singlechain") & long["status"].astype(str).str.lower().eq("ok")].copy()
    sc = sc.merge(sub[["target", "ligand_id", "subset_id", "vina_score"]], on=["target", "ligand_id"], how="left")
    sc["oriented"] = np.where(sc["higher_is_better"].astype(int).eq(1), sc["score"], -sc["score"])
    sc_rows = []
    for t in ["GyrA_Ec", "FabI_Ec"]:
        s = sc[sc["target"] == t]
        y = s["label"].to_numpy(int); a = s["oriented"].to_numpy(float)
        sc_rows.append({"target": t, "method": "diffdock_singlechain", "n_ligands": len(s), "n_analysis": len(s), "n_active": int(y.sum()),
                        "n_decoy": int((y == 0).sum()), "auc": auc(y, a), "auc_ci_lo": np.nan, "auc_ci_hi": np.nan,
                        "top5_recall": top5(y, a), "rho_vina": np.nan, "status": "ok"})
    all_metrics.append(pd.DataFrame(sc_rows))
    metrics = pd.concat(all_metrics, ignore_index=True)
    paired = pd.concat(all_paired, ignore_index=True) if all_paired else pd.DataFrame()
    cov = coverage(long)
    metrics.to_csv(OUT / "metrics_all_methods_final.csv", index=False)
    paired.to_csv(OUT / "paired_all_methods_final.csv", index=False)
    cov.to_csv(OUT / "coverage_all_methods_final.csv", index=False)
    pooled_main = {m: pooled(metrics[metrics["method"] == m], MAIN8) for m in ["vina", "diffdock_main", "boltz2", "drugclip"]}
    pooled_main["unidock"] = pooled(metrics[metrics["method"] == "unidock"], COMMON7)
    common7 = {m: pooled(metrics[metrics["method"] == m], COMMON7) for m in ["vina", "diffdock_main", "boltz2", "drugclip", "unidock"]}
    payload = {"methods": {m: {"label": {"vina": "Vina", "unidock": "Uni-Dock", "diffdock_main": "DiffDock", "boltz2": "Boltz-2", "drugclip": "DrugCLIP"}[m],
                                   "pooled": pooled_main[m], "common7": common7[m]} for m in pooled_main},
               "pooling": {"main8": MAIN8, "unidock7": COMMON7, "common7": COMMON7},
               "coverage": cov.to_dict("records"), "preview_not_used": True}
    (OUT / "metrics_pooled_final.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    # SI tables
    metrics[metrics["method"] == "drugclip"].to_csv(OUT / "table_S29_drugclip_per_target.csv", index=False)
    metrics[metrics["method"] == "unidock"].to_csv(OUT / "table_S30_unidock_per_target.csv", index=False)
    pooled_rows = []
    for m, v in pooled_main.items():
        pooled_rows.append({"method": m, "pooled_scope": "Uni-Dock 7 targets (excl FabI)" if m == "unidock" else "main 8 targets",
                            "mean_auc": v["mean_auc"], "ci_lo": v["ci95"][0], "ci_hi": v["ci95"][1], "n_targets": v["n_targets"]})
    pd.DataFrame(pooled_rows).to_csv(OUT / "table_S31_method_pooled.csv", index=False)
    paired.to_csv(OUT / "table_S31b_method_paired.csv", index=False)
    cov.to_csv(OUT / "table_S31c_method_coverage.csv", index=False)
    metrics[metrics["method"].isin(["vina", "diffdock_main", "boltz2", "drugclip", "unidock"])].to_csv(OUT / "fig7_data.csv", index=False)
    make_figure(metrics, pooled_main, paired, cov)
    return {"metrics": len(metrics), "pooled": {m: pooled_main[m]["mean_auc"] for m in pooled_main}, "paired": len(paired), "coverage": len(cov)}


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, default=str))
