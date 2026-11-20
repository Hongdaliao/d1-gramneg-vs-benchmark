#!/usr/bin/env python
"""Amendment #12 method-axis analysis (DiffDock + single-chain + ECFP4 reference).

Outputs (all under results/trackA/method_axis/analysis):
  metrics_per_target.csv, metrics_pooled.json, paired_bootstrap.json,
  coverage_na.csv, fig7_data.csv, table_S25...S28...
Boltz-2 is reported as pending until the complete 166-row score table exists.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/trackA/method_axis"
OUT = BASE / "analysis"
OUT.mkdir(parents=True, exist_ok=True)
TARGETS = ["LpxC_Ec", "LpxC_Pa", "DHFR_Ec", "Gyrase_Ec", "GyrA_Ec", "GyrB_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]
EVAL_TARGETS = [t for t in TARGETS if t != "LepB_Ec"]
SINGLE_TARGETS = ["GyrA_Ec", "FabI_Ec", "LepB_Ec"]
N_BOOT = 2000
SEED = 42


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def auc(y, score) -> float:
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    if y.sum() == 0 or y.sum() == len(y):
        return float("nan")
    try:
        return float(roc_auc_score(y, score))
    except ValueError:
        return float("nan")


def boot_auc(y, score, rng) -> Tuple[float, float]:
    y = np.asarray(y, dtype=int); score = np.asarray(score, dtype=float)
    pos = np.where(y == 1)[0]; neg = np.where(y == 0)[0]
    if len(pos) < 2 or len(neg) < 2:
        return float("nan"), float("nan")
    draws = []
    for _ in range(N_BOOT):
        idx = np.r_[rng.choice(pos, len(pos), replace=True), rng.choice(neg, len(neg), replace=True)]
        val = auc(y[idx], score[idx])
        if np.isfinite(val):
            draws.append(val)
    return (float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))) if draws else (float("nan"), float("nan"))


def top5_recall(y, score) -> float:
    y = np.asarray(y, dtype=int); score = np.asarray(score, dtype=float)
    k = min(5, len(y))
    order = np.argsort(-score)[:k]
    return float(y[order].sum() / min(5, max(1, y.sum())))


def rho_signed(x, y) -> float:
    x = pd.Series(np.asarray(x, dtype=float)); y = pd.Series(np.asarray(y, dtype=float))
    valid = x.notna() & y.notna()
    if valid.sum() < 3:
        return float("nan")
    return float(x[valid].rank().corr(y[valid].rank()))


def load_diffdock(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    d = pd.read_csv(path)
    d["method"] = "DiffDock"
    return d


def load_subset() -> pd.DataFrame:
    sub = pd.read_csv(BASE / "subsets_diffdock.csv")
    sub["method"] = "DiffDock"
    return sub


def load_vina(target: str) -> pd.DataFrame:
    path = ROOT / "results/trackA/full" / target / f"{target}_scores.csv"
    d = pd.read_csv(path)
    d = d[d["status"].astype(str).str.lower().eq("ok")].copy()
    d["ligand_id_norm"] = d["ligand_id"].astype(str)
    d["vina_raw"] = pd.to_numeric(d["score"], errors="coerce")  # fallback if raw_score absent
    if "raw_score" in d.columns:
        d["vina_raw"] = pd.to_numeric(d["raw_score"], errors="coerce")
    return d[["ligand_id_norm", "vina_raw"]].drop_duplicates("ligand_id_norm")


def load_ecfp4(target: str) -> pd.DataFrame:
    path = ROOT / "results/trackA/ecfp4_baseline" / f"{target}.csv"
    if not path.exists():
        return pd.DataFrame()
    d = pd.read_csv(path)
    d["ligand_id_norm"] = d["ligand_id"].astype(str)
    d["ecfp4_score"] = pd.to_numeric(d["score"], errors="coerce")
    return d[["ligand_id_norm", "ecfp4_score"]].drop_duplicates("ligand_id_norm")


def load_main() -> pd.DataFrame:
    frames = []
    for target in TARGETS:
        path = BASE / "diffdock/main" / target / "scores.csv"
        if path.exists():
            d = load_diffdock(path)
            d["target"] = target
            frames.append(d)
    if not frames:
        raise FileNotFoundError("no DiffDock main score files")
    return pd.concat(frames, ignore_index=True)


def merge_main(subset: pd.DataFrame, diff: pd.DataFrame) -> pd.DataFrame:
    dd = diff[["target", "subset_id", "ligand_id", "status", "confidence_rank1", "retry_used", "n_rank_files", "complex_name"]].rename(
        columns={"status": "dd_status", "complex_name": "dd_complex_name"}
    )
    merged = subset.merge(dd, on=["target", "subset_id", "ligand_id"], how="left")
    if merged["dd_status"].isna().any():
        raise RuntimeError(f"{int(merged['dd_status'].isna().sum())} DiffDock score rows missing")
    merged["dd_ok"] = merged["dd_status"].astype(str).str.lower().eq("ok")
    merged["conf"] = pd.to_numeric(merged["confidence_rank1"], errors="coerce")
    merged["vina_score"] = pd.to_numeric(merged["vina_score"], errors="coerce")
    return merged


def target_metrics(frame: pd.DataFrame, target: str) -> List[Dict[str, Any]]:
    full = frame[frame["target"] == target].copy()
    if full.empty:
        return []
    sub = full[full["dd_ok"] & full["vina_score"].notna()].copy()
    y = sub["label"].to_numpy(dtype=int)
    dd = sub["conf"].to_numpy(dtype=float)
    vina = -sub["vina_score"].to_numpy(dtype=float)  # lower Vina score is better
    out = []
    if target != "LepB_Ec":
        for method, score in (("DiffDock", dd), ("Vina", vina)):
            lo, hi = boot_auc(y, score, np.random.default_rng(SEED + TARGETS.index(target) + (0 if method == "DiffDock" else 100)))
            out.append({"target": target, "method": method, "n_ligands": len(sub), "n_analysis": len(sub),
                        "n_active": int(y.sum()), "n_decoy": int((y == 0).sum()), "n_na_excluded": len(full) - len(sub),
                        "auc": auc(y, score), "auc_ci_lo": lo, "auc_ci_hi": hi,
                        "top5_recall": top5_recall(y, score), "rho_confidence_vina": rho_signed(dd, sub["vina_score"]),
                        "status": "paired same ligands" if method == "Vina" else "ok"})
        # ECFP4 subset reference is intentionally NOT computed: subset actives
        # overlap the temporal-train reference (68-94% for the 8 comparison
        # targets except TEM1 38%; LepB n=3 descriptive).  See Addendum 12a;
        # the full-panel test-only ECFP4 reference remains in §3.3.

    else:
        out.append({"target": target, "method": "DiffDock", "n_ligands": len(sub), "n_analysis": len(sub),
                    "n_active": int(y.sum()), "n_decoy": int((y == 0).sum()), "n_na_excluded": len(full) - len(sub),
                    "auc": float("nan"), "auc_ci_lo": float("nan"), "auc_ci_hi": float("nan"),
                    "top5_recall": float("nan"), "rho_confidence_vina": float("nan"), "status": "descriptive only (n=6)"})
    return out


def paired_main(frame: pd.DataFrame) -> List[Dict[str, Any]]:
    rows = []
    for target in EVAL_TARGETS:
        sub = frame[(frame["target"] == target) & frame["dd_ok"] & frame["vina_score"].notna()].copy()
        pairs = [g for _, g in sub.groupby("subset_id") if len(g) == 2 and set(g["label"].astype(int)) == {0, 1}]
        if not pairs:
            continue
        rng = np.random.default_rng(SEED + TARGETS.index(target))
        deltas = []
        for _ in range(N_BOOT):
            sample = pd.concat([pairs[i] for i in rng.integers(0, len(pairs), len(pairs))], ignore_index=True)
            yy = sample["label"].to_numpy(dtype=int)
            a_dd = auc(yy, sample["conf"].to_numpy(dtype=float))
            a_vina = auc(yy, -sample["vina_score"].to_numpy(dtype=float))
            if np.isfinite(a_dd) and np.isfinite(a_vina):
                deltas.append(a_dd - a_vina)
        point = auc(sub["label"], sub["conf"]) - auc(sub["label"], -sub["vina_score"])
        rows.append({"target": target, "n_pairs": len(pairs), "delta_auc_diffdock_minus_vina": point,
                     "delta_ci_lo": float(np.percentile(deltas, 2.5)), "delta_ci_hi": float(np.percentile(deltas, 97.5)),
                     "n_boot": len(deltas)})
    return rows


def singlechain_analysis() -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, Any]]:
    main = load_main()
    sc_frames = []
    for t in SINGLE_TARGETS:
        p = BASE / "diffdock/singlechain" / t / "scores.csv"
        if p.exists():
            d = pd.read_csv(p); d["target"] = t; sc_frames.append(d)
    sc = pd.concat(sc_frames, ignore_index=True) if sc_frames else pd.DataFrame()
    rows, cov_rows = [], []
    summary = {"status": "ok" if not sc.empty else "pending", "targets": {}}
    for t in SINGLE_TARGETS:
        expected = len(pd.read_csv(BASE / "subsets_diffdock.csv").query("target == @t"))
        m = main[(main["target"] == t) & (main["status"].astype(str).str.lower() == "ok")].copy()
        s = sc[(sc["target"] == t) & (sc["status"].astype(str).str.lower() == "ok")].copy()
        merged = m[["subset_id", "ligand_id", "label", "confidence_rank1"]].merge(
            s[["subset_id", "ligand_id", "confidence_rank1"]], on=["subset_id", "ligand_id"], how="inner", suffixes=("_main", "_single"))
        merged["conf_main"] = pd.to_numeric(merged["confidence_rank1_main"], errors="coerce")
        merged["conf_single"] = pd.to_numeric(merged["confidence_rank1_single"], errors="coerce")
        merged = merged.dropna(subset=["conf_main", "conf_single"])
        y = merged["label"].to_numpy(dtype=int)
        a_main = auc(y, merged["conf_main"].to_numpy(dtype=float))
        a_single = auc(y, merged["conf_single"].to_numpy(dtype=float))
        delta = a_single - a_main
        # paired bootstrap
        pair_groups = [g for _, g in merged.groupby("subset_id") if len(g) == 2 and set(g["label"].astype(int)) == {0, 1}]
        rng = np.random.default_rng(SEED + 333)
        deltas = []
        for _ in range(N_BOOT):
            sample = pd.concat([pair_groups[i] for i in rng.integers(0, len(pair_groups), len(pair_groups))], ignore_index=True)
            yy = sample["label"].to_numpy(dtype=int)
            am = auc(yy, sample["conf_main"].to_numpy(dtype=float))
            as_ = auc(yy, sample["conf_single"].to_numpy(dtype=float))
            if np.isfinite(am) and np.isfinite(as_):
                deltas.append(as_ - am)
        if t != "LepB_Ec" and len(merged) >= 4:
            rows.append({"target": t, "n_expected": expected, "n_pairs": len(pair_groups), "n_analysis": len(merged),
                         "auc_multichain": a_main, "auc_singlechain": a_single, "delta_singlechain_minus_multichain": delta,
                         "delta_ci_lo": float(np.percentile(deltas, 2.5)), "delta_ci_hi": float(np.percentile(deltas, 97.5)),
                         "n_boot": len(deltas), "status": "same pairs; paired bootstrap"})
            summary["targets"][t] = {"auc_multichain": a_main, "auc_singlechain": a_single, "delta": delta,
                                     "delta_ci": [float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5))]}
        else:
            rows.append({"target": t, "n_expected": expected, "n_pairs": len(pair_groups) if "pair_groups" in locals() else 0,
                         "n_analysis": len(merged), "auc_multichain": a_main, "auc_singlechain": a_single,
                         "delta_singlechain_minus_multichain": delta, "delta_ci_lo": float("nan"), "delta_ci_hi": float("nan"),
                         "n_boot": 0, "status": "descriptive only (n<4)"})
        cov_rows.append({"target": t, "n_expected": expected,
                         "n_multichain_ok": int(((main["target"] == t) & (main["status"].astype(str).str.lower() == "ok")).sum()),
                         "n_singlechain_ok": int(((sc["target"] == t) & (sc["status"].astype(str).str.lower() == "ok")).sum()) if not sc.empty else 0,
                         "n_na": expected - int(len(merged)), "note": "same pairs; paired bootstrap" if not sc.empty else "single-chain arm pending"})
    return pd.DataFrame(rows), pd.DataFrame(cov_rows), summary


def coverage_main(frame: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    out = frame.copy()
    out["attempts"] = 1 + pd.to_numeric(out["retry_used"], errors="coerce").fillna(0).astype(int)
    out["failed_at"] = ""
    failed = ~out["dd_ok"]
    no_rank = pd.to_numeric(out["n_rank_files"], errors="coerce").fillna(0).eq(0)
    out.loc[failed & no_rank, "failed_at"] = "no_rank_files"
    out.loc[failed & ~no_rank, "failed_at"] = "diffdock_inference_after_retry"
    keep = ["target", "dd_complex_name", "subset_id", "ligand_id", "label", "dd_status", "confidence_rank1", "attempts", "retry_used", "n_rank_files", "failed_at"]
    cov = out[keep].rename(columns={"dd_status": "status", "dd_complex_name": "complex_name"})
    summary = cov.groupby("target")["status"].value_counts().unstack(fill_value=0).to_dict("index")
    return cov, summary


def boltz_analysis():
    """Boltz-2 `ligand_iptm` analysis; returns None metrics until 166/166 rows exist."""
    expected = pd.read_csv(BASE / "subsets_boltz2.csv")
    score_path = BASE / "boltz2_scores.csv"
    if not score_path.exists():
        return None, None, None, {"status": "pending", "n_observed": 0, "n_expected": len(expected), "reason": "no score table"}
    observed = pd.read_csv(score_path)
    if len(observed) < len(expected):
        return None, None, None, {"status": "pending", "n_observed": len(observed), "n_expected": len(expected), "reason": "score table incomplete"}
    sc = observed[["target", "ligand_id", "status", "ligand_iptm", "confidence_score", "attempts"]].rename(
        columns={"status": "bz_status", "confidence_score": "bz_confidence"})
    m = expected.merge(sc, on=["target", "ligand_id"], how="left")
    if m["bz_status"].isna().any():
        raise RuntimeError("Boltz rows missing after merge")
    m["bz_ok"] = m["bz_status"].astype(str).str.lower().eq("ok")
    m["iptm"] = pd.to_numeric(m["ligand_iptm"], errors="coerce")
    m["vina"] = pd.to_numeric(m["vina_score"], errors="coerce")
    cov = m[["target", "ligand_id", "label", "bz_status", "ligand_iptm", "attempts"]].rename(columns={"bz_status": "status"}).copy()
    cov["failed_at"] = ""
    cov.loc[(cov["status"].astype(str).str.lower() != "ok") & cov["attempts"].fillna(0).le(2), "failed_at"] = "boltz_inference_after_attempt"
    metrics_rows, paired_rows, conf = [], [], []
    for t in TARGETS:
        sub = m[(m["target"] == t) & m["bz_ok"] & m["iptm"].notna()].copy()
        if sub.empty:
            continue
        y = sub["label"].to_numpy(dtype=int); sc_iptm = sub["iptm"].to_numpy(dtype=float)
        rho = rho_signed(sub["iptm"], sub["vina"])
        conf.append({"target": t, "n": len(sub), "mean_ligand_iptm": float(np.mean(sc_iptm)),
                     "median_ligand_iptm": float(np.median(sc_iptm)), "sd_ligand_iptm": float(np.std(sc_iptm, ddof=1)) if len(sub) > 1 else float("nan")})
        if t == "LepB_Ec":
            metrics_rows.append({"target": t, "method": "Boltz-2", "n_ligands": len(sub), "n_analysis": len(sub),
                                 "n_active": int(y.sum()), "n_decoy": int((y == 0).sum()), "n_na_excluded": int(len(m[m["target"] == t]) - len(sub)),
                                 "auc": float("nan"), "auc_ci_lo": float("nan"), "auc_ci_hi": float("nan"),
                                 "top5_recall": float("nan"), "rho_confidence_vina": rho, "status": "descriptive only (n=6)"})
            continue
        lo, hi = boot_auc(y, sc_iptm, np.random.default_rng(SEED + 800 + TARGETS.index(t)))
        metrics_rows.append({"target": t, "method": "Boltz-2", "n_ligands": len(sub), "n_analysis": len(sub),
                             "n_active": int(y.sum()), "n_decoy": int((y == 0).sum()), "n_na_excluded": int(len(m[m["target"] == t]) - len(sub)),
                             "auc": auc(y, sc_iptm), "auc_ci_lo": lo, "auc_ci_hi": hi,
                             "top5_recall": float("nan"), "rho_confidence_vina": rho, "status": "ok"})
        pairs = [g for _, g in sub.groupby("subset_id") if len(g) == 2 and set(g["label"].astype(int)) == {0, 1}]
        if pairs:
            rng = np.random.default_rng(SEED + 900 + TARGETS.index(t)); deltas = []
            for _ in range(N_BOOT):
                sample = pd.concat([pairs[i] for i in rng.integers(0, len(pairs), len(pairs))], ignore_index=True)
                yy = sample["label"].to_numpy(dtype=int)
                aa = auc(yy, sample["iptm"].to_numpy(dtype=float))
                av = auc(yy, -sample["vina"].to_numpy(dtype=float))
                if np.isfinite(aa) and np.isfinite(av):
                    deltas.append(aa - av)
            paired_rows.append({"target": t, "n_pairs": len(pairs),
                                "delta_auc_boltz_minus_vina": auc(y, sc_iptm) - auc(y, -sub["vina"].to_numpy(dtype=float)),
                                "delta_ci_lo": float(np.percentile(deltas, 2.5)) if deltas else float("nan"),
                                "delta_ci_hi": float(np.percentile(deltas, 97.5)) if deltas else float("nan"),
                                "boltz_auc": auc(y, sc_iptm), "vina_auc": auc(y, -sub["vina"].to_numpy(dtype=float)),
                                "n_boot": len(deltas)})
    summary = {
        "status": "ok", "n_expected": len(expected), "n_observed": len(observed), "n_ok": int(m["bz_ok"].sum()),
        "n_na": int((~m["bz_ok"]).sum()), "metrics": metrics_rows, "paired": paired_rows,
        "confidence_distribution": conf,
    }
    if metrics_rows:
        target_auc = [r["auc"] for r in metrics_rows if r["target"] in EVAL_TARGETS and np.isfinite(r["auc"])]
        if target_auc:
            rng = np.random.default_rng(SEED); arr = np.array(target_auc)
            draws = [float(np.mean(arr[rng.integers(0, len(arr), len(arr))])) for _ in range(N_BOOT)]
            summary["pooled_auc_mean"] = float(np.mean(arr)); summary["pooled_auc_ci95"] = [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))]
    if paired_rows:
        target_deltas = np.array([r["delta_auc_boltz_minus_vina"] for r in paired_rows], dtype=float)
        rng = np.random.default_rng(SEED + 5678)
        draws = [float(np.mean(target_deltas[rng.integers(0, len(target_deltas), len(target_deltas))])) for _ in range(N_BOOT)]
        summary["pooled_delta_boltz_minus_vina"] = float(np.mean(target_deltas))
        summary["pooled_delta_ci95"] = [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))]
    return pd.DataFrame(metrics_rows), pd.DataFrame(paired_rows), cov, summary


def make_figure(metrics: pd.DataFrame, pooled: Dict[str, Any], single: pd.DataFrame, boltz: Optional[pd.DataFrame] = None) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 5, figsize=(16.5, 3.5))
    targets = EVAL_TARGETS
    vina = metrics[(metrics.method == "Vina") & metrics.target.isin(targets)].set_index("target").reindex(targets)
    dd = metrics[(metrics.method == "DiffDock") & metrics.target.isin(targets)].set_index("target").reindex(targets)
    for target in targets:
        axes[0].plot([0, 1], [vina.loc[target, "auc"], dd.loc[target, "auc"]], color="0.6", lw=0.8)
        axes[0].scatter([0], [vina.loc[target, "auc"]], color="#225ea8", s=20)
        axes[0].scatter([1], [dd.loc[target, "auc"]], color="#cb181d", s=20)
        axes[0].annotate(target, (1, dd.loc[target, "auc"]), fontsize=5.5, xytext=(3, 0), textcoords="offset points")
    axes[0].set_xticks([0, 1], ["Vina", "DiffDock"]); axes[0].set_ylabel("AUC"); axes[0].set_title("A. Paired per-target AUC")
    axes[0].set_xlim(-0.3, 1.8)
    axes[1].bar([0], [pooled["vina_auc_mean"]], color="#225ea8", width=0.5)
    axes[1].bar([1], [pooled["diffdock_auc_mean"]], color="#cb181d", width=0.5)
    axes[1].set_xticks([0, 1], ["Vina", "DiffDock"]); axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("Pooled mean AUC"); axes[1].set_title("B. Pooled AUC (8 targets)")
    axes[1].text(0.02, 0.98, f"DiffDock {pooled['diffdock_auc_mean']:.3f}\nCI {pooled['diffdock_auc_ci_lo']:.3f}-{pooled['diffdock_auc_ci_hi']:.3f}",
                 transform=axes[1].transAxes, va="top", fontsize=6.5)
    rhos = dd["rho_confidence_vina"].astype(float)
    axes[2].barh(np.arange(len(rhos)), rhos.values, color=["#cb181d" if v < 0 else "#225ea8" for v in rhos.values])
    axes[2].axvline(0, color="0.5", lw=0.7); axes[2].set_yticks(np.arange(len(rhos))); axes[2].set_yticklabels(rhos.index, fontsize=6)
    axes[2].set_xlabel("Spearman rho (confidence vs Vina)"); axes[2].set_title("C. Rank agreement")
    # D single-chain
    if single is not None and len(single):
        sub = single[single.status.astype(str).str.contains("paired")].copy()
        x = np.arange(len(sub))
        axes[3].bar(x - 0.18, sub["auc_multichain"], width=0.36, color="#225ea8", label="multichain")
        axes[3].bar(x + 0.18, sub["auc_singlechain"], width=0.36, color="#cb181d", label="single-chain")
        axes[3].set_xticks(x); axes[3].set_xticklabels(sub["target"], rotation=20, fontsize=6)
        axes[3].set_ylim(0, 1); axes[3].set_ylabel("AUC"); axes[3].set_title("D. Single-chain sensitivity")
        axes[3].legend(frameon=False, fontsize=6)
    else:
        axes[3].axis("off"); axes[3].text(0.5, 0.5, "D. Single-chain\npending", ha="center", va="center")
    if boltz is not None and len(boltz):
        b = boltz[(boltz["method"] == "Boltz-2") & boltz["target"].isin(EVAL_TARGETS)].set_index("target").reindex(EVAL_TARGETS)
        x = np.arange(len(b))
        axes[4].bar(x - 0.18, vina.loc[b.index, "auc"], width=0.36, color="#225ea8", label="Vina")
        axes[4].bar(x + 0.18, b["auc"], width=0.36, color="#cb181d", label="Boltz-2")
        axes[4].set_xticks(x); axes[4].set_xticklabels(b.index, rotation=45, ha="right", fontsize=5.5)
        axes[4].set_ylim(0, 1); axes[4].set_ylabel("AUC"); axes[4].set_title("E. Boltz-2 (ligand_iptm)")
        bm = pooled.get("methods", {}).get("Boltz-2", {})
        if bm:
            axes[4].text(0.02, 0.98, f"mean {bm.get('mean_auc', float('nan')):.3f}", transform=axes[4].transAxes, va="top", fontsize=6)
        axes[4].legend(frameon=False, fontsize=5.5)
    else:
        axes[4].axis("off"); axes[4].text(0.5, 0.5, "Boltz-2 axis\npending\n(no complete score table)", ha="center", va="center", fontsize=8)
        axes[4].set_title("E. Boltz-2")
    fig.tight_layout(); fig.savefig(ROOT / "figures/fig7_method_axis_draft.png", dpi=300, bbox_inches="tight")
    fig.savefig(ROOT / "figures/fig7_method_axis_draft.pdf", bbox_inches="tight"); plt.close(fig)


def main() -> Dict[str, Any]:
    subset = load_subset(); diff = load_main(); merged = merge_main(subset, diff)
    # ECFP4 leak diagnostic (Addendum 12a): subset actives in temporal-train reference
    act = subset[subset["label"] == 1].copy()
    overlap = []
    for t in TARGETS:
        g = act[act["target"] == t]
        n = len(g); tr = int((g["temporal_split"].astype(str).str.lower() == "train").sum())
        overlap.append({"target": t, "n_subset_actives": n, "train_reference_overlap": tr,
                        "overlap_pct": round(100.0 * tr / n, 1) if n else float("nan")})
    pd.DataFrame(overlap).to_csv(OUT / "ecfp4_leakage_overlap.csv", index=False)
    metrics_rows = []
    for t in TARGETS:
        metrics_rows += target_metrics(merged, t)
    metrics = pd.DataFrame(metrics_rows)
    boltz_metrics, boltz_paired, boltz_cov, boltz_summary = boltz_analysis()
    if boltz_metrics is not None:
        metrics = pd.concat([metrics, boltz_metrics], ignore_index=True)
        boltz_cov.to_csv(OUT / "boltz2_coverage_na.csv", index=False)
        boltz_cov.to_csv(OUT / "table_S27b_boltz2_coverage.csv", index=False)
        boltz_paired.to_csv(OUT / "boltz2_paired_bootstrap.csv", index=False)
        boltz_paired.to_csv(OUT / "table_S26b_boltz2_paired.csv", index=False)
    paired = paired_main(merged)
    cov, cov_summary = coverage_main(merged)
    metrics.to_csv(OUT / "metrics_per_target.csv", index=False)
    cov.to_csv(OUT / "coverage_na.csv", index=False)
    pd.DataFrame(paired).to_csv(OUT / "paired_bootstrap_per_target.csv", index=False)
    single, single_cov, single_summary = singlechain_analysis()
    if len(single) and len(single_cov):
        single = single.merge(single_cov.drop(columns=["n_expected"]), on="target", how="left")
    single.to_csv(OUT / "singlechain_comparison.csv", index=False)
    single_cov.to_csv(OUT / "singlechain_coverage.csv", index=False)
    single_cov.assign(arm="diffdock_singlechain").to_csv(OUT / "singlechain_status.csv", index=False)
    # single-chain NA rows (explicit, not hidden)
    sc_frames = []
    for t in SINGLE_TARGETS:
        p = BASE / "diffdock/singlechain" / t / "scores.csv"
        if p.exists():
            d = pd.read_csv(p); d["target"] = t; sc_frames.append(d)
    sc_na = pd.concat(sc_frames, ignore_index=True) if sc_frames else pd.DataFrame()
    if len(sc_na):
        sc_na = sc_na[sc_na["status"].astype(str).str.lower() != "ok"].copy()
        sc_na["attempts"] = 1 + pd.to_numeric(sc_na["retry_used"], errors="coerce").fillna(0).astype(int)
        sc_na.to_csv(OUT / "singlechain_coverage_na.csv", index=False)
    else:
        pd.DataFrame().to_csv(OUT / "singlechain_coverage_na.csv", index=False)

    # pooled means/CI (8 targets)
    p_rows = {}
    for method in ("DiffDock", "Vina", "Boltz-2"):
        vals = metrics[(metrics.method == method) & metrics.target.isin(EVAL_TARGETS)].set_index("target")["auc"]
        vals = vals.dropna()
        if vals.empty:
            continue
        rng = np.random.default_rng(SEED)
        arr = vals.to_numpy(dtype=float)
        draws = [float(np.mean(arr[rng.integers(0, len(arr), len(arr))])) for _ in range(N_BOOT)]
        p_rows[method] = {"mean_auc": float(vals.mean()), "ci95": [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))], "targets": vals.to_dict()}
    target_deltas = np.array([r["delta_auc_diffdock_minus_vina"] for r in paired], dtype=float)
    rng_delta = np.random.default_rng(SEED)
    delta_draws = [float(np.mean(target_deltas[rng_delta.integers(0, len(target_deltas), len(target_deltas))])) for _ in range(N_BOOT)]
    pooled = {
        "n_targets": len(EVAL_TARGETS), "targets": EVAL_TARGETS, "excluded": "LepB_Ec (n=6; descriptive only)",
        "methods": p_rows,
        "diffdock_auc_mean": p_rows.get("DiffDock", {}).get("mean_auc"),
        "diffdock_auc_ci_lo": p_rows.get("DiffDock", {}).get("ci95", [None, None])[0],
        "diffdock_auc_ci_hi": p_rows.get("DiffDock", {}).get("ci95", [None, None])[1],
        "vina_auc_mean": p_rows.get("Vina", {}).get("mean_auc"),
        "delta_auc_diffdock_minus_vina_mean": float(np.mean(target_deltas)) if len(target_deltas) else float("nan"),
        "delta_auc_ci_lo": float(np.percentile(delta_draws, 2.5)) if delta_draws else float("nan"),
        "delta_auc_ci_hi": float(np.percentile(delta_draws, 97.5)) if delta_draws else float("nan"),
        "ci_method": "target-level bootstrap over 8 targets (2,000 draws)",
        "per_target_ci_method": "ligand-level bootstrap within target (2,000 draws)",
        "singlechain": single_summary,
        "boltz2": boltz_summary,
    }
    metadata = {
        "amendment_sha256": sha256(ROOT / "notes/preregistration_amendment12_method_axis_draft.md"),
        "subset_diffdock_sha256": sha256(BASE / "subsets_diffdock.csv"),
        "subset_boltz2_sha256": sha256(BASE / "subsets_boltz2.csv"),
        "n_expected": len(subset), "n_diffdock_ok": int(merged["dd_ok"].sum()), "n_diffdock_na": int((~merged["dd_ok"]).sum()),
        "n_boot": N_BOOT, "seed": SEED,
        "diffdock_main_scores": str(BASE / "diffdock/main"),
        "singlechain_scores": str(BASE / "diffdock/singlechain"),
    }
    (OUT / "metrics_pooled.json").write_text(json.dumps({"task": "Amendment #12", "status": "DRY_RUN_DIFFDOCK_PLUS_SINGLECHAIN", "metadata": metadata, "pooled": pooled, "metrics": metrics.to_dict("records")}, indent=2), encoding="utf-8")
    (OUT / "paired_bootstrap.json").write_text(json.dumps({"task": "Amendment #12", "status": "DRY_RUN", "metadata": metadata, "diffdock_minus_vina": paired, "singlechain_minus_multichain": single.to_dict("records"), "boltz_minus_vina": (boltz_paired.to_dict("records") if boltz_paired is not None else {"status": "pending"})}, indent=2), encoding="utf-8")
    # fig data & SI skeletons
    fig_data = metrics.copy()
    if len(single):
        extra = []
        for _, r in single.iterrows():
            extra.append({"target": r["target"], "method": "DiffDock_multichain", "n_ligands": r.get("n_expected"),
                          "n_analysis": r.get("n_analysis"), "n_active": None, "n_decoy": None, "n_na_excluded": r.get("n_na"),
                          "auc": r.get("auc_multichain"), "auc_ci_lo": None, "auc_ci_hi": None,
                          "top5_recall": None, "rho_confidence_vina": None, "status": "single-chain panel"})
            extra.append({"target": r["target"], "method": "DiffDock_singlechain", "n_ligands": r.get("n_expected"),
                          "n_analysis": r.get("n_analysis"), "n_active": None, "n_decoy": None, "n_na_excluded": r.get("n_na"),
                          "auc": r.get("auc_singlechain"), "auc_ci_lo": None, "auc_ci_hi": None,
                          "top5_recall": None, "rho_confidence_vina": None, "status": "single-chain panel"})
        fig_data = pd.concat([fig_data, pd.DataFrame(extra)], ignore_index=True)
    fig_data.to_csv(OUT / "fig7_data.csv", index=False)
    metrics.to_csv(OUT / "table_S25_method_axis_per_target.csv", index=False)
    pd.DataFrame(paired).to_csv(OUT / "table_S26_method_axis_paired.csv", index=False)
    cov.to_csv(OUT / "table_S27_method_axis_coverage.csv", index=False)
    single.to_csv(OUT / "table_S28_method_axis_singlechain.csv", index=False)
    # S25 footnote: ECFP4 subset reference withdrawn for leakage (Addendum 12a)
    overlap_df = pd.read_csv(OUT / "ecfp4_leakage_overlap.csv")
    lines = ["# Table S25 footnote — ECFP4 subset reference withdrawn", "",
             "Addendum 12a: ECFP4 test-only was not recomputed on the DiffDock balanced subset because",
             "subset actives overlap the temporal-train reference used by the test-only ECFP4 protocol.",
             "The 2D reference remains the full-panel test-only result in §3.3 (median AUC 0.775).", "",
             "| target | subset actives (n) | in temporal-train reference (n) | overlap |",
             "|---|---:|---:|---:|"]
    for _, r in overlap_df.iterrows():
        lines.append(f"| {r['target']} | {int(r['n_subset_actives'])} | {int(r['train_reference_overlap'])} | {r['overlap_pct']:.1f}% |")
    lines += ["", "Across the eight comparison targets overlap is 68–94% (GyrB 68% to LpxC 94%); TEM1 is 38% and LepB 100% (n=3; descriptive).", ""]
    (OUT / "table_S25_method_axis_per_target_README.md").write_text("\n".join(lines), encoding="utf-8")
    (OUT / "table_S26_method_axis_paired_README.md").write_text(
        "# Table S26 footnote\n\nPaired comparisons are DiffDock−Vina (main) and singlechain−multichain (S28). "
        "No ECFP4 subset reference is included: Addendum 12a withdrew it because of temporal-train overlap "
        "(see S25 footnote); the full-panel test-only ECFP4 result is reported in §3.3.\n", encoding="utf-8")
    make_figure(metrics, pooled, single, boltz_metrics)
    return {"out": str(OUT), "pooled": pooled, "n_metrics": len(metrics), "figure": str(ROOT / "figures/fig7_method_axis_draft.png")}


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, default=str))
