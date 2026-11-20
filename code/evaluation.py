"""Evaluation utilities: discrimination, calibration, ablations, sensitivity.

Implements the metrics required by the reporting checklist:

* AUROC, Brier score, ECE + reliability curve            (checklist M1/C3)
* Spearman rho(MAES, RMSD) with the expected 0.5-0.8 band (checklist C4)
* grouped-by-target cross-validation                    (checklist C2)
* ablation (full / no M / no S / no strain / uncalibrated) (checklist A1)
* entropy + Dirichlet random weight sensitivity, Kendall tau (checklist A2/A3)
* five-class failure taxonomy                           (checklist F1/F2)

Pure numpy/scipy so that it also runs where sklearn is unavailable.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from .calibration import MAESModel, build_evaluation_labels, make_calibrator
from .combine import combine_subscores, dirichlet_weights, entropy_weights
from .config import MAESConfig

__all__ = [
    "ablation_study",
    "auroc",
    "brier_score",
    "ece",
    "evaluate_predictions",
    "failure_count_matrix",
    "failure_taxonomy",
    "grouped_cv",
    "kendall_tau",
    "reliability_curve",
    "spearman_rho",
    "weight_sensitivity",
]


# ------------------------------------------------------------------ basic metrics
def auroc(y_true: Sequence[int], y_score: Sequence[float]) -> float:
    """Area under the ROC curve (Mann-Whitney U form, ties handled)."""
    y = np.asarray(y_true, dtype=int)
    s = np.asarray(y_score, dtype=float)
    if y.size != s.size:
        raise ValueError("y_true and y_score differ in length")
    mask = np.isfinite(s)
    y, s = y[mask], s[mask]
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(s.size, dtype=float)
    sorted_s = s[order]
    i = 0
    while i < s.size:
        j = i
        while j + 1 < s.size and sorted_s[j + 1] == sorted_s[i]:
            j += 1
        avg_rank = 0.5 * (i + j) + 1.0
        ranks[order[i:j + 1]] = avg_rank
        i = j + 1
    u = ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2.0
    return float(u / (n_pos * n_neg))


def brier_score(y_true: Sequence[int], p: Sequence[float]) -> float:
    y = np.asarray(y_true, dtype=float)
    p = np.clip(np.asarray(p, dtype=float), 0.0, 1.0)
    return float(np.mean((p - y) ** 2))


def reliability_curve(y_true: Sequence[int], p: Sequence[float], n_bins: int = 10) -> pd.DataFrame:
    """Equal-width reliability diagram data."""
    y = np.asarray(y_true, dtype=float)
    p = np.clip(np.asarray(p, dtype=float), 0.0, 1.0)
    bins = np.clip(np.floor(p * n_bins).astype(int), 0, n_bins - 1)
    rows = []
    for b in range(n_bins):
        mask = bins == b
        count = int(mask.sum())
        rows.append(
            {
                "bin": b,
                "bin_lo": b / n_bins,
                "bin_hi": (b + 1) / n_bins,
                "count": count,
                "mean_pred": float(p[mask].mean()) if count else float("nan"),
                "frac_pos": float(y[mask].mean()) if count else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def ece(y_true: Sequence[int], p: Sequence[float], n_bins: int = 10) -> float:
    """Expected calibration error (equal-width bins, count weighted)."""
    table = reliability_curve(y_true, p, n_bins=n_bins)
    total = float(table["count"].sum())
    if total <= 0:
        return float("nan")
    gaps = (table["mean_pred"] - table["frac_pos"]).abs().fillna(0.0)
    return float((gaps * table["count"]).sum() / total)


def _rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=float)
    sorted_v = values[order]
    i = 0
    while i < values.size:
        j = i
        while j + 1 < values.size and sorted_v[j + 1] == sorted_v[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return ranks


def spearman_rho(a: Sequence[float], b: Sequence[float]) -> float:
    """Spearman rank correlation (numpy implementation, ties averaged)."""
    x = np.asarray(a, dtype=float)
    y = np.asarray(b, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if x.size < 3:
        return float("nan")
    rx = _rankdata(x) - (_rankdata(x).mean())
    ry = _rankdata(y) - (_rankdata(y).mean())
    denom = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    if denom == 0:
        return float("nan")
    return float((rx * ry).sum() / denom)


def kendall_tau(a: Sequence[float], b: Sequence[float]) -> float:
    """Kendall tau-b (via scipy when available, else O(n^2) fallback)."""
    x = np.asarray(a, dtype=float)
    y = np.asarray(b, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if x.size < 3:
        return float("nan")
    try:
        from scipy.stats import kendalltau

        return float(kendalltau(x, y).statistic)
    except Exception:
        n = x.size
        concordant = discordant = 0
        for i in range(n - 1):
            for j in range(i + 1, n):
                dx = x[i] - x[j]
                dy = y[i] - y[j]
                if dx * dy > 0:
                    concordant += 1
                elif dx * dy < 0:
                    discordant += 1
        total = n * (n - 1) / 2
        return float((concordant - discordant) / total)


# ------------------------------------------------------------------ model scoring
def evaluate_predictions(
    df: pd.DataFrame,
    model: MAESModel,
    extra_scores: Optional[Mapping[str, Sequence[float]]] = None,
) -> Dict[str, Any]:
    """Evaluate a fitted model on an independent label set (never the calib label)."""
    cfg = model.config
    maes = model.predict(df)
    eval_primary = build_evaluation_labels(df, cfg.labels, which="primary")
    eval_secondary = build_evaluation_labels(df, cfg.labels, which="secondary")
    n_bins = 10

    result: Dict[str, Any] = {
        "n": int(len(df)),
        "n_eval_positive_primary": int(eval_primary.sum()),
        "n_eval_positive_secondary": int(eval_secondary.sum()),
        "auroc_primary": auroc(eval_primary, maes),
        "auroc_secondary": auroc(eval_secondary, maes),
        "brier_primary": brier_score(eval_primary, maes),
        "ece_primary": ece(eval_primary, maes, n_bins=n_bins),
        "mean_maes": float(np.mean(maes)),
    }
    if "rmsd" in df.columns:
        result["rho_maes_rmsd"] = spearman_rho(maes, df["rmsd"].to_numpy(dtype=float))
        # C4: MAES must carry information beyond RMSD.  |rho| in [0.5, 0.8] is the
        # target band; |rho| > 0.95 means MAES is an RMSD proxy (attack A2).
        abs_rho = abs(result["rho_maes_rmsd"]) if np.isfinite(result["rho_maes_rmsd"]) else float("nan")
        result["abs_rho_maes_rmsd"] = abs_rho
        result["rho_ok_band"] = bool(np.isfinite(abs_rho) and 0.5 <= abs_rho <= 0.8)
        result["rho_rmsd_proxy_warning"] = bool(np.isfinite(abs_rho) and abs_rho > 0.95)
    if "plif_recovery" in df.columns:
        result["rho_maes_plif"] = spearman_rho(maes, df["plif_recovery"].to_numpy(dtype=float))
    if extra_scores:
        for name, scores in extra_scores.items():
            result[f"auroc_{name}"] = auroc(eval_primary, scores)
    return result


def _clone_with(df: pd.DataFrame, overrides: Mapping[str, Any]) -> pd.DataFrame:
    out = df.copy()
    for column, value in overrides.items():
        out[column] = value
    return out


def ablation_study(df: pd.DataFrame, config: MAESConfig) -> pd.DataFrame:
    """Fit/evaluate the full model and its ablations on the same data.

    Variants: ``full``, ``no_M``, ``no_S``, ``no_strain`` (V := V_PB only when a
    ``pb_pass_rate`` column exists), ``uncalibrated`` (raw score as probability).
    """
    if "V" not in df.columns:
        raise KeyError("ablation_study expects pre-computed V/M/S columns")
    eval_primary = build_evaluation_labels(df, config.labels, which="primary")
    rows: list = []

    def _run(name: str, frame: pd.DataFrame, cfg: MAESConfig) -> None:
        model = MAESModel(config=cfg).fit(frame)
        scores = model.predict(frame)
        rows.append(
            {
                "variant": name,
                "auroc": auroc(eval_primary, scores),
                "brier": brier_score(eval_primary, scores),
                "ece": ece(eval_primary, scores),
                "rho_rmsd": spearman_rho(scores, frame["rmsd"].to_numpy(dtype=float))
                if "rmsd" in frame.columns
                else float("nan"),
                "n": int(len(frame)),
            }
        )

    _run("full", df, config)

    neutral = df.copy()
    neutral["M"] = 1.0
    _run("no_M", neutral, config)

    neutral_s = df.copy()
    neutral_s["S"] = 1.0
    _run("no_S", neutral_s, config)

    if "pb_pass_rate" in df.columns:
        no_strain = df.copy()
        no_strain["V"] = np.clip(df["pb_pass_rate"].to_numpy(dtype=float), 0.0, 1.0)
        _run("no_strain", no_strain, config)
    else:
        _run("no_strain", df, config)

    cfg_uncal = MAESConfig.from_dict(config.to_dict())
    cfg_uncal.calibrator.kind = "none"
    _run("uncalibrated", df, cfg_uncal)
    return pd.DataFrame(rows)


def grouped_cv(
    df: pd.DataFrame,
    config: MAESConfig,
    group_col: str = "system_id",
    n_splits: int = 5,
) -> pd.DataFrame:
    """Grouped (by target/system) K-fold CV of MAES discrimination and calibration.

    The calibration label is used for fitting within each fold; the *primary
    evaluation* label is used only for scoring, preserving label independence.
    """
    if group_col not in df.columns:
        raise KeyError(f"missing group column {group_col!r}")
    groups = df[group_col].astype(str).to_numpy()
    unique = np.array(sorted(set(groups)))
    if unique.size < 2:
        raise ValueError("need at least 2 groups for grouped CV")
    k = int(min(max(2, n_splits), unique.size))
    folds = [unique[i::k] for i in range(k)]
    rows = []
    for fold_id, held_out in enumerate(folds):
        test_mask = np.isin(groups, held_out)
        train_df = df.loc[~test_mask]
        test_df = df.loc[test_mask]
        if len(train_df) == 0 or len(test_df) == 0:
            continue
        try:
            model = MAESModel(config=config).fit(train_df)
        except Exception as exc:  # single-class fold etc.
            rows.append({"fold": fold_id, "held_out": ",".join(map(str, held_out)), "error": str(exc)})
            continue
        y_test = build_evaluation_labels(test_df, config.labels, which="primary")
        scores = model.predict(test_df)
        rows.append(
            {
                "fold": fold_id,
                "held_out": ",".join(map(str, held_out)),
                "n_train": int(len(train_df)),
                "n_test": int(len(test_df)),
                "auroc": auroc(y_test, scores),
                "brier": brier_score(y_test, scores),
                "ece": ece(y_test, scores),
                "error": "",
            }
        )
    out = pd.DataFrame(rows)
    merged_auroc = out["auroc"].dropna()
    if len(merged_auroc):
        out.attrs["mean_auroc"] = float(merged_auroc.mean())
        out.attrs["sd_auroc"] = float(merged_auroc.std(ddof=1)) if len(merged_auroc) > 1 else 0.0
    return out


def weight_sensitivity(
    df: pd.DataFrame,
    config: MAESConfig,
    base_weights: Optional[Mapping[str, float]] = None,
    n_samples: int = 1000,
    alpha: float = 1.0,
    seed: int = 42,
    top_k: int = 3,
) -> Dict[str, Any]:
    """Ranking stability under Dirichlet-random sub-score weights (checklist A2)."""
    base_weights = dict(base_weights or config.scores.weights)
    base_model = MAESModel(config=config, weights=base_weights).fit(df)
    base_rank = base_model.predict(df)
    base_top = set(pd.Series(base_rank).nlargest(top_k).index)

    samples = dirichlet_weights(n_samples=n_samples, n_components=3, alpha=alpha, seed=seed)
    taus, overlaps = [], []
    for w_vec in samples:
        weights = {"V": float(w_vec[0]), "M": float(w_vec[1]), "S": float(w_vec[2])}
        try:
            cfg = MAESConfig.from_dict(config.to_dict())
            model = MAESModel(config=cfg, weights=weights).fit(df)
            scores = model.predict(df)
        except Exception:
            continue
        taus.append(kendall_tau(base_rank, scores))
        top = set(pd.Series(scores).nlargest(top_k).index)
        overlaps.append(len(top & base_top) / max(1, top_k))
    taus_arr = np.asarray([t for t in taus if np.isfinite(t)], dtype=float)
    return {
        "n_samples": int(len(taus_arr)),
        "base_weights": {k: float(v) for k, v in base_weights.items()},
        "entropy_weights": entropy_weights(df),
        "kendall_tau_mean": float(taus_arr.mean()) if taus_arr.size else float("nan"),
        "kendall_tau_p05": float(np.quantile(taus_arr, 0.05)) if taus_arr.size else float("nan"),
        "kendall_tau_p95": float(np.quantile(taus_arr, 0.95)) if taus_arr.size else float("nan"),
        f"top{top_k}_overlap_mean": float(np.mean(overlaps)) if overlaps else float("nan"),
        f"top{top_k}_stable": bool(np.mean(overlaps) >= 0.99) if overlaps else False,
    }


# ------------------------------------------------------------------ failure taxonomy
FAILURE_CLASSES = (
    "acceptable",
    "chemically_invalid",
    "wrong_pocket",
    "pose_wrong_state_right",
    "state_wrong_pose_right",
    "membrane_mislocalized",
)


def failure_taxonomy(
    df: pd.DataFrame,
    config: MAESConfig,
    v_threshold: float = 0.5,
    m_threshold: float = 0.5,
    s_threshold: float = 0.5,
    rmsd_pose_bad: float = 2.0,
    rmsd_pose_good: float = 1.0,
) -> pd.DataFrame:
    """Assign one of five actionable failure classes to each pose (checklist F1)."""
    out = df.copy()
    rmsd = out["rmsd"].to_numpy(dtype=float) if "rmsd" in out.columns else np.full(len(out), np.nan)
    plif = (
        out["plif_recovery"].to_numpy(dtype=float)
        if "plif_recovery" in out.columns
        else np.full(len(out), np.nan)
    )
    V = out["V"].to_numpy(dtype=float)
    M = out["M"].to_numpy(dtype=float)
    S = out["S"].to_numpy(dtype=float)
    plif_ok = config.labels.eval_plif_threshold
    classes = []
    for i in range(len(out)):
        if V[i] < v_threshold:
            classes.append("chemically_invalid")
        elif np.isfinite(rmsd[i]) and rmsd[i] > rmsd_pose_bad and np.isfinite(plif[i]) and plif[i] < plif_ok:
            classes.append("wrong_pocket")
        elif np.isfinite(rmsd[i]) and rmsd[i] > rmsd_pose_bad and S[i] >= s_threshold:
            classes.append("pose_wrong_state_right")
        elif np.isfinite(rmsd[i]) and rmsd[i] <= rmsd_pose_good and S[i] < s_threshold:
            classes.append("state_wrong_pose_right")
        elif M[i] < m_threshold:
            classes.append("membrane_mislocalized")
        else:
            classes.append("acceptable")
    out["failure_class"] = classes
    return out


def failure_count_matrix(
    df: pd.DataFrame,
    method_col: str = "method",
    class_col: str = "failure_class",
) -> pd.DataFrame:
    """Method x failure-class count matrix (checklist F2)."""
    if method_col not in df.columns or class_col not in df.columns:
        raise KeyError(f"need columns {method_col!r} and {class_col!r}")
    table = pd.crosstab(df[method_col], df[class_col])
    for name in FAILURE_CLASSES:
        if name not in table.columns:
            table[name] = 0
    return table[list(FAILURE_CLASSES)]


# ============================================================================
# W2 statistical rules (Lead update, 2026-11-16)
# ----------------------------------------------------------------------------
# R1. EF1% is reported ONLY for targets with n_active >= 150; smaller targets
#     report EF5% / BEDROC(alpha=20) / AUC instead.
# R2. Track B has exactly three PREREGISTERED focus comparisons; everything
#     else is descriptive + failure taxonomy.
# R3. AcrB poses (55/68) are not independent samples (same crystal form /
#     target clustering): use cluster bootstrap by target/crystal form and a
#     separate AcrB group term in IRT; report n_eff ~ 20-35 in the main text.
# ============================================================================

EF1_MIN_ACTIVES = 150
BEDROC_ALPHA = 20.0

PREREGISTERED_TRACKB_COMPARISONS = (
    {
        "id": "C1",
        "name": "Vina vs Boltz-2 acceptable-pose rate",
        "unit": "pose",
        "statistic": "difference in acceptable rate (MAES >= 0.5)",
        "resampling": "cluster bootstrap by crystal form / target",
    },
    {
        "id": "C2",
        "name": "holo vs predicted receptor",
        "unit": "ligand",
        "statistic": "difference in acceptable rate (MAES >= 0.5)",
        "resampling": "cluster bootstrap by target",
    },
    {
        "id": "C3",
        "name": "water vs membrane MD pose stability",
        "unit": "system",
        "statistic": "difference in pose RMSD drift over the trajectory",
        "resampling": "cluster bootstrap by system",
    },
)


def ef1_allowed(n_active: int) -> bool:
    """R1: EF1% is only reportable when the target has at least 150 actives."""
    return int(n_active) >= EF1_MIN_ACTIVES


def primary_enrichment_metric(n_active: int) -> str:
    """Name of the headline enrichment metric under R1."""
    return "EF1%" if ef1_allowed(n_active) else "EF5%"


def enrichment_factor(y_true: Sequence[int], y_score: Sequence[float], fraction: float) -> float:
    """Enrichment factor at the top ``fraction`` of the ranked list."""
    if not 0.0 < float(fraction) <= 1.0:
        raise ValueError("fraction must be in (0, 1]")
    y = np.asarray(y_true, dtype=float)
    s = np.asarray(y_score, dtype=float)
    mask = np.isfinite(s)
    y, s = y[mask], s[mask]
    n = y.size
    n_pos = float(y.sum())
    if n == 0 or n_pos == 0:
        return float("nan")
    k = max(1, int(round(float(fraction) * n)))
    top = np.argsort(-s, kind="mergesort")[:k]
    hits = float(y[top].sum())
    expected = k * n_pos / n
    return float(hits / expected) if expected > 0 else float("nan")


def bedroc(y_true: Sequence[int], y_score: Sequence[float], alpha: float = BEDROC_ALPHA) -> float:
    """BEDROC (Truchon & Bayly 2004) with early-recognition weight ``alpha``."""
    y = np.asarray(y_true, dtype=float)
    s = np.asarray(y_score, dtype=float)
    mask = np.isfinite(s)
    y, s = y[mask], s[mask]
    n = y.size
    n_pos = int(round(float(y.sum())))
    if n == 0 or n_pos == 0 or n_pos == n:
        return float("nan")
    order = np.argsort(-s, kind="mergesort")
    y_ord = y[order]
    ranks_pos = np.flatnonzero(y_ord == 1) + 1.0
    weight = alpha / (1.0 - np.exp(-alpha))
    rie = weight * float(np.sum(np.exp(-alpha * ranks_pos / n))) / n_pos
    all_ranks = np.arange(1.0, n + 1.0)
    rie_max = weight * float(np.sum(np.exp(-alpha * all_ranks[:n_pos] / n))) / n_pos
    rie_min = weight * float(np.sum(np.exp(-alpha * all_ranks[-n_pos:] / n))) / n_pos
    if abs(rie_max - rie_min) < 1e-12:
        return float("nan")
    return float((rie - rie_min) / (rie_max - rie_min))


def enrichment_metrics(
    y_true: Sequence[int],
    y_score: Sequence[float],
    n_active: Optional[int] = None,
    alpha: float = BEDROC_ALPHA,
) -> Dict[str, Any]:
    """R1-compliant enrichment summary for one target.

    ``ef1`` is ``None`` (and ``ef1_suppressed=True``) whenever ``n_active < 150``.
    """
    y = np.asarray(y_true, dtype=int)
    n_active = int(y.sum()) if n_active is None else int(n_active)
    out: Dict[str, Any] = {
        "n_active": n_active,
        "ef1_allowed": ef1_allowed(n_active),
        "primary_metric": primary_enrichment_metric(n_active),
        "ef1": None,
        "ef1_suppressed": not ef1_allowed(n_active),
        "ef5": enrichment_factor(y, y_score, 0.05),
        "bedroc": bedroc(y, y_score, alpha=alpha),
        "auroc": auroc(y, y_score),
    }
    if ef1_allowed(n_active):
        out["ef1"] = enrichment_factor(y, y_score, 0.01)
    return out


# ----------------------------------------------------------- R3: clustered stats
def effective_sample_size(df: pd.DataFrame, cluster_col: str) -> Dict[str, Any]:
    """Report n_obs / n_clusters / design diagnostics (the number to put in text).

    For strongly clustered designs (AcrB: 34 of 40 primary poses), report the
    cluster count, the Kish effective sample size ``(sum m)^2 / sum m^2`` and —
    when an outcome column is available — the ICC-corrected n_eff from
    :func:`design_effect_n_eff`.
    """
    if cluster_col not in df.columns:
        raise KeyError(f"missing cluster column {cluster_col!r}")
    sizes = df[cluster_col].astype(str).value_counts()
    kish = float((sizes.sum() ** 2) / np.square(sizes.to_numpy(dtype=float)).sum())
    return {
        "cluster_col": cluster_col,
        "n_obs": int(len(df)),
        "n_clusters": int(sizes.size),
        "mean_cluster_size": float(sizes.mean()),
        "max_cluster_size": int(sizes.max()),
        "kish_n_eff": kish,
        "n_eff_range": [int(sizes.size), int(len(df))],
    }


def design_effect_n_eff(
    outcome: Sequence[float],
    clusters: Sequence[Any],
    level: float = 0.95,
) -> Dict[str, Any]:
    """ICC-corrected effective sample size for a clustered outcome.

    ``n_eff = n / (1 + (m_bar - 1) * ICC)`` with the one-way ANOVA estimate of
    the intraclass correlation (Kish design effect).  Use this with the MAES
    acceptable indicator on the inhibitor subset to quote a defensible n_eff.
    """
    y = np.asarray(outcome, dtype=float)
    c = np.asarray([str(v) for v in clusters], dtype=object)
    mask = np.isfinite(y)
    y, c = y[mask], c[mask]
    uniq, inverse = np.unique(c, return_inverse=True)
    k = uniq.size
    n = y.size
    if k < 2 or n <= k:
        return {"n_eff": float(n), "icc": 0.0, "n_clusters": int(k), "n_obs": int(n), "design_effect": 1.0}
    sizes = np.bincount(inverse).astype(float)
    m_bar = float(n / k)
    grand = float(y.mean())
    means = np.array([y[inverse == j].mean() for j in range(k)])
    ss_between = float(np.sum(sizes * (means - grand) ** 2))
    ss_within = float(np.sum([np.sum((y[inverse == j] - means[j]) ** 2) for j in range(k)]))
    ms_between = ss_between / (k - 1)
    ms_within = ss_within / (n - k)
    if ms_between + (m_bar - 1.0) * ms_within <= 0:
        icc = 0.0
    else:
        icc = float(max(0.0, (ms_between - ms_within) / (ms_between + (m_bar - 1.0) * ms_within)))
    design_effect = float(1.0 + (m_bar - 1.0) * icc)
    return {
        "n_eff": float(n / design_effect),
        "icc": icc,
        "design_effect": design_effect,
        "n_clusters": int(k),
        "n_obs": int(n),
        "mean_cluster_size": m_bar,
        "level": float(level),
    }


def n_eff_report(
    df: pd.DataFrame,
    cluster_col: str,
    outcome_col: Optional[str] = None,
) -> Dict[str, Any]:
    """Static (cluster count / Kish) + outcome-based (ICC) effective sample size."""
    report = effective_sample_size(df, cluster_col)
    if outcome_col is not None and outcome_col in df.columns:
        report["outcome_based"] = design_effect_n_eff(
            df[outcome_col].to_numpy(dtype=float), df[cluster_col].to_numpy()
        )
    return report


def _cluster_indices(clusters: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    unique = np.unique(clusters)
    sampled = rng.choice(unique, size=unique.size, replace=True)
    return np.concatenate([np.flatnonzero(clusters == c) for c in sampled])


def cluster_bootstrap_ci(
    y_true: Sequence[int],
    y_score: Sequence[float],
    clusters: Sequence[Any],
    metric=auroc,
    n_boot: int = 2000,
    seed: int = 42,
    level: float = 0.95,
) -> Dict[str, Any]:
    """Cluster bootstrap CI (resampling unit = target / crystal form, not pose)."""
    y = np.asarray(y_true, dtype=int)
    s = np.asarray(y_score, dtype=float)
    c = np.asarray([str(v) for v in clusters], dtype=object)
    if not (y.size == s.size == c.size):
        raise ValueError("y_true, y_score and clusters must have the same length")
    rng = np.random.default_rng(seed)
    stats = []
    for _ in range(int(n_boot)):
        idx = _cluster_indices(c, rng)
        try:
            value = metric(y[idx], s[idx])
        except Exception:
            continue
        if np.isfinite(value):
            stats.append(float(value))
    arr = np.asarray(stats, dtype=float)
    lo_q = (1.0 - level) / 2.0 * 100.0
    hi_q = (1.0 - (1.0 - level) / 2.0) * 100.0
    return {
        "estimate": float(metric(y, s)),
        "lo": float(np.percentile(arr, lo_q)) if arr.size else float("nan"),
        "hi": float(np.percentile(arr, hi_q)) if arr.size else float("nan"),
        "n_boot_valid": int(arr.size),
        "n_clusters": int(np.unique(c).size),
        "n_obs": int(y.size),
    }


def cluster_bootstrap_mean(
    values: Sequence[float],
    clusters: Sequence[Any],
    n_boot: int = 2000,
    seed: int = 42,
    level: float = 0.95,
) -> Dict[str, Any]:
    """Cluster bootstrap CI for a mean / paired difference (e.g. acceptance-rate gap)."""
    v = np.asarray(values, dtype=float)
    c = np.asarray([str(x) for x in clusters], dtype=object)
    mask = np.isfinite(v)
    v, c = v[mask], c[mask]
    if v.size == 0:
        return {"estimate": float("nan"), "lo": float("nan"), "hi": float("nan"), "n_clusters": 0}
    rng = np.random.default_rng(seed)
    stats = []
    for _ in range(int(n_boot)):
        idx = _cluster_indices(c, rng)
        stats.append(float(v[idx].mean()))
    arr = np.asarray(stats, dtype=float)
    lo_q = (1.0 - level) / 2.0 * 100.0
    return {
        "estimate": float(v.mean()),
        "lo": float(np.percentile(arr, lo_q)),
        "hi": float(np.percentile(arr, 100.0 - lo_q)),
        "n_boot_valid": int(arr.size),
        "n_clusters": int(np.unique(c).size),
        "n_obs": int(v.size),
    }


def wild_cluster_bootstrap_auc(
    y_true: Sequence[int],
    y_score: Sequence[float],
    clusters: Sequence[Any],
    n_boot: int = 4999,
    seed: int = 42,
) -> Dict[str, Any]:
    """Wild cluster bootstrap CI for AUROC (DeLong influence + Rademacher weights).

    Robustness check for the small-cluster regime (15 clusters in v1): ordinary
    cluster/ligand bootstrap CIs can be too narrow; this adds cluster-level
    Rademacher perturbations to the DeLong influence function.
    """
    y = np.asarray(y_true, dtype=int)
    s = np.asarray(y_score, dtype=float)
    c = np.asarray([str(value) for value in clusters])
    valid = np.isfinite(s) & np.isfinite(y.astype(float))
    y, s, c = y[valid], s[valid], c[valid]
    if len(y) == 0 or y.sum() < 2 or (len(y) - y.sum()) < 2:
        return {"auc": float("nan"), "ci_lo": float("nan"), "ci_hi": float("nan"),
                "n_boot": 0, "n_clusters": int(len(np.unique(c))),
                "method": "wild cluster bootstrap (DeLong influence, Rademacher)"}
    auc_value = float(auroc(y, s))
    pos = s[y == 1]
    neg = s[y == 0]
    n_pos, n_neg = len(pos), len(neg)
    v10 = np.array([
        (float(np.sum(p > neg)) + 0.5 * float(np.sum(p == neg))) / n_neg for p in pos
    ])
    v01 = np.array([
        (float(np.sum(pos > n)) + 0.5 * float(np.sum(pos == n))) / n_pos for n in neg
    ])
    influence = np.zeros(len(y), dtype=float)
    influence[y == 1] = (v10 - auc_value) / n_pos
    influence[y == 0] = (auc_value - v01) / n_neg
    cluster_sums = pd.Series(influence).groupby(pd.Series(c)).sum()
    rng = np.random.default_rng(int(seed))
    draws = np.empty(int(n_boot), dtype=float)
    values = cluster_sums.to_numpy(dtype=float)
    for idx in range(int(n_boot)):
        signs = rng.choice((-1.0, 1.0), size=len(values))
        draws[idx] = auc_value + float(np.sum(signs * values))
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return {
        "auc": auc_value,
        "ci_lo": float(lo),
        "ci_hi": float(hi),
        "n_boot": int(n_boot),
        "n_clusters": int(cluster_sums.size),
        "method": "wild cluster bootstrap (DeLong influence, Rademacher)",
    }


def acceptable_rate(df: pd.DataFrame, maes_col: str = "MAES", threshold: float = 0.5) -> float:
    """Fraction of poses whose MAES exceeds the acceptable threshold."""
    if maes_col not in df.columns:
        raise KeyError(f"missing {maes_col!r} column")
    return float((df[maes_col].to_numpy(dtype=float) >= float(threshold)).mean())


def paired_acceptance_difference(
    df: pd.DataFrame,
    method_a: str,
    method_b: str,
    method_col: str = "method",
    maes_col: str = "MAES",
    cluster_col: str = "crystal_form",
    pair_col: str = "pair_id",
    threshold: float = 0.5,
    n_boot: int = 2000,
    seed: int = 42,
) -> Dict[str, Any]:
    """C1-style paired difference in acceptable rate with a cluster-bootstrap CI.

    Expects one row per (pair, method); ``pair_id`` identifies the paired pose.
    """
    needed = {method_col, maes_col, cluster_col, pair_col}
    missing = needed - set(df.columns)
    if missing:
        raise KeyError(f"missing columns: {sorted(missing)}")
    sub = df[df[method_col].astype(str).isin([str(method_a), str(method_b)])].copy()
    sub["_acceptable"] = (sub[maes_col].to_numpy(dtype=float) >= float(threshold)).astype(float)
    pivot = sub.pivot_table(index=[pair_col, cluster_col], columns=method_col, values="_acceptable")
    if method_a not in pivot.columns or method_b not in pivot.columns:
        raise ValueError("both methods must be present for every pair")
    diff = (pivot[method_a] - pivot[method_b]).dropna()
    clusters = [idx[1] for idx in diff.index]
    stats = cluster_bootstrap_mean(diff.to_numpy(dtype=float), clusters, n_boot=n_boot, seed=seed)
    stats.update(
        {
            "comparison_id": "C1",
            "method_a": method_a,
            "method_b": method_b,
            "threshold": float(threshold),
            "description": "paired difference in acceptable rate (a - b)",
        }
    )
    return stats


def compare_arms(
    df: pd.DataFrame,
    arm_col: str,
    value_col: str,
    cluster_col: str,
    n_boot: int = 2000,
    seed: int = 42,
) -> Dict[str, Any]:
    """Generic C2/C3-style paired arm comparison (holo vs predicted; water vs membrane)."""
    needed = {arm_col, value_col, cluster_col}
    missing = needed - set(df.columns)
    if missing:
        raise KeyError(f"missing columns: {sorted(missing)}")
    pivot = df.pivot_table(index=cluster_col, columns=arm_col, values=value_col, aggfunc="mean")
    if pivot.shape[1] != 2:
        raise ValueError(f"expected exactly 2 arms in {arm_col!r}, got {pivot.shape[1]}")
    arms = list(pivot.columns)
    diff = (pivot[arms[0]] - pivot[arms[1]]).dropna()
    clusters = [str(x) for x in diff.index]
    stats = cluster_bootstrap_mean(diff.to_numpy(dtype=float), clusters, n_boot=n_boot, seed=seed)
    stats.update({"arm_a": arms[0], "arm_b": arms[1], "value_col": value_col, "cluster_col": cluster_col})
    return stats


# ============================================================================
# Metal-dependence stratification (Lead 2026-11-16)
# ----------------------------------------------------------------------------
# G1 pilot found a Vina enrichment INVERSION on the Zn-dependent target LpxC
# (raw AUC 0.212, EF1/EF5 = 0, active ranks 82-208) with the catalytic Zn
# stripped; Vina has no metal-coordination term.  Pilot v2 then showed that a
# Zn-RETAINED receptor does NOT rescue it (AUC 0.220 vs 0.222), and the DHFR
# non-metal control (n=213, same protocol) gives AUC 0.706, EF5 = 2.0.  This
# REPRODUCES a KNOWN limitation of physics-based scoring on metalloenzymes
# (literature: metalloenzyme DUE evaluation, metal-aware scoring reviews,
# Vina metal-as-H-bond-donor artifact) as a construct-validity/sanity control —
# it is NOT a new finding.  `metal_dependent` is a *target-level metadata*
# covariate for stratified reporting; it does NOT enter the MAES formula.
# ============================================================================

TARGET_METAL_DEPENDENCE = {
    # Track A (potency-driven targets)
    "LpxC_Ec": "Zn",
    "LpxC_Pa": "Zn",
    "DHFR_Ec": "none",
    "GyrA_Ec": "Mg",
    "Gyrase_Ec": "Mg",
    "GyrB_Ec": "Mg",
    "FabI_Ec": "none",
    "LepB_Ec": "none",
    "TEM1_Ec": "none",
    # Track B (structure-rich membrane targets; metadata only)
    "MsbA_Ec": "none",
    "AcrB_Ec": "none",
    "BamA_Ec": "none",
    "LptD_Ec": "none",
    "LptB_Ec": "Mg",
    "MexB_Pa": "none",
    "LptD_Pa": "none",
    "BamA_Pa": "none",
    "PBP3_Pa": "none",
}
METAL_STRATA = ("metal", "non-metal", "unknown")


def metal_stratum(metal_dependent) -> str:
    """Map a ``metal_dependent`` value to the reporting stratum."""
    if metal_dependent is None:
        return "unknown"
    if isinstance(metal_dependent, float) and np.isnan(metal_dependent):
        return "unknown"
    text = str(metal_dependent).strip().lower()
    if text in {"", "nan", "none"}:
        return "non-metal"
    if text in {"no", "false", "non-metal", "nonmetal", "non_metal"}:
        return "non-metal"
    return "metal"  # Zn / Mg / Fe / other metal


def add_metal_dependence(
    df: pd.DataFrame,
    target_col: str = "target_id",
    target_map: Optional[Mapping[str, str]] = None,
) -> pd.DataFrame:
    """Attach ``metal_dependent_raw`` and ``metal_stratum`` to a results table."""
    if target_col not in df.columns:
        raise KeyError(f"missing target column {target_col!r}")
    mapping = dict(target_map) if target_map is not None else dict(TARGET_METAL_DEPENDENCE)
    targets = df[target_col].astype(str).str.strip()
    raw = targets.map(lambda t: mapping.get(t, "UNKNOWN"))
    out = df.copy()
    out["metal_dependent_raw"] = raw
    out["metal_stratum"] = [metal_stratum(None if v == "UNKNOWN" else v) for v in raw]
    return out


def stratified_metric_table(
    df: pd.DataFrame,
    stratum_col: str = "metal_stratum",
    metric_cols: Optional[Sequence[str]] = None,
    strata: Sequence[str] = METAL_STRATA,
) -> pd.DataFrame:
    """Mean / median / n per stratum for the requested (or all numeric) metrics.

    Track A main tables report enrichment metrics (EF1% gated by n_active >= 150,
    else EF5%/BEDROC/AUC) split by ``metal_stratum``; this helper keeps the
    split consistent across targets and figures.
    """
    if stratum_col not in df.columns:
        raise KeyError(f"missing stratum column {stratum_col!r}")
    if metric_cols is None:
        metric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    cols = [c for c in metric_cols if c in df.columns]
    if not cols:
        raise ValueError("no numeric metric columns found")
    subset = df[df[stratum_col].isin(list(strata))]
    table = subset.groupby(stratum_col)[cols].agg(["mean", "median", "count"])
    return table
