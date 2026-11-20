"""IRT / hierarchical-Bayes skeleton for MAES v2 (design doc §3.3 method C).

Measurement model (psychometrics): poses are *items*, methods are *examinees*,
targets carry *item difficulty*.

    Rasch:  logit P(y_ij = 1) = theta_i - b_j
    2PL:    logit P(y_ij = 1) = a_i * (theta_i - b_j)

with ``theta_i`` the ability of method ``i`` and ``b_j`` the difficulty of
target ``j``.  A light L2 prior plays the role of pooling; :func:`fit_hierarchical_rasch`
additionally performs empirical-Bayes shrinkage (partial pooling surrogate) and
is the v2 skeleton to be replaced by full PyMC when the environment allows
(``pymc_available()`` is provided for that decision).

Python-only (numpy/scipy); no GPU, no internet.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

__all__ = [
    "IRTResult",
    "build_responses",
    "fit_2pl",
    "fit_hierarchical_rasch",
    "fit_rasch",
    "pymc_available",
]


def pymc_available() -> bool:
    try:  # pragma: no cover - depends on environment
        import pymc  # noqa: F401

        return True
    except Exception:
        return False


@dataclass
class IRTResult:
    methods: list
    targets: list
    theta: np.ndarray
    b: np.ndarray
    model: str
    n_obs: int
    loglik: float
    discrim: Optional[np.ndarray] = None
    theta_ci: Optional[np.ndarray] = None
    b_ci: Optional[np.ndarray] = None
    backend: str = "custom-optimizer"
    extra: Dict[str, Any] = field(default_factory=dict)

    # ------------------------------------------------------------------ reports
    def ability_frame(self) -> pd.DataFrame:
        data = {"method": self.methods, "theta": self.theta}
        if self.theta_ci is not None:
            data["theta_lo"], data["theta_hi"] = self.theta_ci[:, 0], self.theta_ci[:, 1]
        if self.discrim is not None:
            data["discrimination"] = self.discrim
        return pd.DataFrame(data).sort_values("theta", ascending=False).reset_index(drop=True)

    def difficulty_frame(self) -> pd.DataFrame:
        data = {"target": self.targets, "b": self.b}
        if self.b_ci is not None:
            data["b_lo"], data["b_hi"] = self.b_ci[:, 0], self.b_ci[:, 1]
        return pd.DataFrame(data).sort_values("b", ascending=False).reset_index(drop=True)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "model": self.model,
            "backend": self.backend,
            "n_obs": int(self.n_obs),
            "loglik": float(self.loglik),
            "methods": list(self.methods),
            "targets": list(self.targets),
            "theta": self.theta.tolist(),
            "b": self.b.tolist(),
            "discrim": None if self.discrim is None else self.discrim.tolist(),
            "theta_ci": None if self.theta_ci is None else self.theta_ci.tolist(),
            "b_ci": None if self.b_ci is None else self.b_ci.tolist(),
            "extra": dict(self.extra),
        }


def _build_design(
    responses: pd.DataFrame,
    method_col: str,
    target_col: str,
    response_col: str,
) -> tuple:
    for column in (method_col, target_col, response_col):
        if column not in responses.columns:
            raise KeyError(f"responses needs a {column!r} column")
    frame = responses[[method_col, target_col, response_col]].dropna()
    y = frame[response_col].to_numpy(dtype=float)
    if not set(np.unique(y)).issubset({0.0, 1.0}):
        raise ValueError("response column must be binary (0/1)")
    methods = sorted(frame[method_col].astype(str).unique().tolist())
    targets = sorted(frame[target_col].astype(str).unique().tolist())
    m_index = {m: i for i, m in enumerate(methods)}
    t_index = {t: j for j, t in enumerate(targets)}
    mi = frame[method_col].astype(str).map(m_index).to_numpy()
    ti = frame[target_col].astype(str).map(t_index).to_numpy()
    return y, mi, ti, methods, targets


def build_responses(
    df: pd.DataFrame,
    label: str = "calibration",
    method_col: str = "method",
    target_col: str = "system_id",
    eval_fn=None,
) -> pd.DataFrame:
    """Turn a pose table into a binary (method, target) response table.

    ``label='calibration'`` uses RMSD<2 A AND PoseBusters; ``label='evaluation'``
    uses the independent primary label.  ``eval_fn(df) -> np.ndarray`` overrides
    both (must return 0/1).
    """
    if eval_fn is not None:
        y = np.asarray(eval_fn(df), dtype=int)
    else:
        from .calibration import build_calibration_labels, build_evaluation_labels
        from .config import default_config

        cfg = default_config()
        if label == "calibration":
            y = build_calibration_labels(df, cfg.labels)
        elif label == "evaluation":
            y = build_evaluation_labels(df, cfg.labels, which="primary")
        else:
            raise ValueError(f"unknown label: {label!r}")
    if method_col not in df.columns or target_col not in df.columns:
        raise KeyError(f"need {method_col!r} and {target_col!r} columns")
    return pd.DataFrame(
        {
            "method": df[method_col].astype(str).to_numpy(),
            "target": df[target_col].astype(str).to_numpy(),
            "y": y,
        }
    )


def _loglik(y: np.ndarray, z: np.ndarray) -> float:
    return float(np.sum(y * -np.logaddexp(0.0, -z) + (1.0 - y) * -np.logaddexp(0.0, z)))


def _fit_rasch_core(y: np.ndarray, mi: np.ndarray, ti: np.ndarray, n_m: int, n_t: int, l2: float) -> tuple:
    from scipy.optimize import minimize

    def unpack(params: np.ndarray) -> tuple:
        return params[:n_m], params[n_m:n_m + n_t]

    def loss(params: np.ndarray) -> float:
        theta, b = unpack(params)
        z = theta[mi] - b[ti]
        return float(-_loglik(y, z) + l2 * (np.sum(theta ** 2) + np.sum(b ** 2)))

    x0 = np.zeros(n_m + n_t)
    res = minimize(loss, x0, method="L-BFGS-B")
    theta, b = unpack(res.x)
    # identifiability: shift both by a constant so that mean(theta) = 0
    shift = float(theta.mean())
    theta = theta - shift
    b = b - shift
    return theta, b, float(-loss(res.x))


def _bootstrap_cis(
    y: np.ndarray,
    mi: np.ndarray,
    ti: np.ndarray,
    n_m: int,
    n_t: int,
    l2: float,
    n_boot: int,
    seed: int,
    two_pl: bool = False,
) -> tuple:
    """Cluster bootstrap over targets (the correct resampling unit for difficulty)."""
    rng = np.random.default_rng(seed)
    theta_draws, b_draws = [], []
    unique_targets = np.unique(ti)
    if n_boot <= 0 or unique_targets.size < 2:
        return None, None
    for _ in range(int(n_boot)):
        sampled = rng.choice(unique_targets, size=unique_targets.size, replace=True)
        idx = np.concatenate([np.flatnonzero(ti == t) for t in sampled])
        if idx.size == 0:
            continue
        try:
            if two_pl:
                theta, b, disc, _ = _fit_2pl_core(y[idx], mi[idx], ti[idx], n_m, n_t, l2)
            else:
                theta, b, _ = _fit_rasch_core(y[idx], mi[idx], ti[idx], n_m, n_t, l2)
            theta_draws.append(theta)
            b_draws.append(b)
        except Exception:
            continue
    if not theta_draws:
        return None, None
    theta_arr = np.vstack(theta_draws)
    b_arr = np.vstack(b_draws)
    return (
        np.percentile(theta_arr, [2.5, 97.5], axis=0).T,
        np.percentile(b_arr, [2.5, 97.5], axis=0).T,
    )


def fit_rasch(
    responses: pd.DataFrame,
    method_col: str = "method",
    target_col: str = "target",
    response_col: str = "y",
    l2: float = 1.0,
    n_boot: int = 0,
    seed: int = 42,
) -> IRTResult:
    """Rasch model via L2-regularised maximum likelihood."""
    y, mi, ti, methods, targets = _build_design(responses, method_col, target_col, response_col)
    theta, b, loglik = _fit_rasch_core(y, mi, ti, len(methods), len(targets), l2=float(l2))
    theta_ci, b_ci = _bootstrap_cis(
        y, mi, ti, len(methods), len(targets), float(l2), n_boot, seed, two_pl=False
    )
    return IRTResult(
        methods=methods,
        targets=targets,
        theta=theta,
        b=b,
        model="rasch",
        n_obs=int(y.size),
        loglik=loglik,
        theta_ci=theta_ci,
        b_ci=b_ci,
    )


def _fit_2pl_core(y: np.ndarray, mi: np.ndarray, ti: np.ndarray, n_m: int, n_t: int, l2: float) -> tuple:
    from scipy.optimize import minimize

    def unpack(params: np.ndarray) -> tuple:
        theta = params[:n_m]
        logit_a = params[n_m:n_m + n_m]
        b = params[n_m + n_m:]
        disc = np.exp(logit_a)  # positive discrimination
        return theta, b, disc

    def loss(params: np.ndarray) -> float:
        theta, b, disc = unpack(params)
        z = disc[mi] * (theta[mi] - b[ti])
        return float(
            -_loglik(y, z)
            + l2 * (np.sum(theta ** 2) + np.sum(b ** 2) + np.sum((np.log(disc)) ** 2))
        )

    x0 = np.zeros(2 * n_m + n_t)
    res = minimize(loss, x0, method="L-BFGS-B")
    theta, b, disc = unpack(res.x)
    shift = float(theta.mean())
    theta = theta - shift
    b = b - shift
    return theta, b, disc, float(-loss(res.x))


def fit_2pl(
    responses: pd.DataFrame,
    method_col: str = "method",
    target_col: str = "target",
    response_col: str = "y",
    l2: float = 1.0,
    n_boot: int = 0,
    seed: int = 42,
) -> IRTResult:
    """2PL model (method-specific discrimination) via L2-regularised likelihood."""
    y, mi, ti, methods, targets = _build_design(responses, method_col, target_col, response_col)
    theta, b, disc, loglik = _fit_2pl_core(y, mi, ti, len(methods), len(targets), l2=float(l2))
    theta_ci, b_ci = _bootstrap_cis(
        y, mi, ti, len(methods), len(targets), float(l2), n_boot, seed, two_pl=True
    )
    return IRTResult(
        methods=methods,
        targets=targets,
        theta=theta,
        b=b,
        discrim=disc,
        model="2pl",
        n_obs=int(y.size),
        loglik=loglik,
        theta_ci=theta_ci,
        b_ci=b_ci,
    )


def fit_hierarchical_rasch(
    responses: pd.DataFrame,
    method_col: str = "method",
    target_col: str = "target",
    response_col: str = "y",
    l2: float = 1.0,
    kappa: float = 5.0,
    n_boot: int = 0,
    seed: int = 42,
) -> IRTResult:
    """Empirical-Bayes (partial pooling) Rasch skeleton for MAES v2.

    Ability/difficulty estimates from :func:`fit_rasch` are shrunk toward the
    population mean with a Dirichlet-style pseudo-count ``kappa``::

        b_j^EB = (n_j * b_j + kappa * mean(b)) / (n_j + kappa)

    which is the cheap surrogate for the hierarchical-bayes model in the design
    doc.  If PyMC is installed the caller can swap in a full model later; this
    function is deliberately dependency-free.
    """
    base = fit_rasch(
        responses, method_col=method_col, target_col=target_col,
        response_col=response_col, l2=l2, n_boot=n_boot, seed=seed,
    )
    frame = responses[[method_col, target_col, response_col]].dropna()
    n_target = frame.groupby(target_col)[response_col].size().to_dict()
    n_method = frame.groupby(method_col)[response_col].size().to_dict()
    mean_b = float(np.mean(base.b))
    b_eb = np.array(
        [
            (float(n_target.get(t, 0)) * base.b[j] + float(kappa) * mean_b)
            / (float(n_target.get(t, 0)) + float(kappa))
            for j, t in enumerate(base.targets)
        ]
    )
    theta_eb = np.array(
        [
            (float(n_method.get(m, 0)) * base.theta[i]) / (float(n_method.get(m, 0)) + float(kappa))
            for i, m in enumerate(base.methods)
        ]
    )
    return IRTResult(
        methods=base.methods,
        targets=base.targets,
        theta=theta_eb,
        b=b_eb,
        model="rasch-empirical-bayes",
        n_obs=base.n_obs,
        loglik=base.loglik,
        theta_ci=base.theta_ci,
        b_ci=base.b_ci,
        backend="custom-optimizer (pymc "
        + ("available → upgrade path open" if pymc_available() else "not installed → upgrade path documented")
        + ")",
        extra={"kappa": float(kappa)},
    )


# ============================================================================
# R3 (Lead update 2026-11-16): separate group term for clustered targets.
# AcrB contributes 55/68 poses but only a handful of independent crystal
# forms, so its poses share a group-level offset:
#
#     logit P(y_ij = 1) = theta_i - b_j + g_{group(j)}
#
# where g is a zero-mean group effect (e.g. per crystal form or per target
# family).  With one group per target this reduces to the plain Rasch model.
# ============================================================================
def _build_design_with_groups(
    responses: pd.DataFrame,
    method_col: str,
    target_col: str,
    response_col: str,
    group_col: str,
) -> tuple:
    for column in (method_col, target_col, response_col, group_col):
        if column not in responses.columns:
            raise KeyError(f"responses needs a {column!r} column")
    frame = responses[[method_col, target_col, response_col, group_col]].dropna().copy()
    y = frame[response_col].to_numpy(dtype=float)
    if not set(np.unique(y)).issubset({0.0, 1.0}):
        raise ValueError("response column must be binary (0/1)")
    methods = sorted(frame[method_col].astype(str).unique().tolist())
    targets = sorted(frame[target_col].astype(str).unique().tolist())
    groups = sorted(frame[group_col].astype(str).unique().tolist())
    m_index = {m: i for i, m in enumerate(methods)}
    t_index = {t: j for j, t in enumerate(targets)}
    g_index = {g: k for k, g in enumerate(groups)}
    mi = frame[method_col].astype(str).map(m_index).to_numpy()
    ti = frame[target_col].astype(str).map(t_index).to_numpy()
    # group effect is a per-row offset, so one target may span several crystal forms
    gi = frame[group_col].astype(str).map(g_index).to_numpy()
    return y, mi, ti, gi, groups, methods, targets


def fit_rasch_grouped(
    responses: pd.DataFrame,
    group_col: str = "crystal_form",
    method_col: str = "method",
    target_col: str = "target",
    response_col: str = "y",
    l2: float = 1.0,
    n_boot: int = 0,
    seed: int = 42,
) -> IRTResult:
    """Rasch model with an additive group effect ``g`` (R3 cluster correction)."""
    from scipy.optimize import minimize

    y, mi, ti, gi, groups, methods, targets = _build_design_with_groups(
        responses, method_col, target_col, response_col, group_col
    )
    n_m, n_t, n_g = len(methods), len(targets), len(groups)

    def unpack(params: np.ndarray) -> tuple:
        theta = params[:n_m]
        b = params[n_m:n_m + n_t]
        g = params[n_m + n_t:]
        return theta, b, g

    def loss(params: np.ndarray) -> float:
        theta, b, g = unpack(params)
        z = theta[mi] - b[ti] + g[gi]
        penalty = l2 * (np.sum(theta ** 2) + np.sum(b ** 2) + np.sum(g ** 2))
        return float(-_loglik(y, z) + penalty)

    x0 = np.zeros(n_m + n_t + n_g)
    res = minimize(loss, x0, method="L-BFGS-B")
    theta, b, g = unpack(res.x)
    # identifiability: centre theta and g (adding a constant to both cancels)
    shift = float(theta.mean())
    theta = theta - shift
    b = b - shift
    g = g - float(g.mean())
    z = theta[mi] - b[ti] + g[gi]
    return IRTResult(
        methods=methods,
        targets=targets,
        theta=theta,
        b=b,
        model="rasch-grouped",
        n_obs=int(y.size),
        loglik=float(_loglik(y, z)),
        extra={
            "group_col": group_col,
            "groups": groups,
            "n_groups": int(n_g),
            "group_effect": g.tolist(),
            "note": (
                "group effect is a per-row offset so AcrB crystal forms can share a "
                "target difficulty but keep separate cluster effects; report n_eff ~ 20-35"
            ),
        },
    )


def group_effect_frame(result: IRTResult) -> pd.DataFrame:
    """Group-level offsets from :func:`fit_rasch_grouped` (empty otherwise)."""
    groups = result.extra.get("groups", [])
    effects = result.extra.get("group_effect", [])
    if not groups or not effects:
        return pd.DataFrame(columns=["group", "g"])
    return pd.DataFrame({"group": groups, "g": effects}).sort_values("g", ascending=False)


def fit_hierarchical_rasch_grouped(
    responses: pd.DataFrame,
    group_col: str = "crystal_form",
    method_col: str = "method",
    target_col: str = "target",
    response_col: str = "y",
    l2: float = 1.0,
    kappa: float = 5.0,
    seed: int = 42,
) -> IRTResult:
    """Grouped Rasch + empirical-Bayes shrinkage of theta (v2 skeleton)."""
    base = fit_rasch_grouped(
        responses, group_col=group_col, method_col=method_col,
        target_col=target_col, response_col=response_col, l2=l2, seed=seed,
    )
    counts = responses.groupby(method_col)[response_col].size().to_dict()
    theta_eb = np.array(
        [
            (float(counts.get(m, 0)) * base.theta[i]) / (float(counts.get(m, 0)) + float(kappa))
            for i, m in enumerate(base.methods)
        ]
    )
    base.theta = theta_eb
    base.model = "rasch-grouped-empirical-bayes"
    base.extra["kappa"] = float(kappa)
    return base
