"""Sub-score combination rules and weighting utilities.

* geometric mean (design doc §3.3 method A, the v1.0 primary rule)
* soft-min (method B; penalises a single failing axis; sensitivity analysis)
* arithmetic mean (the naive baseline a reviewer will ask for)
* entropy weights (objective weights; design doc §4)
* Dirichlet random weights (robustness; reporting checklist A2)
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "SUBSCORES",
    "normalize_weights",
    "geometric_combine",
    "softmin_combine",
    "arithmetic_combine",
    "combine_subscores",
    "entropy_weights",
    "dirichlet_weights",
]

SUBSCORES = ("V", "M", "S")


def _subscore_matrix(df: pd.DataFrame, columns=SUBSCORES) -> np.ndarray:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise KeyError(f"missing sub-score columns: {missing}")
    arr = df[list(columns)].to_numpy(dtype=float)
    if not np.isfinite(arr).all():
        raise ValueError("sub-scores contain non-finite values")
    return arr


def normalize_weights(weights) -> np.ndarray:
    """Accept dict/list/array of 3 weights; return a normalised (3,) array."""
    if isinstance(weights, dict):
        vals = np.array([float(weights.get(k, 0.0)) for k in SUBSCORES])
    else:
        vals = np.asarray(weights, dtype=float).reshape(-1)
    if vals.size != 3:
        raise ValueError("expected 3 weights (V, M, S)")
    if (vals < 0).any():
        raise ValueError("weights must be non-negative")
    total = float(vals.sum())
    if total <= 0:
        raise ValueError("weights must not all be zero")
    return vals / total


def geometric_combine(df: pd.DataFrame, weights=None, eps: float = 1e-6) -> pd.Series:
    """Weighted geometric mean: ``prod(x_i ** w_i)`` with ``x`` clipped to [eps, 1]."""
    arr = np.clip(_subscore_matrix(df), eps, 1.0)
    w = normalize_weights(weights if weights is not None else {"V": 1 / 3, "M": 1 / 3, "S": 1 / 3})
    raw = np.exp(np.log(arr) @ w)
    return pd.Series(np.clip(raw, 0.0, 1.0), index=df.index, name="raw")


def softmin_combine(df: pd.DataFrame, weights=None, p: float = 8.0, eps: float = 1e-6) -> pd.Series:
    """Weighted soft-min: ``(sum(w_i * x_i ** -p)) ** (-1/p)`` (lower than the min side)."""
    if p <= 0:
        raise ValueError("p must be positive")
    arr = np.clip(_subscore_matrix(df), eps, 1.0)
    w = normalize_weights(weights if weights is not None else {"V": 1 / 3, "M": 1 / 3, "S": 1 / 3})
    raw = (np.sum(w.reshape(1, -1) * arr ** (-p), axis=1)) ** (-1.0 / p)
    return pd.Series(np.clip(raw, 0.0, 1.0), index=df.index, name="raw")


def arithmetic_combine(df: pd.DataFrame, weights=None) -> pd.Series:
    """Weighted arithmetic mean (naive baseline)."""
    arr = np.clip(_subscore_matrix(df), 0.0, 1.0)
    w = normalize_weights(weights if weights is not None else {"V": 1 / 3, "M": 1 / 3, "S": 1 / 3})
    raw = arr @ w
    return pd.Series(np.clip(raw, 0.0, 1.0), index=df.index, name="raw")


def combine_subscores(
    df: pd.DataFrame,
    kind: str = "geometric",
    weights=None,
    p: float = 8.0,
    eps: float = 1e-6,
) -> pd.Series:
    """Dispatch to the requested combination rule."""
    kind = str(kind).lower()
    if kind == "geometric":
        return geometric_combine(df, weights=weights, eps=eps)
    if kind == "softmin":
        return softmin_combine(df, weights=weights, p=p, eps=eps)
    if kind == "arithmetic":
        return arithmetic_combine(df, weights=weights)
    raise ValueError(f"unknown combination rule: {kind!r}")


def entropy_weights(df: pd.DataFrame, columns=SUBSCORES, eps: float = 1e-12) -> dict:
    """Entropy-weight method (objective weights, design doc §4).

    Columns are min-max scaled, converted to per-column probability mass, and
    weighted by ``1 - normalised entropy``.  Falls back to uniform weights when
    every column is constant.
    """
    arr = _subscore_matrix(df, columns)
    n = arr.shape[0]
    if n == 0:
        raise ValueError("cannot compute entropy weights on an empty table")
    lo = arr.min(axis=0)
    hi = arr.max(axis=0)
    span = hi - lo
    scaled = np.where(span > 0, (arr - lo) / np.where(span > 0, span, 1.0), 0.5)
    mass = scaled + eps
    mass = mass / mass.sum(axis=0, keepdims=True)
    entropy = -(mass * np.log(mass)).sum(axis=0) / np.log(n) if n > 1 else np.zeros(3)
    divergence = 1.0 - entropy
    if float(divergence.sum()) <= 1e-12:
        w = np.full(3, 1.0 / 3.0)
    else:
        w = divergence / divergence.sum()
    return {name: float(w[i]) for i, name in enumerate(columns)}


def dirichlet_weights(
    n_samples: int = 1000,
    n_components: int = 3,
    alpha: float = 1.0,
    seed: int = 42,
) -> np.ndarray:
    """Draw ``n_samples`` weight vectors from ``Dirichlet(alpha * ones(k))``."""
    if n_samples <= 0:
        raise ValueError("n_samples must be positive")
    if n_components <= 0 or alpha <= 0:
        raise ValueError("n_components and alpha must be positive")
    rng = np.random.default_rng(seed)
    return rng.dirichlet(np.full(n_components, float(alpha)), size=int(n_samples))
