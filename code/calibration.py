"""Calibration of the raw MAES score into ``P(pose acceptable)``.

Implements the v1.0 protocol from ``notes/maes_design.md`` §3.3:

    raw = V ** w_V * M ** w_M * S ** w_S          (weighted geometric mean)
    MAES = sigmoid(a + b * logit(raw))            (Platt)  or  isotonic(raw)

**Label separation is enforced at this layer** (red-team A2):

* the calibrator is fitted on ``build_calibration_labels`` (RMSD<2 A AND PoseBusters)
* the model is then evaluated on ``build_evaluation_labels`` (RMSD<1 A / PLIF)
* :class:`MAESModel` stores SHA-256 fingerprints of both label vectors, so the
  "same label in disguise" attack can be checked mechanically.

Implemented without sklearn for isotonic (PAVA fallback) so the package works
on a bare numpy/scipy environment; sklearn is used when available.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from .combine import combine_subscores
from .config import LabelConfig, LabelSeparationError, MAESConfig, default_config

#: isotonic calibration is enabled only above this many calibration positives
#: (small-n audit M3: n=36 / 8 positives -> Platt only)
ISOTONIC_MIN_POSITIVES = 100

__all__ = [
    "ISOTONIC_MIN_POSITIVES",
    "IsotonicCalibrator",
    "MAESModel",
    "PlattCalibrator",
    "assert_labels_independent",
    "build_calibration_labels",
    "build_evaluation_labels",
    "label_fingerprint",
    "make_calibrator",
]


def _sigmoid(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    exp_x = np.exp(x[~pos])
    out[~pos] = exp_x / (1.0 + exp_x)
    return out


def _logit(p: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    p = np.clip(np.asarray(p, dtype=float), eps, 1.0 - eps)
    return np.log(p / (1.0 - p))


def label_fingerprint(y: Sequence[int]) -> str:
    """SHA-256 of a label vector (allows a mechanical independence check)."""
    arr = np.asarray(y, dtype=np.int8)
    return hashlib.sha256(arr.tobytes()).hexdigest()


# --------------------------------------------------------------- label builders
def build_calibration_labels(df: pd.DataFrame, labels: LabelConfig) -> np.ndarray:
    """Calibration label: ``RMSD < calib_rmsd`` AND (optionally) PoseBusters pass."""
    if "rmsd" not in df.columns:
        raise KeyError("calibration labels need an 'rmsd' column")
    rmsd = pd.to_numeric(df["rmsd"], errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(rmsd) & (rmsd < float(labels.calib_rmsd))
    if labels.calib_require_posebusters:
        if "posebusters_pass" in df.columns:
            pb = df["posebusters_pass"].astype("boolean").fillna(False).astype(bool).to_numpy()
        elif "pb_pass_rate" in df.columns:
            pb = pd.to_numeric(df["pb_pass_rate"], errors="coerce").fillna(0.0).to_numpy(dtype=float) >= 0.999
        else:
            raise KeyError(
                "calib_require_posebusters=True needs a 'posebusters_pass' or 'pb_pass_rate' column"
            )
        ok = ok & pb
    return ok.astype(int)


def build_evaluation_labels(df: pd.DataFrame, labels: LabelConfig, which: str = "primary") -> np.ndarray:
    """Evaluation labels (never used for fitting).

    ``which='primary'`` uses the stricter RMSD cut; ``which='secondary'`` uses
    PLIF / contact recovery.
    """
    family = labels.eval_primary if which == "primary" else labels.eval_secondary
    if family == "rmsd_lt":
        if "rmsd" not in df.columns:
            raise KeyError("rmsd label needs a 'rmsd' column")
        rmsd = pd.to_numeric(df["rmsd"], errors="coerce").to_numpy(dtype=float)
        return (np.isfinite(rmsd) & (rmsd < float(labels.eval_rmsd))).astype(int)
    if family == "plif_recovery":
        if "plif_recovery" not in df.columns:
            raise KeyError("PLIF label needs a 'plif_recovery' column")
        plif = pd.to_numeric(df["plif_recovery"], errors="coerce").to_numpy(dtype=float)
        return (np.isfinite(plif) & (plif >= float(labels.eval_plif_threshold))).astype(int)
    raise LabelSeparationError(f"unknown evaluation label family: {family!r}")


def assert_labels_independent(calib: np.ndarray, evaluation: np.ndarray) -> None:
    """Refuse to proceed when calibration and evaluation labels are identical."""
    if np.asarray(calib).shape != np.asarray(evaluation).shape:
        raise LabelSeparationError("calibration and evaluation labels differ in length")
    if label_fingerprint(calib) == label_fingerprint(evaluation):
        raise LabelSeparationError(
            "calibration and evaluation label vectors are identical — MAES would be "
            "circular (red-team A2). Build the evaluation label from a stricter RMSD "
            "or from PLIF/contact recovery."
        )


# ------------------------------------------------------------------- calibrators
@dataclass
class PlattCalibrator:
    """1-D logistic calibration ``P = sigmoid(a + b * logit(raw))``."""

    l2: float = 1e-3
    a: float = 0.0
    b: float = 1.0
    fitted: bool = False

    def fit(self, raw: np.ndarray, y: np.ndarray) -> "PlattCalibrator":
        from scipy.optimize import minimize

        x = _logit(np.asarray(raw, dtype=float))
        y = np.asarray(y, dtype=float)
        if x.size == 0:
            raise ValueError("cannot fit a calibrator on an empty set")
        if x.size != y.size:
            raise ValueError("raw and y differ in length")

        def loss(params: np.ndarray) -> float:
            a, b = float(params[0]), float(params[1])
            z = a + b * x
            # numerically stable logistic log-likelihood
            ll = np.sum(y * -np.logaddexp(0.0, -z) + (1.0 - y) * -np.logaddexp(0.0, z))
            penalty = float(self.l2) * (a * a + (b - 1.0) ** 2)
            monotone_penalty = 1e3 * max(0.0, -b) ** 2
            return float(-ll + penalty + monotone_penalty)

        res = minimize(loss, np.array([0.0, 1.0]), method="L-BFGS-B")
        self.a, self.b = float(res.x[0]), float(res.x[1])
        if self.b <= 1e-6:
            # numerical guard: MAES must stay non-decreasing in the raw score
            def loss_a(a_only: np.ndarray) -> float:
                z = float(a_only[0]) + 1.0 * x
                ll = np.sum(y * -np.logaddexp(0.0, -z) + (1.0 - y) * -np.logaddexp(0.0, z))
                return float(-ll)

            res_a = minimize(loss_a, np.array([0.0]), method="L-BFGS-B")
            self.a, self.b = float(res_a.x[0]), 1.0
        self.fitted = True
        return self

    def predict(self, raw: np.ndarray) -> np.ndarray:
        if not self.fitted:
            raise RuntimeError("PlattCalibrator is not fitted")
        return np.clip(_sigmoid(self.a + self.b * _logit(np.asarray(raw, dtype=float))), 0.0, 1.0)

    def to_dict(self) -> Dict[str, Any]:
        return {"kind": "platt", "a": self.a, "b": self.b, "l2": self.l2, "fitted": self.fitted}

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "PlattCalibrator":
        return cls(l2=float(data.get("l2", 1e-3)), a=float(data["a"]), b=float(data["b"]), fitted=True)


def _pava(values: np.ndarray, weights: Optional[np.ndarray] = None) -> np.ndarray:
    """Pool-adjacent-violators (weighted, non-decreasing) for isotonic regression."""
    y = np.asarray(values, dtype=float)
    w = np.ones_like(y) if weights is None else np.asarray(weights, dtype=float)
    if y.size == 0:
        return y
    # blocks: [sum_wy, sum_w, count, start_idx]
    blocks: list = []
    for i, (yi, wi) in enumerate(zip(y, w)):
        blocks.append([yi * wi, wi, 1, i])
        while len(blocks) >= 2 and (blocks[-2][0] / blocks[-2][1]) > (blocks[-1][0] / blocks[-1][1]):
            wy2, w2, n2, s2 = blocks.pop()
            wy1, w1, n1, s1 = blocks.pop()
            blocks.append([wy1 + wy2, w1 + w2, n1 + n2, s1])
    fitted = np.empty_like(y)
    for wy, wsum, count, start in blocks:
        fitted[start:start + count] = wy / wsum
    return fitted


@dataclass
class IsotonicCalibrator:
    """Monotone piecewise-constant calibration (isotonic regression)."""

    x_: Optional[np.ndarray] = None
    y_: Optional[np.ndarray] = None
    fitted: bool = False

    def fit(self, raw: np.ndarray, y: np.ndarray) -> "IsotonicCalibrator":
        raw = np.asarray(raw, dtype=float)
        y = np.asarray(y, dtype=float)
        if raw.size != y.size:
            raise ValueError("raw and y differ in length")
        try:  # prefer sklearn when available
            from sklearn.isotonic import IsotonicRegression

            iso = IsotonicRegression(out_of_bounds="clip", increasing=True)
            iso.fit(raw, y)
            grid = np.unique(raw)
            self.x_ = grid
            self.y_ = np.clip(iso.predict(grid), 0.0, 1.0)
        except Exception:
            order = np.argsort(raw)
            xs = raw[order]
            ys = y[order]
            uniq, inverse = np.unique(xs, return_inverse=True)
            sums = np.zeros(uniq.size)
            counts = np.zeros(uniq.size)
            np.add.at(sums, inverse, ys)
            np.add.at(counts, inverse, 1.0)
            fitted = _pava(sums / counts, counts)
            self.x_ = uniq
            self.y_ = np.clip(fitted, 0.0, 1.0)
        self.fitted = True
        return self

    def predict(self, raw: np.ndarray) -> np.ndarray:
        if not self.fitted:
            raise RuntimeError("IsotonicCalibrator is not fitted")
        values = np.interp(
            np.asarray(raw, dtype=float), self.x_, self.y_, left=float(self.y_[0]), right=float(self.y_[-1])
        )
        return np.clip(values, 0.0, 1.0)

    def to_dict(self) -> Dict[str, Any]:
        if not self.fitted:
            raise RuntimeError("IsotonicCalibrator is not fitted")
        return {"kind": "isotonic", "x": self.x_.tolist(), "y": self.y_.tolist(), "fitted": True}

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "IsotonicCalibrator":
        return cls(
            x_=np.asarray(data["x"], dtype=float),
            y_=np.asarray(data["y"], dtype=float),
            fitted=True,
        )


class _IdentityCalibrator:
    """No calibration (ablation only): MAES == clipped raw score."""

    kind = "none"

    def fit(self, raw: np.ndarray, y: np.ndarray) -> "_IdentityCalibrator":
        return self

    def predict(self, raw: np.ndarray) -> np.ndarray:
        return np.clip(np.asarray(raw, dtype=float), 0.0, 1.0)

    def to_dict(self) -> Dict[str, Any]:
        return {"kind": "none"}

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "_IdentityCalibrator":
        return cls()


def make_calibrator(kind: str = "platt", l2: float = 1e-3):
    kind = str(kind).lower()
    if kind == "platt":
        return PlattCalibrator(l2=l2)
    if kind == "isotonic":
        return IsotonicCalibrator()
    if kind == "none":
        return _IdentityCalibrator()
    raise ValueError(f"unknown calibrator kind: {kind!r}")


# ------------------------------------------------------------------- MAES model
@dataclass
class MAESModel:
    """Fitted MAES scorer: sub-score weights + calibrator + frozen metadata."""

    config: MAESConfig = field(default_factory=default_config)
    weights: Dict[str, float] = field(default_factory=dict)
    calibrator: Any = None
    fitted: bool = False
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.config.validate()
        if not self.weights:
            self.weights = dict(self.config.scores.weights)

    # ------------------------------------------------------------------ training
    def raw_scores(self, df: pd.DataFrame) -> pd.Series:
        return combine_subscores(
            df,
            kind=self.config.scores.combine,
            weights=self.weights,
            p=self.config.scores.softmin_p,
            eps=self.config.scores.logit_eps,
        )

    def fit(
        self,
        df: pd.DataFrame,
        y: Optional[Sequence[int]] = None,
        label_column: Optional[str] = None,
    ) -> "MAESModel":
        """Fit the calibrator on the **calibration** label only.

        Weights are fixed a-priori config weights (no data-derived weighting in
        the fit); Platt is the primary calibrator and isotonic is gated to
        ``n_pos >= ISOTONIC_MIN_POSITIVES`` (small-n audit M3).
        """
        cfg = self.config
        if y is None and label_column is not None:
            y = df[label_column].to_numpy()
        calib = build_calibration_labels(df, cfg.labels) if y is None else np.asarray(y, dtype=int)
        eval_primary = build_evaluation_labels(df, cfg.labels, which="primary")
        eval_secondary = build_evaluation_labels(df, cfg.labels, which="secondary")
        assert_labels_independent(calib, eval_primary)
        # primary label independence is a hard gate; the secondary label may be
        # numerically identical on a small/synthetic set (then it carries no extra
        # information and is reported as unusable instead of aborting the fit)
        secondary_independent = True
        try:
            assert_labels_independent(calib, eval_secondary)
        except LabelSeparationError:
            secondary_independent = False

        raw = self.raw_scores(df).to_numpy()
        kind = cfg.calibrator.kind
        n_pos = int(calib.sum())
        n_neg = int(len(calib) - n_pos)
        fallback_reason = ""
        isotonic_gate_applied = False
        if kind == "isotonic" and n_pos < ISOTONIC_MIN_POSITIVES:
            fallback_reason = (
                f"isotonic requires >= {ISOTONIC_MIN_POSITIVES} calibration positives "
                f"(n_pos={n_pos}); falling back to Platt"
            )
            isotonic_gate_applied = True
            kind = "platt"
        if kind != "none" and (n_pos == 0 or n_neg == 0):
            fallback_reason = "single-class calibration set; falling back to identity calibrator"
            kind = "none"
            isotonic_gate_applied = False
        self.calibrator = make_calibrator(kind, l2=cfg.calibrator.l2).fit(raw, calib)
        self.fitted = True
        self.metadata = {
            "maes_version": cfg.version,
            "n_fit": int(len(df)),
            "n_calibration_positive": n_pos,
            "n_calibration_negative": n_neg,
            "calibration_positive_rate": float(n_pos / max(1, len(calib))),
            "calibrator": kind,
            "calibrator_fallback_reason": fallback_reason,
            "calibrator_policy": {
                "primary": "platt",
                "requested": cfg.calibrator.kind,
                "isotonic_min_positives": ISOTONIC_MIN_POSITIVES,
                "isotonic_gate_applied": bool(isotonic_gate_applied),
            },
            "weights": {k: float(v) for k, v in self.weights.items()},
            "weights_protocol": {
                "scheme": "fixed a-priori config weights (V=M=S=1/3 in v1.0)",
                "entropy_weights_used_in_fit": False,
                "frozen_at_inference": True,
            },
            "combine": cfg.scores.combine,
            "label_protocol": cfg.labels.independence_report(),
            "label_fingerprints": {
                "calibration": label_fingerprint(calib),
                "evaluation_primary": label_fingerprint(eval_primary),
                "evaluation_secondary": label_fingerprint(eval_secondary),
            },
            "evaluation_secondary_independent": bool(secondary_independent),
            "monotone": bool(
                isinstance(self.calibrator, PlattCalibrator) and self.calibrator.b > 0
            ) or isinstance(self.calibrator, IsotonicCalibrator),
        }
        return self

    # ----------------------------------------------------------------- inference
    def predict(self, df: pd.DataFrame) -> np.ndarray:
        if not self.fitted:
            raise RuntimeError("MAESModel is not fitted")
        raw = self.raw_scores(df).to_numpy()
        return np.clip(self.calibrator.predict(raw), 0.0, 1.0)

    def predict_frame(self, df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        raw = self.raw_scores(df)
        out["raw"] = raw.to_numpy()
        out["MAES"] = self.predict(df)
        return out

    # ------------------------------------------------------------------ persist
    def to_dict(self) -> Dict[str, Any]:
        if not self.fitted:
            raise RuntimeError("MAESModel is not fitted")
        return {
            "maes_version": self.config.version,
            "config": self.config.to_dict(),
            "weights": {k: float(v) for k, v in self.weights.items()},
            "calibrator": self.calibrator.to_dict(),
            "fitted": True,
            "metadata": self.metadata,
        }

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")
        return path

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "MAESModel":
        cfg = MAESConfig.from_dict(data["config"])
        kind = data["calibrator"]["kind"]
        calibrator = make_calibrator(kind).from_dict(data["calibrator"])
        model = cls(config=cfg, weights=dict(data.get("weights", {})), calibrator=calibrator)
        model.fitted = True
        model.metadata = dict(data.get("metadata", {}))
        return model

    @classmethod
    def load(cls, path: str | Path) -> "MAESModel":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))
