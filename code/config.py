"""Configuration objects for MAES v1.0.

Two ideas are load-bearing here:

1. **Label separation** (red-team item A2, the most fatal attack on the paper).
   The calibration label and the evaluation label are *different measurements*
   and live in different fields.  :class:`LabelConfig.validate` refuses to run
   if they coincide, unless the caller explicitly opts into ``allow_shared_threshold``
   (which then shows up as ``strictly_independent = False`` in the model
   metadata and should only be used for debugging).

      calibration label : RMSD < 2 A  AND  PoseBusters pass
      evaluation label  : RMSD < 1 A  (primary) / PLIF-contact recovery (secondary)

2. **Every knob is serialisable**, so a fitted model and the protocol can be
   frozen into one JSON file (reporting checklist D10 / O3).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict

__all__ = [
    "LabelSeparationError",
    "LabelConfig",
    "ScoreConfig",
    "CalibratorConfig",
    "MAESConfig",
    "default_config",
]

_EPS = 1e-9


class LabelSeparationError(ValueError):
    """Raised when calibration and evaluation labels are not independent."""


@dataclass
class LabelConfig:
    """Definitions of the calibration label and the independent evaluation labels."""

    # --- calibration label (used to fit Platt / isotonic parameters) ---
    calib_rmsd: float = 2.0
    calib_require_posebusters: bool = True

    # --- evaluation labels (never used for fitting) ---
    eval_rmsd: float = 1.0
    eval_plif_threshold: float = 0.8
    eval_primary: str = "rmsd_lt"          # 'rmsd_lt' | 'plif_recovery'
    eval_secondary: str = "plif_recovery"  # 'rmsd_lt' | 'plif_recovery'

    # --- escape hatch (debug only) ---
    allow_shared_threshold: bool = False

    # ------------------------------------------------------------------ checks
    def validate(self) -> "LabelConfig":
        if self.calib_rmsd <= 0 or self.eval_rmsd <= 0:
            raise LabelSeparationError("RMSD thresholds must be positive")
        if self.eval_primary == self.eval_secondary:
            raise LabelSeparationError(
                "eval_primary and eval_secondary must be different label families"
            )
        for name in (self.eval_primary, self.eval_secondary):
            if name not in {"rmsd_lt", "plif_recovery"}:
                raise LabelSeparationError(f"unknown evaluation label family: {name!r}")
        if self.shared_threshold and not self.allow_shared_threshold:
            raise LabelSeparationError(
                "calibration and evaluation RMSD thresholds are identical "
                f"({self.calib_rmsd} A). This makes MAES circular (red-team A2). "
                "Use eval_rmsd < calib_rmsd, or set allow_shared_threshold=True "
                "for debugging only."
            )
        if self.eval_rmsd > self.calib_rmsd:
            raise LabelSeparationError(
                "evaluation threshold must be at least as strict as the calibration "
                f"threshold ({self.eval_rmsd} A > {self.calib_rmsd} A)"
            )
        return self

    @property
    def shared_threshold(self) -> bool:
        return abs(float(self.calib_rmsd) - float(self.eval_rmsd)) < _EPS

    def independence_report(self) -> Dict[str, Any]:
        """Machine-readable description of the label protocol (goes into metadata)."""
        calib = f"RMSD < {self.calib_rmsd} A"
        if self.calib_require_posebusters:
            calib += " AND PoseBusters pass"
        eval_primary = (
            f"RMSD < {self.eval_rmsd} A"
            if self.eval_primary == "rmsd_lt"
            else f"PLIF/contact recovery >= {self.eval_plif_threshold}"
        )
        eval_secondary = (
            f"RMSD < {self.eval_rmsd} A"
            if self.eval_secondary == "rmsd_lt"
            else f"PLIF/contact recovery >= {self.eval_plif_threshold}"
        )
        return {
            "calibration_label": calib,
            "evaluation_primary": eval_primary,
            "evaluation_secondary": eval_secondary,
            "shared_threshold": bool(self.shared_threshold),
            "strictly_independent": not self.shared_threshold,
            "protocol": (
                "A2 fix: calibration label is built from RMSD+PoseBusters, evaluation "
                "labels are built from a stricter RMSD cut and from PLIF/contact "
                "recovery; the two label families are never derived from each other."
            ),
        }

    # ------------------------------------------------------------- (de)serialise
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "LabelConfig":
        return cls(**data)


@dataclass
class ScoreConfig:
    """Sub-score definitions (design doc section 3.2)."""

    # --- V: validity ---
    alpha: float = 0.2                    # weight of the strain term inside V
    strain_temperature_k: float = 298.15  # kT for the Boltzmann strain term
    strain_mode: str = "exp"             # 'exp' | 'squared'
    strain_sigma_kcal: float = 3.0       # only used when strain_mode == 'squared'

    # --- M: membrane consistency ---
    membrane_tau_angstrom: float = 4.0    # soft location tolerance
    membrane_lambda_kcal: float = 2.0     # hydrophobic-match tolerance
    membrane_angle_tau_deg: float = 30.0  # orientation tolerance
    interface_width_angstrom: float = 8.0
    membrane_weights: Dict[str, float] = field(
        default_factory=lambda: {"loc": 1 / 3, "hyd": 1 / 3, "ori": 1 / 3}
    )
    default_thickness_angstrom: float = 30.0

    # --- S: state consistency ---
    state_sigma: float = 1.0
    n_gmm_components: int = 2

    # --- combination & calibration plumbing ---
    combine: str = "geometric"           # 'geometric' | 'softmin' | 'arithmetic'
    softmin_p: float = 8.0
    weights: Dict[str, float] = field(
        default_factory=lambda: {"V": 1 / 3, "M": 1 / 3, "S": 1 / 3}
    )
    logit_eps: float = 1e-6
    random_state: int = 42

    def validate(self) -> "ScoreConfig":
        if not 0.0 <= float(self.alpha) <= 1.0:
            raise ValueError("alpha must be in [0, 1]")
        if self.combine not in {"geometric", "softmin", "arithmetic"}:
            raise ValueError(f"unknown combination rule: {self.combine!r}")
        if self.strain_mode not in {"exp", "squared"}:
            raise ValueError(f"unknown strain_mode: {self.strain_mode!r}")
        if self.strain_mode == "exp" and self.strain_temperature_k <= 0:
            raise ValueError("strain_temperature_k must be positive")
        if self.strain_mode == "squared" and self.strain_sigma_kcal <= 0:
            raise ValueError("strain_sigma_kcal must be positive")
        total = float(sum(self.weights.values()))
        if abs(total - 1.0) > 1e-6:
            raise ValueError(f"sub-score weights must sum to 1 (got {total})")
        if any(w <= 0 for w in self.weights.values()):
            raise ValueError("sub-score weights must be positive")
        if set(self.weights) != {"V", "M", "S"}:
            raise ValueError("weights must have keys V, M, S")
        m_total = float(sum(self.membrane_weights.values()))
        if abs(m_total - 1.0) > 1e-6:
            raise ValueError(f"membrane_weights must sum to 1 (got {m_total})")
        if self.n_gmm_components < 1:
            raise ValueError("n_gmm_components must be >= 1")
        return self

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ScoreConfig":
        return cls(**data)


@dataclass
class CalibratorConfig:
    """Calibration family used to map raw score -> P(acceptable)."""

    kind: str = "platt"   # 'platt' | 'isotonic' | 'none'
    l2: float = 1e-3

    def validate(self) -> "CalibratorConfig":
        if self.kind not in {"platt", "isotonic", "none"}:
            raise ValueError(f"unknown calibrator kind: {self.kind!r}")
        return self

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "CalibratorConfig":
        return cls(**data)


@dataclass
class MAESConfig:
    """Top-level configuration: labels + scores + calibrator."""

    labels: LabelConfig = field(default_factory=LabelConfig)
    scores: ScoreConfig = field(default_factory=ScoreConfig)
    calibrator: CalibratorConfig = field(default_factory=CalibratorConfig)
    version: str = "1.0.0"

    def validate(self) -> "MAESConfig":
        self.labels.validate()
        self.scores.validate()
        self.calibrator.validate()
        return self

    def to_dict(self) -> Dict[str, Any]:
        return {
            "version": self.version,
            "labels": self.labels.to_dict(),
            "scores": self.scores.to_dict(),
            "calibrator": self.calibrator.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "MAESConfig":
        return cls(
            labels=LabelConfig.from_dict(data.get("labels", {})),
            scores=ScoreConfig.from_dict(data.get("scores", {})),
            calibrator=CalibratorConfig.from_dict(data.get("calibrator", {})),
            version=data.get("version", "1.0.0"),
        )

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")
        return path

    @classmethod
    def load(cls, path: str | Path) -> "MAESConfig":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


def default_config() -> MAESConfig:
    """The frozen v1.0 protocol (label separation ON, equal sub-score weights)."""
    return MAESConfig().validate()
