"""MAES — Membrane-Aware Evaluation Score (reference implementation, v1.0.0).

MAES turns a membrane-protein ligand pose into a *calibrated probability* that
the pose is acceptable, by combining three sub-scores:

    V  validity          PoseBusters pass rate + RDKit ligand strain energy
    M  membrane consistency   membrane-plane projection + hydrophobic term + orientation
    S  state consistency  conformational-state descriptors + GMM classification

Design document: ``notes/maes_design.md`` (v0.1, 2026-11-16).
Red-team A2 (calibration/evaluation label circularity) is fixed by hard label
separation, enforced in :mod:`maes.config` and :mod:`maes.calibration`.

Public API (see :mod:`maes.api` for convenience wrappers)::

    from maes import default_config, MAESModel, load_placeholders
    cfg = default_config()
    model = MAESModel(cfg).fit(df)
    pred = model.predict(df)

Everything runs on CPU; no GPU or network access is required.
"""

from __future__ import annotations

__version__ = "1.0.0"
__all__ = [
    "__version__",
    "LabelConfig",
    "LabelSeparationError",
    "MAESConfig",
    "MAESModel",
    "ScoreConfig",
    "default_config",
    "load_placeholders",
]


def __getattr__(name: str):  # pragma: no cover - lazy re-exports
    if name in {"LabelConfig", "LabelSeparationError", "MAESConfig", "ScoreConfig", "default_config"}:
        from . import config
        return getattr(config, name)
    if name == "MAESModel":
        from .calibration import MAESModel
        return MAESModel
    if name == "load_placeholders":
        from .placeholders import load_placeholders
        return load_placeholders
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
