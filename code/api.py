"""Convenience API for MAES v1.0 (pipeline glue, smoke test, figure).

Typical use for Track B::

    from maes.api import build_features, fit_evaluate_smoke
    feats = build_features(pose_table, config=cfg, membrane_map=we_plane_table)
    model = MAESModel(cfg).fit(feats)
    report = evaluate_predictions(feats, model)

The smoke test (:func:`run_smoke`) exercises the whole chain on the three
placeholder bodies from :mod:`maes.placeholders`.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from .calibration import (
    MAESModel,
    assert_labels_independent,
    build_calibration_labels,
    build_evaluation_labels,
)
from .combine import entropy_weights
from .config import MAESConfig, default_config
from .evaluation import (
    ablation_study,
    evaluate_predictions,
    failure_count_matrix,
    failure_taxonomy,
    grouped_cv,
    weight_sensitivity,
)
from .irt import build_responses, fit_hierarchical_rasch, fit_rasch
from .membrane_map import (
    MembraneMap,
    NullMembraneMap,
    TrackBMembraneMap,
    build_trackb_membrane_map,
    load_frame_transforms,
)
from .placeholders import load_placeholders, placeholder_membrane_map
from .subscores import build_subscores

__all__ = [
    "build_features",
    "check_monotonicity",
    "fit_evaluate_smoke",
    "label_audit",
    "run_smoke",
]


def build_features(
    poses: pd.DataFrame,
    config: Optional[MAESConfig] = None,
    membrane_map: Optional[MembraneMap] = None,
    require_membrane: bool = False,
) -> pd.DataFrame:
    """Compute V/M/S for a Track-B pose table (see ``subscores.build_subscores``).

    Parameters
    ----------
    poses : DataFrame
        Primitive pose table; see :func:`maes.subscores.build_subscores` for the
        column contract (``pb_pass_rate``, ``strain_kcal``, ``z_lig_mean``,
        ``z_ref_mean``, ``angle_*``, ``dG_*``, ``phi_lig``/``phi_ref``).
    config : MAESConfig, optional
        Defaults to the frozen v1.0 protocol.
    membrane_map : MembraneMap, optional
        Plane provider; defaults to :class:`NullMembraneMap` (M neutral).
    require_membrane : bool
        If True, raise when any row lacks a usable membrane plane.
    """
    cfg = (config or default_config()).validate()
    provider = membrane_map if membrane_map is not None else NullMembraneMap()
    if "system_id" not in poses.columns:
        raise KeyError("pose table needs a 'system_id' column")
    features = build_subscores(poses, cfg.scores, membrane_map=provider)
    if require_membrane:
        missing = features.loc[features["m_neutral"], "system_id"].unique().tolist()
        if missing:
            raise ValueError(f"rows without a usable membrane plane: {missing}")
    return features


def label_audit(features: pd.DataFrame, config: Optional[MAESConfig] = None) -> Dict[str, Any]:
    """Explicit audit of the label-separation protocol (red-team A2 evidence)."""
    cfg = (config or default_config()).validate()
    calib = build_calibration_labels(features, cfg.labels)
    primary = build_evaluation_labels(features, cfg.labels, which="primary")
    secondary = build_evaluation_labels(features, cfg.labels, which="secondary")
    assert_labels_independent(calib, primary)
    assert_labels_independent(calib, secondary)
    from .calibration import label_fingerprint

    return {
        "protocol": cfg.labels.independence_report(),
        "n_calibration_positive": int(np.asarray(calib).sum()),
        "n_eval_primary_positive": int(np.asarray(primary).sum()),
        "n_eval_secondary_positive": int(np.asarray(secondary).sum()),
        "fingerprints": {
            "calibration": label_fingerprint(calib),
            "evaluation_primary": label_fingerprint(primary),
            "evaluation_secondary": label_fingerprint(secondary),
        },
        "independent": True,
    }


def check_monotonicity(
    model: MAESModel,
    n_steps: int = 11,
    low: float = 0.05,
    high: float = 1.0,
) -> Dict[str, Any]:
    """MAES must be non-decreasing in each sub-score when the others are held fixed."""
    grid = []
    levels = np.linspace(low, high, int(n_steps))
    for column in ("V", "M", "S"):
        for value in levels:
            row = {"V": 0.5, "M": 0.5, "S": 0.5}
            row[column] = float(value)
            grid.append(row)
    frame = pd.DataFrame(grid)
    scores = model.predict(frame)
    per_subscore: Dict[str, bool] = {}
    for i, column in enumerate(("V", "M", "S")):
        block = scores[i * n_steps:(i + 1) * n_steps]
        per_subscore[column] = bool(np.all(np.diff(block) >= -1e-12))
    return {
        "ok": bool(all(per_subscore.values())),
        "per_subscore": per_subscore,
        "maes_min": float(np.min(scores)),
        "maes_max": float(np.max(scores)),
        "in_unit_interval": bool(np.all((scores >= 0.0) & (scores <= 1.0))),
    }


def _smoke_figure(frame: pd.DataFrame, scores: np.ndarray, path: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    cfg = default_config()
    y_primary = build_evaluation_labels(frame, cfg.labels, which="primary")
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), dpi=200)

    # (a) reliability diagram
    from .evaluation import reliability_curve

    curve = reliability_curve(y_primary, scores, n_bins=8)
    ax = axes[0]
    ax.plot([0, 1], [0, 1], ls="--", lw=1, color="grey", label="perfect")
    valid = curve["count"] > 0
    ax.plot(curve.loc[valid, "mean_pred"], curve.loc[valid, "frac_pos"], "o-", ms=4, label="MAES")
    ax.set_xlabel("predicted P(acceptable)")
    ax.set_ylabel("observed fraction")
    ax.set_title("reliability (placeholder data)")
    ax.legend(frameon=False, fontsize=8)

    # (b) MAES vs RMSD
    ax = axes[1]
    for method, sub in frame.assign(MAES=scores).groupby("method"):
        ax.scatter(sub["rmsd"], sub["MAES"], s=18, alpha=0.8, label=str(method))
    ax.axvline(cfg.labels.eval_rmsd, color="grey", ls=":", lw=1)
    ax.axvline(cfg.labels.calib_rmsd, color="black", ls=":", lw=1)
    ax.set_xlabel("ligand RMSD (A)")
    ax.set_ylabel("MAES")
    ax.set_title("MAES vs RMSD (placeholder data)")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def run_smoke(
    figure_path: Optional[str | Path] = None,
    n_poses: int = 24,
    n_boot: int = 0,
    n_weight_samples: int = 200,
    seed: int = 42,
    verbose: bool = True,
) -> Dict[str, Any]:
    """End-to-end smoke test on the three placeholder bodies.

    Returns a JSON-serialisable report with: MAES range, monotonicity, label
    audit, held-out evaluation, ablation table, grouped CV, weight sensitivity,
    failure taxonomy counts and IRT ability/difficulty estimates.
    """
    cfg = default_config()
    poses = load_placeholders(n_poses=n_poses, seed=seed)
    provider = placeholder_membrane_map()
    features = build_features(poses, config=cfg, membrane_map=provider)

    audit = label_audit(features, cfg)
    model = MAESModel(cfg).fit(features)
    scores = model.predict(features)
    mono = check_monotonicity(model)

    evaluation = evaluate_predictions(features, model)
    ablation = ablation_study(features, cfg)
    try:
        cv = grouped_cv(features, cfg, group_col="system_id", n_splits=3)
        cv_summary = {
            "folds": cv.drop(columns=["error"], errors="ignore").to_dict(orient="records"),
            "mean_auroc": cv.attrs.get("mean_auroc"),
            "sd_auroc": cv.attrs.get("sd_auroc"),
        }
    except Exception as exc:  # pragma: no cover - defensive
        cv_summary = {"error": f"{type(exc).__name__}: {exc}"}

    sensitivity = weight_sensitivity(
        features, cfg, n_samples=int(n_weight_samples), seed=seed, top_k=3
    )
    taxonomy = failure_taxonomy(features, cfg)
    failure_matrix = failure_count_matrix(taxonomy).to_dict()

    responses = build_responses(features, label="calibration")
    irt_rasch = fit_rasch(responses, n_boot=int(n_boot), seed=seed)
    irt_eb = fit_hierarchical_rasch(responses, n_boot=0, seed=seed)

    report: Dict[str, Any] = {
        "maes_version": cfg.version,
        "n_poses": int(len(features)),
        "n_systems": int(features["system_id"].nunique()),
        "subscore_summary": {
            column: {
                "min": float(features[column].min()),
                "mean": float(features[column].mean()),
                "max": float(features[column].max()),
            }
            for column in ("V", "M", "S")
        },
        "m_neutral_fraction": float(features["m_neutral"].mean()),
        "s_neutral_fraction": float(features["s_neutral"].mean()),
        "maes_range": {"min": float(scores.min()), "max": float(scores.max())},
        "monotonicity": mono,
        "label_audit": audit,
        "model_metadata": model.metadata,
        "evaluation": evaluation,
        "ablation": ablation.to_dict(orient="records"),
        "grouped_cv": cv_summary,
        "weight_sensitivity": sensitivity,
        "entropy_weights": entropy_weights(features),
        "failure_counts": failure_matrix,
        "irt": {
            "rasch_abilities": irt_rasch.ability_frame().to_dict(orient="records"),
            "rasch_difficulties": irt_rasch.difficulty_frame().to_dict(orient="records"),
            "empirical_bayes_abilities": irt_eb.ability_frame().to_dict(orient="records"),
            "backend": irt_eb.backend,
        },
    }
    if figure_path is not None:
        out = _smoke_figure(features, scores, Path(figure_path))
        report["figure"] = str(out)
    if verbose:
        print(
            "MAES smoke: n=%d  MAES=[%.3f, %.3f]  monotone=%s  AUROC=%.3f  Brier=%.3f  ECE=%.3f  rho(RMSD)=%.3f"
            % (
                report["n_poses"],
                scores.min(),
                scores.max(),
                mono["ok"],
                evaluation.get("auroc_primary", float("nan")),
                evaluation.get("brier_primary", float("nan")),
                evaluation.get("ece_primary", float("nan")),
                evaluation.get("rho_maes_rmsd", float("nan")),
            )
        )
    return report


def fit_evaluate_smoke(config: Optional[MAESConfig] = None, seed: int = 42) -> Dict[str, Any]:
    """Small convenience wrapper returning (evaluation metrics) for tests/CI."""
    cfg = config or default_config()
    poses = load_placeholders(seed=seed)
    features = build_features(poses, config=cfg, membrane_map=placeholder_membrane_map())
    model = MAESModel(cfg).fit(features)
    return {"model": model, "features": features, "metrics": evaluate_predictions(features, model)}


# ============================================================================
# Track B data interface (C-line manifest, 2026-11-16)
# ----------------------------------------------------------------------------
# ``data/trackB/pose_manifest_curated.csv`` (68 rows) currently has:
#   target, pdb_id, ligand_id, smiles, mw, resolution, chain,
#   state_placeholder, ligand_class, notes, source, date
# It is a *reference-ligand* manifest.  MAES additionally needs one row per
# (ligand, method) predicted pose with the columns listed in
# ``TRACKB_PREDICTION_COLUMNS``.  ``build_features_from_manifest`` joins the two
# and fails loudly (with the exact missing columns) until B-line predictions land.
# ============================================================================

POSE_MANIFEST_COLUMNS = (
    "target", "pdb_id", "ligand_id", "smiles", "ligand_class",
    "state_placeholder", "resolution", "chain",
    # C-line W2 additions (verified on disk 2026-11-16):
    "crystal_form_cluster", "cluster_key", "cluster_size",
    "membrane_plane_source", "include_main",
)

# --- W2 primary / sensitivity analysis sets (Lead update 2026-11-16) --------
# main analysis = inhibitor-class poses (n=40, "just at threshold"), with all
# 68 poses (lipid + substrate_analogue included) as the sensitivity analysis.
PRIMARY_LIGAND_CLASS = "inhibitor"
N_PRIMARY_POSES = 40
PRIMARY_MEETS_THRESHOLD = True   # 40 >= 40; keep the benchmark wording

# --- three-tier analysis sets (Lead update 2026-11-16) ---------------------
# primary_main     : ligand_class == 'inhibitor'                 (currently 40)
# primary_extended : relaxed (<=3.5 A) extended inhibitor set   (target ~48)
# all              : every pose in the manifest                  (currently 68)
EXTENDED_MAX_RESOLUTION = 3.5
EXTRA_EXTENDED_CLASSES = ("inhibitor",)
ANALYSIS_SUBSETS = {
    "primary_main": "ligand_class == 'inhibitor' (benchmark threshold, n=40)",
    "primary_extended": (
        f"ligand_class in {EXTRA_EXTENDED_CLASSES} AND resolution <= {EXTENDED_MAX_RESOLUTION} A; "
        "or the C-line 'include_extended' flag when present (target n~48)"
    ),
    "all": "every pose in the manifest, incl. lipid/substrate_analogue (sensitivity set)",
}
DEFAULT_ANALYSIS_SUBSET = os.environ.get("MAES_ANALYSIS_SUBSET", "primary_main")
if DEFAULT_ANALYSIS_SUBSET not in ANALYSIS_SUBSETS:
    DEFAULT_ANALYSIS_SUBSET = "primary_main"

#: C-line W3 extended manifest (85 rows; ``inhibitor_main == TRUE`` -> n=48)
DEFAULT_EXTENDED_MANIFEST = (
    Path(__file__).resolve().parents[2] / "data" / "trackB" / "pose_manifest_extended_curated.csv"
)


def load_extended_manifest(path: Optional[str | Path] = None) -> Optional[pd.DataFrame]:
    """Load the C-line extended manifest (n=48 inhibitors) when present."""
    target = Path(path) if path is not None else DEFAULT_EXTENDED_MANIFEST
    if not target.exists():
        return None
    frame = pd.read_csv(target)
    return attach_pose_ids(frame) if "is_primary" not in frame.columns else frame


#: C-line Track B membrane planes + RCSB->OPM frame transforms (2026-11-16)
DEFAULT_PLANES_PATH = Path(__file__).resolve().parents[2] / "data" / "trackB" / "membrane_planes.json"
DEFAULT_TRANSFORMS_PATH = Path(__file__).resolve().parents[2] / "data" / "trackB" / "frame_transforms.csv"


def load_trackb_membrane_map(
    planes_path: Optional[str | Path] = None,
    transforms_path: Optional[str | Path] = None,
) -> Optional[TrackBMembraneMap]:
    """Dual-frame membrane map (OPM + RCSB-mapped) for Track B, or None if absent."""
    planes = Path(planes_path) if planes_path is not None else DEFAULT_PLANES_PATH
    if not planes.exists():
        return None
    transforms = Path(transforms_path) if transforms_path is not None else DEFAULT_TRANSFORMS_PATH
    return build_trackb_membrane_map(planes, transforms if transforms.exists() else None)


def _frame_norm(value) -> Optional[str]:
    text = str(value).strip().lower()
    if text in {"opm", "opm_oriented", "opm-oriented", "oriented"} or text.startswith("opm"):
        return "opm"
    if (
        text in {"rcsb", "rcsb_original", "rcsb-original", "original", "pdb"}
        or text.startswith("rcsb")
        or text.startswith("original")
    ):
        return "rcsb"
    return None


def resolve_membrane_keys(poses: pd.DataFrame, provider) -> pd.Series:
    """Resolve per-pose membrane keys ``f'{pdb_id}|{frame}'`` for a dual-frame map.

    Frame convention (C-line `frame_transforms.csv`): ``x_opm = R x_rcsb + t``.
    Poses whose receptor is in the OPM frame use ``pdb_id|opm``; poses in the
    RCSB original frame use ``pdb_id|rcsb`` (plane mapped by the loader).  When
    the pose frame is missing, the first available frame is used; when no plane
    is available at all the bare pdb_id is returned (M stays neutral).
    """
    if "pdb_id_pred" in poses.columns:
        pids = poses["pdb_id_pred"].astype(str)
    elif "pdb_id" in poses.columns:
        pids = poses["pdb_id"].astype(str)
    else:
        raise KeyError("need 'pdb_id_pred' or 'pdb_id' to resolve membrane keys")
    if "coord_frame" in poses.columns:
        frames = poses["coord_frame"]
    else:
        frames = pd.Series([None] * len(poses), index=poses.index)
    metadata = getattr(provider, "metadata", {}) or {}
    keys: list = []
    for pid, raw in zip(pids, frames):
        frame = None if pd.isna(raw) else _frame_norm(raw)
        info = metadata.get(str(pid), {})
        available = list(info.get("available_frames", []))
        chosen = None
        if frame is not None and hasattr(provider, "get_for_frame") and provider.get_for_frame(str(pid), frame) is not None:
            chosen = frame
        elif available:
            chosen = frame if frame in available else available[0]
        keys.append(f"{pid}|{chosen}" if chosen else str(pid))
    return pd.Series(keys, index=poses.index)


def membrane_coverage(
    manifest: pd.DataFrame,
    planes_path: Optional[str | Path] = None,
    transforms_path: Optional[str | Path] = None,
    membrane_map=None,
) -> Dict[str, Any]:
    """M-term coverage for a manifest's PDB entries (OPM / geometric / unavailable)."""
    provider = membrane_map or load_trackb_membrane_map(planes_path, transforms_path)
    if provider is None or "pdb_id" not in manifest.columns:
        return {"available": False}
    pids = sorted({str(v) for v in manifest["pdb_id"].astype(str) if str(v).strip()})
    metadata = getattr(provider, "metadata", {})
    by_source: Dict[str, int] = {}
    by_confidence: Dict[str, int] = {}
    frames: Dict[str, int] = {}
    no_plane: list = []
    for pid in pids:
        info = metadata.get(pid, {})
        source = str(info.get("membrane_source", "unknown"))
        by_source[source] = by_source.get(source, 0) + 1
        conf = str(info.get("membrane_confidence", "none"))
        by_confidence[conf] = by_confidence.get(conf, 0) + 1
        available = list(info.get("available_frames", []))
        for frame in available:
            frames[frame] = frames.get(frame, 0) + 1
        if not available:
            no_plane.append(pid)
    return {
        "available": True,
        "n_pdb_entries": len(pids),
        "by_membrane_source": by_source,
        "by_confidence": by_confidence,
        "by_available_frame": frames,
        "n_without_any_plane": len(no_plane),
        "pdb_ids_without_plane": no_plane[:20],
        "note": (
            "coverage over the manifest's PDB entries; per-pose coverage is carried in "
            "the feature table via m_neutral / m_confidence / m_frame_mismatch"
        ),
    }

# Cluster definitions are deliberately switchable.  `target_spacegroup` becomes
# available automatically once the C-line exports a `space_group` column (from
# PDB CRYST1); until then the text can quote the per-entry / state / target range.
CLUSTER_DEFINITIONS = {
    "per_entry": ("target", "pdb_id"),                       # every PDB entry independent (upper bound)
    "crystal_entry": ("target", "crystal_form_cluster"),     # C-line per-entry cluster label
    "target_spacegroup": ("target", "space_group"),          # target x space-group
    "target_crystal_form": ("target", "crystal_form"),       # target x CF (space group + cell clustering)
    "target_cf_subsite": ("target", "crystal_form", "sub_site"),  # target x CF x sub-site (AcrB CF01: 5 pockets)
    "target_state": ("target", "state_placeholder"),         # crystal-form proxy where annotated
    "target": ("target",),                                    # fully clustered (lower bound)
}
DEFAULT_CLUSTER_DEFINITION = "crystal_entry"
# Lead 2026-11-16 (final n_eff rule): cluster = target x crystal form x sub-site
# (`acrb_cf01_pose_pockets.csv` splits CF01 into 5 pockets -> ~15 clusters for
# the 40 inhibitor poses); falls back to target x CF, then space-group, then
# per-entry when the finer columns are unavailable.
PREFERRED_CLUSTER_ORDER = (
    "target_cf_subsite",
    "target_crystal_form",
    "target_spacegroup",
    "crystal_entry",
)

#: C-line sub-site clustering for AcrB CF01 (17 poses -> 5 pockets, 2026-11-16)
DEFAULT_ACRB_POCKETS_PATH = Path(__file__).resolve().parents[2] / "data" / "trackB" / "acrb_cf01_pose_pockets.csv"
_SUBSITE_CACHE: Dict[str, Optional[pd.DataFrame]] = {}


def load_acrb_subsites(path: Optional[str | Path] = None) -> Optional[pd.DataFrame]:
    """Load the CF01 sub-site assignment table (pdb_id, ligand_id -> sub_site)."""
    target = Path(path) if path is not None else DEFAULT_ACRB_POCKETS_PATH
    key = str(target)
    if key in _SUBSITE_CACHE:
        return _SUBSITE_CACHE[key]
    if not target.exists():
        _SUBSITE_CACHE[key] = None
        return None
    frame = pd.read_csv(target)
    if "centroid_cluster_10A" in frame.columns:
        frame = frame.rename(columns={"centroid_cluster_10A": "sub_site"})
    cols = [c for c in ("pdb_id", "ligand_id", "sub_site") if c in frame.columns]
    if len(cols) != 3:
        raise KeyError("sub-site table needs pdb_id / ligand_id / centroid_cluster_10A(sub_site)")
    frame = frame[cols].drop_duplicates(subset=["pdb_id", "ligand_id"]).copy()
    frame["pdb_id"] = frame["pdb_id"].astype(str)
    frame["ligand_id"] = frame["ligand_id"].astype(str)
    frame["sub_site"] = frame["sub_site"].astype(str)
    _SUBSITE_CACHE[key] = frame
    return frame


def attach_subsites(manifest: pd.DataFrame, path: Optional[str | Path] = None) -> pd.DataFrame:
    """Attach the ``sub_site`` column (AcrB CF01 pocket id; 'unknown' elsewhere)."""
    out = manifest.copy()
    if "sub_site" in out.columns:
        return out
    pockets = load_acrb_subsites(path)
    if pockets is None or len(pockets) == 0:
        out["sub_site"] = "unknown"
        return out
    out["pdb_id"] = out["pdb_id"].astype(str)
    out["ligand_id"] = out["ligand_id"].astype(str)
    out = out.merge(pockets, on=["pdb_id", "ligand_id"], how="left")
    out["sub_site"] = out["sub_site"].fillna("unknown").astype(str)
    return out


def set_default_analysis_subset(subset: str) -> str:
    """Switch the default tier (Lead decision point); returns the active tier."""
    if subset not in ANALYSIS_SUBSETS:
        raise ValueError(f"unknown subset {subset!r}; choose from {sorted(ANALYSIS_SUBSETS)}")
    global DEFAULT_ANALYSIS_SUBSET
    DEFAULT_ANALYSIS_SUBSET = subset
    return DEFAULT_ANALYSIS_SUBSET

TRACKB_PREDICTION_COLUMNS = (
    # --- required ---------------------------------------------------------
    "pose_key",          # str  f"{target}|{pdb_id}|{ligand_id}|{method}" (unique)
    "method",            # str  vina | unidock | gnina | boltz2 | af3 | ...
    "rmsd",              # float Angstrom, symmetry-aware best match vs co-crystal reference
    "plif_recovery",     # float [0,1], reference contacts recovered (evaluation label family)
    "posebusters_pass",  # bool, full PoseBusters pass (calibration label family)
    # --- strongly recommended --------------------------------------------
    "pb_pass_rate",      # float [0,1], PoseBusters partial pass rate (else leave NA)
    "strain_kcal",       # float kcal/mol, RDKit UFF/MMFF strain vs relaxed conformer
    "z_lig_mean",        # float Angstrom, mean projected ligand z in the OPM-oriented frame
    "z_ref_mean",        # float Angstrom, reference ligand z (same frame)
    "angle_lig_deg",     # float degrees, ligand principal axis vs membrane normal
    "angle_ref_deg",     # float degrees, reference value
    "dG_lig",            # float, simplified insertion-energy proxy (kcal/mol)
    "dG_ref",            # float, reference value
    "phi_lig",           # JSON list, state descriptors e.g. [tm_angle, pocket_void_A3, iface_min_A]
    "phi_ref",           # JSON list, reference descriptors
    # --- provenance / QC (recommended) -----------------------------------
    "crystal_form",      # str  target x CF cluster id (falls back to manifest)
    "receptor_source",   # str  holo | apo | predicted
    "score",             # float native method score (affinity/pLDDT/ipTM ...)
    "pose_rank",         # int   1 = best pose of that ligand
    "method_version",    # str   e.g. vina 1.2.5, boltz-2.2
    "seed",              # int   docking/co-folding seed
    "status",            # str   ok | oom | failed
    "error_message",     # str   required when status != ok
)

TRACKB_REQUIRED_COLUMNS = ("pose_key", "method", "rmsd", "plif_recovery", "posebusters_pass")

# Lead G0 decision (2026-11-16): MAES main analysis runs on Vina docking poses
# for all 40 inhibitor systems; Boltz-2 (<=400 tokens) is a separate secondary arm.
PRIMARY_METHODS = ("vina",)
SECONDARY_METHODS = ("boltz2", "af3", "boltz-2")


def select_method_subset(predictions: pd.DataFrame, methods: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Filter a prediction table to one method arm (default: the Vina primary arm)."""
    if "method" not in predictions.columns:
        raise KeyError("prediction table needs a 'method' column")
    wanted = tuple(str(m).strip().lower() for m in (methods if methods is not None else PRIMARY_METHODS))
    mask = predictions["method"].astype(str).str.strip().str.lower().isin(wanted)
    return predictions.loc[mask].reset_index(drop=True)


# --- H2 predicted-receptor QC (preregistration amendments #1-#4, 2026-11-16)
H2_PREDICTED_SOURCES = ("predicted_afdb",)
H2_QC_MIN_TM_SCORE = 0.70
H2_QC_MAX_EXCLUDED = 4          # >4/40 systems -> H2 downgraded to exploratory
AFDB_VERSION = "v6"             # Amendment #3: current release at analysis time
# Amendment #4: AcrB domain constructs are aligned on two non-contiguous
# segments of the full-length AFDB model; the H2 main table allows
# {region, multi_segment}, while `fullchain` is SI-only.
H2_ALIGNMENT_MODES = ("region", "multi_segment")
H2_SI_ALIGNMENT_MODES = ("fullchain",)
H2_MULTI_SEGMENT_COLUMNS = ("tm_score_segment1", "tm_score_segment2")
H2_RECEPTOR_QC_COLUMNS = (
    "receptor_source", "receptor_id", "tm_score", "pocket_rmsd", "plddt_pocket",
    # Amendment #2: AFDB models are restricted to the SIFTS-mapped residue range
    "uniprot_start", "uniprot_end", "afdb_region_id",
)
H2_RECEPTOR_PROVENANCE_COLUMNS = (
    "afdb_version", "afdb_entry_id", "afdb_model_created", "afdb_model_sha256",
    "afdb_url", "plddt_global", "afdb_alignment_mode", "sifts_multisegment",
    "coord_frame",
)


def check_h2_receptor_qc(
    predictions: pd.DataFrame,
    min_tm_score: float = H2_QC_MIN_TM_SCORE,
    max_excluded: int = H2_QC_MAX_EXCLUDED,
    source: str = "predicted_afdb",
) -> Dict[str, Any]:
    """QC gate for the H2 predicted-receptor arm (preregistration amendments #1-#4).

    Unit = one inhibitor system (one ``pose_key``).  Two alignment modes are
    allowed in the H2 main table (Amendment #4):

    * ``region`` — AFDB v6 trimmed to the SIFTS-mapped residue range; fails when
      ``tm_score < 0.70``;
    * ``multi_segment`` — AcrB domain constructs spanning two non-contiguous
      segments; the full-length AFDB model is aligned on both segments and the
      system fails only when BOTH segment TM-scores are < 0.70 (per-segment
      values in ``tm_score_segment1/2`` are required).

    Additional failures: ``excluded_predicted`` flag / ``exclusion_reason``,
    missing QC columns, missing/invalid SIFTS range, or a ``fullchain`` mode
    (SI-only).  If more than ``max_excluded`` systems fail, H2 must be
    downgraded to exploratory (preregistered stop rule).
    """
    if "pose_key" not in predictions.columns:
        raise KeyError("predictions need a 'pose_key' column")
    missing_cols = [c for c in H2_RECEPTOR_QC_COLUMNS if c not in predictions.columns]
    missing_provenance = [c for c in H2_RECEPTOR_PROVENANCE_COLUMNS if c not in predictions.columns]
    pred = predictions.copy()
    if "receptor_source" in pred.columns:
        pred = pred[pred["receptor_source"].astype(str).str.strip().str.lower() == str(source).lower()]
    versions_seen = (
        sorted(pred["afdb_version"].dropna().astype(str).str.strip().unique().tolist())
        if "afdb_version" in pred.columns
        else []
    )
    coord_frames_seen = (
        sorted(pred["coord_frame"].dropna().astype(str).str.strip().unique().tolist())
        if "coord_frame" in pred.columns
        else []
    )
    n_missing_coord_frame = (
        int(pred["coord_frame"].isna().sum()) if "coord_frame" in pred.columns else int(len(pred))
    )
    systems = pred["pose_key"].astype(str).drop_duplicates().tolist()
    rows: list = []
    for key in systems:
        sub = pred[pred["pose_key"].astype(str) == key]
        mode = "region"
        if "afdb_alignment_mode" in sub.columns:
            values = sub["afdb_alignment_mode"].dropna().astype(str).str.strip().str.lower()
            if len(values):
                mode = str(values.iloc[0])
        tm = pd.to_numeric(sub["tm_score"], errors="coerce").max() if "tm_score" in sub.columns else float("nan")
        seg1 = (
            pd.to_numeric(sub["tm_score_segment1"], errors="coerce").max()
            if "tm_score_segment1" in sub.columns else float("nan")
        )
        seg2 = (
            pd.to_numeric(sub["tm_score_segment2"], errors="coerce").max()
            if "tm_score_segment2" in sub.columns else float("nan")
        )
        flagged = bool(sub["excluded_predicted"].fillna(False).astype(bool).any()) if "excluded_predicted" in sub.columns else False
        reason = ""
        if "exclusion_reason" in sub.columns:
            reasons = sub["exclusion_reason"].dropna().astype(str).str.strip()
            reasons = reasons[reasons != ""]
            if len(reasons):
                reason = "; ".join(reasons.unique()[:2])
        if missing_cols:
            passed, reason = False, reason or ("missing QC columns: " + ", ".join(missing_cols))
        elif mode in H2_SI_ALIGNMENT_MODES:
            passed = False
            reason = reason or f"afdb_alignment_mode={mode!r} is SI-only (not allowed in H2 main table)"
        elif mode not in H2_ALIGNMENT_MODES:
            passed = False
            reason = reason or f"invalid afdb_alignment_mode={mode!r}"
        elif flagged:
            passed = False
        elif mode == "multi_segment":
            # Amendment #4: two-segment alignment; only exclude when BOTH
            # segment TM-scores are below threshold.
            seg_missing = [c for c in H2_MULTI_SEGMENT_COLUMNS if c not in sub.columns]
            if seg_missing:
                passed = False
                reason = reason or ("missing multi-segment TM columns: " + ", ".join(seg_missing))
            elif not (np.isfinite(seg1) and np.isfinite(seg2)):
                passed = False
                reason = reason or "multi_segment requires finite tm_score_segment1/2"
            elif max(seg1, seg2) < float(min_tm_score):
                passed = False
                reason = reason or f"both segment TM-scores < {min_tm_score} (seg1={seg1}, seg2={seg2})"
            else:
                ranges = str(sub["uniprot_ranges"].iloc[0]).strip() if "uniprot_ranges" in sub.columns else ""
                region_id = str(sub["afdb_region_id"].iloc[0]).strip() if "afdb_region_id" in sub.columns else ""
                ranges_ok = bool(ranges and ranges.lower() not in {"nan", "none", ""})
                region_ok = bool(region_id and region_id.lower() not in {"nan", "none", ""})
                passed = bool(ranges_ok and region_ok)
                if not passed:
                    reason = reason or (
                        "missing or invalid SIFTS multi-segment range "
                        "(uniprot_ranges / afdb_region_id)"
                    )
        else:  # region mode (Amendment #2)
            if (not np.isfinite(tm)) or float(tm) < float(min_tm_score):
                passed = False
                reason = reason or (f"tm_score={tm}" if np.isfinite(tm) else "tm_score missing")
            else:
                try:
                    start = float(pd.to_numeric(sub["uniprot_start"], errors="coerce").iloc[0])
                    end = float(pd.to_numeric(sub["uniprot_end"], errors="coerce").iloc[0])
                    region_id = str(sub["afdb_region_id"].iloc[0]).strip()
                except Exception:
                    start = end = float("nan")
                    region_id = ""
                region_ok = bool(
                    np.isfinite(start) and np.isfinite(end) and end >= start
                    and region_id and region_id.lower() not in {"nan", "none", ""}
                )
                if region_ok:
                    passed = True
                else:
                    passed = False
                    reason = reason or (
                        "missing or invalid SIFTS residue range "
                        "(uniprot_start / uniprot_end / afdb_region_id)"
                    )
        rows.append(
            {
                "pose_key": key,
                "alignment_mode": mode,
                "tm_score": (None if not np.isfinite(tm) else float(tm)),
                "tm_score_segment1": (None if not np.isfinite(seg1) else float(seg1)),
                "tm_score_segment2": (None if not np.isfinite(seg2) else float(seg2)),
                "qc_passed": bool(passed),
                "exclusion_reason": reason,
            }
        )
    qc = pd.DataFrame(rows)
    n_pred = int(len(qc))
    n_passed = int(qc["qc_passed"].sum()) if n_pred else 0
    n_excluded = n_pred - n_passed
    without_region = int(
        qc["exclusion_reason"].astype(str).str.contains("SIFTS").sum()
    ) if n_pred else 0
    modes = qc["alignment_mode"].astype(str) if n_pred else pd.Series(dtype=str)
    return {
        "n_predicted_systems": n_pred,
        "n_qc_passed": n_passed,
        "n_excluded": n_excluded,
        "excluded_fraction": float(n_excluded / n_pred) if n_pred else float("nan"),
        "max_excluded_allowed": int(max_excluded),
        "downgrade_h2_exploratory": bool(n_excluded > int(max_excluded)),
        "n_region_mode": int((modes == "region").sum()),
        "n_multi_segment_mode": int((modes == "multi_segment").sum()),
        "alignment_modes_allowed": list(H2_ALIGNMENT_MODES),
        "missing_qc_columns": missing_cols,
        "missing_provenance_columns": missing_provenance,
        "afdb_versions_seen": versions_seen,
        "expected_afdb_version": AFDB_VERSION,
        "coord_frames_seen": coord_frames_seen,
        "n_missing_coord_frame": n_missing_coord_frame,
        "n_without_sifts_range": without_region,
        "min_tm_score": float(min_tm_score),
        "region_fields": ("uniprot_start", "uniprot_end", "afdb_region_id"),
        "multi_segment_fields": ("tm_score_segment1", "tm_score_segment2", "uniprot_ranges"),
        "per_system": qc.to_dict(orient="records"),
    }

MISSING_TRACKB_INTERFACE = {
    "provided_by_C": ["target", "pdb_id", "ligand_id", "smiles", "ligand_class",
                      "state_placeholder", "resolution", "chain",
                      "membrane planes (system_id -> normal/center/thickness)"],
    "needed_from_B": list(TRACKB_PREDICTION_COLUMNS),
    "n_eff_note": (
        "AcrB contributes 55/68 reference poses but few independent crystal forms; "
        "cluster bootstrap on crystal_form and the grouped IRT term "
        "(irt.fit_rasch_grouped) are the R3 mitigations; report n_eff ~ 20-35."
    ),
}


class MissingTrackBDataError(RuntimeError):
    """Raised when the Track B pose predictions are not available yet."""


def load_pose_manifest(path: str | Path) -> pd.DataFrame:
    """Load ``data/trackB/pose_manifest_curated.csv`` and attach join keys."""
    manifest = pd.read_csv(Path(path))
    return attach_pose_ids(manifest)


def attach_pose_ids(manifest: pd.DataFrame) -> pd.DataFrame:
    """Add analysis-set flags, join keys and cluster ids to a reference manifest."""
    out = manifest.copy()
    for column in ("target", "pdb_id", "ligand_id"):
        if column not in out.columns:
            raise KeyError(f"pose manifest needs a {column!r} column")
    if "ligand_class" not in out.columns:
        raise KeyError("pose manifest needs a 'ligand_class' column for the W2 analysis sets")
    ligand_class = out["ligand_class"].astype(str).str.strip().str.lower()
    out["ligand_class_norm"] = ligand_class
    out["is_primary"] = ligand_class == PRIMARY_LIGAND_CLASS
    out["system_id"] = out["target"].astype(str)
    out["pair_id"] = (
        out["target"].astype(str) + "|" + out["pdb_id"].astype(str) + "|" + out["ligand_id"].astype(str)
    )
    if "pose_id" in out.columns:
        out["pose_id"] = out["pose_id"].astype(str)
    else:
        out["pose_id"] = out["pair_id"]
    # C-line CF01 sub-site assignment (AcrB): 17 CF01 poses -> 5 pockets; other
    # entries stay 'unknown' so the cluster definition degrades gracefully
    out = attach_subsites(out)
    # cluster ids for every *available* definition (NaN-safe); definitions whose
    # columns are missing (e.g. space_group until C-line exports it) are skipped
    built = []
    for name, columns in CLUSTER_DEFINITIONS.items():
        if not all(column in out.columns for column in columns):
            continue
        parts = []
        for column in columns:
            values = out[column].astype(str).str.strip()
            values = values.where(~values.isin({"nan", "None", ""}), "unknown")
            parts.append(values)
        key = parts[0]
        for extra in parts[1:]:
            key = key + "|" + extra
        out[f"cluster_{name}"] = key
        built.append(name)
    out.attrs["cluster_definitions_available"] = built
    if f"cluster_{DEFAULT_CLUSTER_DEFINITION}" in out.columns:
        out["cluster_id"] = out[f"cluster_{DEFAULT_CLUSTER_DEFINITION}"]
    elif "cluster_per_entry" in out.columns:
        out["cluster_id"] = out["cluster_per_entry"]
    out["cluster_key_c"] = out.get("cluster_key", out["cluster_id"]).astype(str)
    return out


def select_analysis_set(
    manifest: pd.DataFrame,
    subset: Optional[str] = None,
    extended_manifest: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Three-tier W2 analysis sets (C-line W3 schema, Lead 2026-11-16).

    primary_main     ``inhibitor_main == TRUE`` (n=40; C-line flag, fallback ligand_class)
    primary_extended C-line extended manifest ``inhibitor_main == TRUE`` (n=48)
    all              every pose in the main manifest (n=68 sensitivity set)
    """
    out = manifest if "is_primary" in manifest.columns else attach_pose_ids(manifest)
    tier = subset or DEFAULT_ANALYSIS_SUBSET
    if tier == "primary_main":
        if "inhibitor_main" in out.columns:
            return out[out["inhibitor_main"].fillna(False).astype(bool)].reset_index(drop=True)
        return out[out["is_primary"]].reset_index(drop=True)
    if tier == "primary_extended":
        extended = extended_manifest if extended_manifest is not None else load_extended_manifest()
        if extended is not None and len(extended):
            if "inhibitor_main" in extended.columns:
                return extended[extended["inhibitor_main"].fillna(False).astype(bool)].reset_index(drop=True)
            return extended[extended["is_primary"]].reset_index(drop=True)
        # fallback: relax the resolution cut on the main manifest only
        if "include_extended" in out.columns:
            mask = out["is_primary"] | out["include_extended"].fillna(False).astype(bool)
        elif "resolution" in out.columns:
            resolution = pd.to_numeric(out["resolution"], errors="coerce")
            relaxed = resolution <= float(EXTENDED_MAX_RESOLUTION)
            extended_classes = [str(c).strip().lower() for c in EXTRA_EXTENDED_CLASSES]
            classes = out["ligand_class_norm"].isin(extended_classes)
            mask = out["is_primary"] | (classes & relaxed)
        else:
            mask = out["is_primary"]
        return out[mask].reset_index(drop=True)
    if tier in {"all", "sensitivity"}:
        return out.reset_index(drop=True)
    raise ValueError(f"unknown subset: {tier!r}; choose from {sorted(ANALYSIS_SUBSETS)}")


def cluster_keys(manifest: pd.DataFrame, definition: str = DEFAULT_CLUSTER_DEFINITION) -> pd.Series:
    """Cluster key series for one of :data:`CLUSTER_DEFINITIONS` (must be available)."""
    if definition not in CLUSTER_DEFINITIONS:
        raise ValueError(f"unknown cluster definition {definition!r}")
    column = f"cluster_{definition}"
    if column not in manifest.columns:
        manifest = attach_pose_ids(manifest)
    if column not in manifest.columns:
        available = manifest.attrs.get("cluster_definitions_available", [])
        raise KeyError(
            f"cluster definition {definition!r} is unavailable: the manifest lacks column(s) "
            f"{CLUSTER_DEFINITIONS[definition]}; available definitions: {available}"
        )
    return manifest[column].astype(str)


def _kish_n_eff(sizes: np.ndarray) -> float:
    sizes = np.asarray(sizes, dtype=float)
    if sizes.size == 0:
        return float("nan")
    return float((sizes.sum() ** 2) / np.square(sizes).sum())


def analysis_set_report(
    manifest: pd.DataFrame,
    cluster_definition: str = DEFAULT_CLUSTER_DEFINITION,
    membrane_map=None,
) -> Dict[str, Any]:
    """n_obs / n_clusters / n_eff for every analysis tier and cluster definition.

    Every tier reports its own sample size, cluster count and Kish n_eff, plus
    the tier rule, so the manuscript can quote "n_obs=40, n_clusters=K,
    n_eff=... (target x space-group cluster bootstrap)" without ambiguity.
    """
    out = attach_pose_ids(manifest)
    extended_manifest = load_extended_manifest()
    available = list(out.attrs.get("cluster_definitions_available", []))
    preferred = cluster_definition
    for candidate in PREFERRED_CLUSTER_ORDER:
        if candidate in available:
            preferred = candidate
            break
    report: Dict[str, Any] = {
        "preferred_cluster_definition": preferred,
        "available_cluster_definitions": available,
        "default_subset": DEFAULT_ANALYSIS_SUBSET,
        "main_set_flag": "inhibitor_main" if "inhibitor_main" in out.columns else "ligand_class",
        "extended_source": (
            "C-line pose_manifest_extended_curated.csv (inhibitor_main)"
            if extended_manifest is not None
            else "fallback resolution rule on the main manifest"
        ),
        "subsets": {},
        "analysis_sets": {tier: rule for tier, rule in ANALYSIS_SUBSETS.items()},
        "warnings": [],
    }
    for tier in ANALYSIS_SUBSETS:
        frame = select_analysis_set(out, tier, extended_manifest=extended_manifest)
        tier_report: Dict[str, Any] = {
            "subset": tier,
            "rule": ANALYSIS_SUBSETS[tier],
            "n_obs": int(len(frame)),
            "by_target": {str(k): int(v) for k, v in frame["target"].value_counts().items()},
            "acrb_share": float((frame["target"].astype(str) == "AcrB_Ec").mean()) if len(frame) else float("nan"),
            "clusters": {},
        }
        if "include_main" in frame.columns:
            tier_report["include_main_flag_count"] = int(frame["include_main"].fillna(False).sum())
        for definition in available:
            keys = cluster_keys(frame, definition).value_counts()
            tier_report["clusters"][definition] = {
                "n_clusters": int(keys.size),
                "max_cluster": int(keys.max()) if keys.size else 0,
                "kish_n_eff": _kish_n_eff(keys.to_numpy()),
            }
        preferred_stats = tier_report["clusters"].get(preferred)
        if preferred_stats is not None:
            tier_report["n_clusters"] = preferred_stats["n_clusters"]
            tier_report["n_eff"] = preferred_stats["kish_n_eff"]
            tier_report["n_eff_note"] = (
                "Kish n_eff = (sum m)^2 / sum m^2 over clusters; report ICC-corrected n_eff "
                "(evaluation.design_effect_n_eff) once MAES acceptable indicators exist."
            )
        report["subsets"][tier] = tier_report

    main = report["subsets"]["primary_main"]
    extended = report["subsets"]["primary_extended"]
    report["n_primary"] = main["n_obs"]
    report["n_sensitivity"] = report["subsets"]["all"]["n_obs"]
    report["n_extended"] = extended["n_obs"]
    report["membrane_coverage"] = membrane_coverage(out, membrane_map=membrane_map)
    report["primary_ligand_class"] = PRIMARY_LIGAND_CLASS
    report["primary_meets_benchmark_threshold"] = bool(main["n_obs"] >= N_PRIMARY_POSES)
    report["threshold"] = int(N_PRIMARY_POSES)
    report["primary_by_target"] = dict(main["by_target"])
    report["acrb_share_primary"] = float(main["acrb_share"])
    report["cluster_definition"] = preferred
    report["clusters"] = main["clusters"]  # backward-compatible structure
    report["include_main_flag_count"] = main.get("include_main_flag_count")
    # convenience headline numbers for the manuscript text
    report["n_obs"] = main["n_obs"]
    report["n_clusters"] = main.get("n_clusters")
    report["n_eff"] = main.get("n_eff")
    report["n_eff_note"] = main.get("n_eff_note")
    if main.get("include_main_flag_count", main["n_obs"]) != main["n_obs"]:
        report["warnings"].append(
            f"C-line include_main flags only {main['include_main_flag_count']} of {main['n_obs']} "
            "inhibitor poses; the W2 primary set follows the Lead instruction (all inhibitor-class poses)."
        )
    if extended["n_obs"] < N_PRIMARY_POSES:
        report["warnings"].append(
            f"primary_extended currently has {extended['n_obs']} poses (< {N_PRIMARY_POSES}): with the "
            "present manifest the <=3.5 A rule only relaxes resolution, it does not reclassify ligands. "
            "C-line must either expand the inhibitor classification to ~48 poses or export an "
            "'include_extended' flag before this tier can replace primary_main."
        )
    if "target_spacegroup" not in available:
        report["warnings"].append(
            "target x space_group clustering is not available yet (no space_group column). "
            "Waiting on C-line to parse PDB CRYST1; meanwhile report target_state and target "
            "as lower bounds and crystal_entry as the upper bound."
        )
    if extended["n_obs"] <= main["n_obs"]:
        report["warnings"].append(
            f"primary_extended currently adds no poses over primary_main ({extended['n_obs']} vs "
            f"{main['n_obs']}): the <=3.5 A rule cannot reclassify ligands by itself. C-line must "
            "expand the inhibitor classification to ~48 poses or export an 'include_extended' flag "
            "before this tier replaces primary_main."
        )
    return report


def check_trackb_inputs(
    manifest: pd.DataFrame,
    predictions: Optional[pd.DataFrame] = None,
    min_ok_with_rmsd: Optional[int] = None,
) -> Dict[str, Any]:
    """Report what is still missing before MAES can run on Track B.

    ``min_ok_with_rmsd`` optionally requires at least that many successful rows
    with a non-missing ``rmsd`` (final Track B table: 40 rows, 2 failed and 4
    ``rmsd=NA(no_reference_coords)`` rows are retained for audit but excluded
    from calibration/evaluation; ingest passes 30).  Missing/empty-column checks
    are unchanged.
    """
    reported: Dict[str, Any] = {
        "manifest_rows": int(len(manifest)),
        "manifest_columns": list(manifest.columns),
        "n_targets": int(manifest["target"].nunique()) if "target" in manifest else None,
        "n_ligand_classes": int(manifest["ligand_class"].nunique())
        if "ligand_class" in manifest.columns else None,
    }
    if predictions is None:
        reported["ready"] = False
        reported["missing"] = list(TRACKB_REQUIRED_COLUMNS)
        reported["missing_recommended"] = [
            c for c in TRACKB_PREDICTION_COLUMNS if c not in TRACKB_REQUIRED_COLUMNS
        ]
        return reported
    missing_required = [c for c in TRACKB_REQUIRED_COLUMNS if c not in predictions.columns]
    missing_recommended = [
        c for c in TRACKB_PREDICTION_COLUMNS
        if c not in TRACKB_REQUIRED_COLUMNS and c not in predictions.columns
    ]
    # "present but empty" required columns (metrics-only re-delivery, e.g. a
    # prediction table whose rmsd / posebusters_pass values are still pending)
    empty_required = []
    for column in TRACKB_REQUIRED_COLUMNS:
        if column in predictions.columns:
            series = predictions[column]
            if len(series) and int(series.notna().sum()) == 0:
                empty_required.append(column)
    reported["prediction_rows"] = int(len(predictions))
    reported["missing"] = missing_required
    reported["empty_required"] = empty_required
    reported["missing_recommended"] = missing_recommended
    # how many successful rows actually carry an RMSD (NA + reason rows are
    # retained for audit but cannot enter calibration/evaluation)
    n_ok_with_rmsd: Optional[int] = None
    if "rmsd" in predictions.columns:
        if "status" in predictions.columns:
            ok_mask = predictions["status"].astype(str).str.strip().str.lower().eq("ok")
        else:
            ok_mask = pd.Series(True, index=predictions.index)
        n_ok_with_rmsd = int((ok_mask & predictions["rmsd"].notna()).sum())
    reported["n_ok_with_rmsd"] = n_ok_with_rmsd
    reported["min_ok_with_rmsd"] = min_ok_with_rmsd
    rmsd_coverage_ok = (
        True if min_ok_with_rmsd is None or n_ok_with_rmsd is None
        else bool(n_ok_with_rmsd >= int(min_ok_with_rmsd))
    )
    reported["rmsd_coverage_ok"] = bool(rmsd_coverage_ok)
    reported["ready"] = (not missing_required) and (not empty_required) and reported["rmsd_coverage_ok"]
    if "pose_key" in predictions.columns:
        reported["n_methods"] = int(predictions["method"].nunique()) if "method" in predictions else None
    return reported


def build_features_from_manifest(
    manifest: pd.DataFrame,
    predictions: Optional[pd.DataFrame] = None,
    config: Optional[MAESConfig] = None,
    membrane_map: Optional[MembraneMap] = None,
    subset: Optional[str] = None,
) -> pd.DataFrame:
    """Join the C-line manifest with B-line pose predictions and compute V/M/S.

    ``subset`` selects the analysis tier (``primary_main`` / ``primary_extended``
    / ``all``); ``None`` uses the module default (Lead-switchable).  Raises
    :class:`MissingTrackBDataError` with the exact missing columns when the
    B-line predictions (or individual fields) are not available.
    """
    tier = subset or DEFAULT_ANALYSIS_SUBSET
    manifest = attach_pose_ids(manifest)
    analysis = select_analysis_set(manifest, tier)
    audit = check_trackb_inputs(manifest, predictions)
    if not audit["ready"]:
        detail = "missing columns: " + ", ".join(map(str, audit["missing"]))
        if audit.get("empty_required"):
            detail += "; present but empty: " + ", ".join(map(str, audit["empty_required"]))
        raise MissingTrackBDataError("Track B predictions are not ready; " + detail)
    preds = predictions.copy()
    # pose_key is the only required join key: target|pdb_id|ligand_id|method
    parts = preds["pose_key"].astype(str).str.split("|", expand=True)
    if parts.shape[1] != 4:
        raise ValueError("pose_key must be f'{target}|{pdb_id}|{ligand_id}|{method}'")
    preds["system_id"] = parts[0].astype(str)
    preds["pdb_id_pred"] = parts[1].astype(str)
    preds["ligand_id_pred"] = parts[2].astype(str)
    preds["pair_id"] = parts[0].astype(str) + "|" + parts[1].astype(str) + "|" + parts[2].astype(str)
    preds["pose_id"] = preds["pose_key"]
    if "cluster_id" not in preds.columns:
        preds["cluster_id"] = preds["pdb_id_pred"]
    # M sub-score: resolve the frame-correct membrane plane.  Plane frames and
    # pose frames must match (C-line convention: x_opm = R x_rcsb + t); the
    # dual-frame TrackB map exposes both pdb_id|opm and pdb_id|rcsb variants.
    provider = membrane_map
    if provider is None:
        provider = load_trackb_membrane_map()
    if provider is not None and hasattr(provider, "metadata"):
        preds["membrane_key"] = resolve_membrane_keys(preds, provider)
    else:
        preds["membrane_key"] = preds["pdb_id_pred"]
    keep = set(analysis["pair_id"].astype(str))
    n_before = len(preds)
    preds = preds[preds["pair_id"].astype(str).isin(keep)].reset_index(drop=True)
    if len(preds) == 0:
        raise MissingTrackBDataError(
            f"no predictions left after subset={subset!r} filtering "
            f"({n_before} input rows; expected pose_key pairs from the manifest)"
        )
    # primary/sensitivity labels come from the C-line manifest, not from user input
    lookup = analysis.set_index("pair_id")[["ligand_class_norm", "cluster_id", "target"]]
    features = build_features(preds, config=config, membrane_map=provider)
    features = features.merge(
        lookup.rename(columns={"ligand_class_norm": "ligand_class", "target": "manifest_target"}),
        left_on="pair_id", right_index=True, how="left",
    )
    features["analysis_set"] = tier
    return features
