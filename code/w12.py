"""W12 calibration gate (task-36).

Inputs (B-line task-37): per-pose tables under
``results/maes/w12/poses/<method>/poses.csv`` plus an enriched
``manifest_features.csv`` (or primitive MAES fields in the pose tables).

Pipeline: V/M/S (or precomputed) -> grouped CV + Platt/isotonic gate ->
evaluation labels (RMSD<1 primary; PLIF>=0.8 secondary; calibration label
RMSD<2 & PB) -> wild cluster bootstrap AUROC (system/CF clusters) ->
target/method/receptor-source strata -> Go/No-go (AUROC>=0.65 and |rho|>=0.40).

The dry-run mode explicitly marks features as partial (M/S neutral) and never
publishes calibration numbers; the formal mode refuses to run when the minimal
fields are missing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from .api import build_features, load_trackb_membrane_map
from .calibration import MAESModel, build_evaluation_labels
from .config import default_config
from .evaluation import (
    auroc,
    brier_score,
    ece,
    evaluate_predictions,
    spearman_rho,
    wild_cluster_bootstrap_auc,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
W12 = PROJECT_ROOT / "results" / "maes" / "w12"
POSE_DIR = W12 / "poses"
METHODS = ("vina", "vinardo", "boltz2")
FEATURES_MANIFEST = W12 / "manifest_features.csv"
REQUIRED = ["pose_key", "method", "system", "target", "pdb_id", "ligand_id", "pose_path", "rmsd", "posebusters_pass"]
PRIMITIVE_OR_PRECOMPUTED = ["V", "M", "S"]
M_PRIMITIVES = ["z_lig_mean", "z_ref_mean", "angle_lig_deg", "angle_ref_deg"]


def _read_pose_tables() -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    for method in METHODS:
        path = POSE_DIR / method / "poses.csv"
        if not path.exists():
            continue
        frame = pd.read_csv(path)
        frame["method"] = method
        frames.append(frame)
    extra = POSE_DIR / "afdb_vina" / "poses.csv"
    if extra.exists():
        frame = pd.read_csv(extra)
        frame["method"] = "vina_afdb"
        frames.append(frame)
    if not frames:
        raise FileNotFoundError(f"no W12 pose manifests under {POSE_DIR}")
    poses = pd.concat(frames, ignore_index=True)
    poses["pose_key"] = poses["pose_key"].astype(str)
    return poses


def _sanitize_phi(value: Any) -> Any:
    """Drop JSON nulls from phi dicts; subscores treat absent keys as NaN."""
    if not isinstance(value, str):
        return value
    try:
        payload = json.loads(value)
    except Exception:
        return value
    if isinstance(payload, dict):
        payload = {key: item for key, item in payload.items() if item is not None}
    return json.dumps(payload)


def _merge_features(poses: pd.DataFrame, features_manifest: Optional[Path] = None) -> pd.DataFrame:
    # B's per-method pose tables already carry the feature columns; only merge a
    # separate manifest when the features are genuinely absent (avoid row blow-up
    # from repeated pose_key rows before the duplicate check).
    if {"z_lig_mean", "z_ref_mean", "phi_lig", "phi_source"}.issubset(poses.columns):
        return poses
    if features_manifest is None or not Path(features_manifest).exists():
        features_manifest = FEATURES_MANIFEST
    if not Path(features_manifest).exists():
        return poses
    features = pd.read_csv(features_manifest)
    keys = ["pose_key"] + (["method"] if "method" in features.columns and "method" in poses.columns else [])
    if not set(keys).issubset(features.columns):
        return poses
    return poses.merge(features, on=keys, how="left", suffixes=("", "_feat"))


def _prepare(poses: pd.DataFrame) -> tuple[pd.DataFrame, Dict[str, Any]]:
    missing = [column for column in REQUIRED if column not in poses.columns]
    if missing:
        raise KeyError(f"W12 pose tables missing required columns: {missing}")
    frame = poses.copy()
    frame["system_id"] = frame["system"].astype(str)
    frame["method"] = frame["method"].astype(str)
    if "phi_source" in frame.columns:
        frame["phi_source_raw"] = frame["phi_source"].astype(str)
        frame["phi_source"] = frame["phi_source"].astype(str).str.lower().map(
            lambda value: "reference_structure" if "reference" in value else ("pose" if "pose" in value else "unknown")
        )
    for column in ("phi_lig", "phi_ref"):
        if column in frame.columns:
            frame[column] = frame[column].map(_sanitize_phi)
    if "status" in frame.columns:
        frame = frame[frame["status"].astype(str).str.lower().str.startswith("ok")].copy()
    frame["rmsd"] = pd.to_numeric(frame["rmsd"], errors="coerce")
    frame["plif_recovery"] = pd.to_numeric(frame.get("plif_recovery"), errors="coerce")
    frame["posebusters_pass"] = frame["posebusters_pass"].map(
        lambda value: str(value).strip().lower() in {"true", "1", "1.0", "yes"}
    )
    precomputed = all(column in frame.columns and frame[column].notna().any() for column in PRIMITIVE_OR_PRECOMPUTED)
    has_m = all(column in frame.columns and frame[column].notna().any() for column in M_PRIMITIVES)
    has_s = all(column in frame.columns and frame[column].notna().any() for column in ("phi_lig", "phi_ref"))
    has_v = (
        ("V" in frame.columns and frame["V"].notna().any())
        or ("pb_pass_rate" in frame.columns and frame["pb_pass_rate"].notna().any())
        or frame["posebusters_pass"].notna().any()
    )
    provider = load_trackb_membrane_map()
    complete = bool(precomputed or (has_m and has_s and has_v and provider is not None))
    if complete:
        frame["system_id"] = frame["system"].astype(str)
        if "pb_pass_rate" not in frame.columns and "V" not in frame.columns:
            frame["pb_pass_rate"] = frame["posebusters_pass"].astype(float)
        features = build_features(frame, config=default_config(), membrane_map=provider)
        partial = False
        missing_features = []
    else:
        # dry-run fallback: V from PoseBusters, M/S neutral; never used for formal claims
        if "pb_pass_rate" not in frame.columns:
            frame["pb_pass_rate"] = frame["posebusters_pass"].astype(float)
        features = build_features(frame, config=default_config(), membrane_map=provider)
        partial = True
        missing_features = []
        if not precomputed:
            if not has_m:
                missing_features.append("M primitives (z_lig_mean/z_ref_mean/angle_* + membrane_key)")
            if not has_s:
                missing_features.append("S primitives (phi_lig/phi_ref) or precomputed S")
        if provider is None:
            missing_features.append("Track B membrane map for M")
    # carry labels/identity into the feature frame
    for column in ("rmsd", "plif_recovery", "posebusters_pass", "method", "system", "target",
                   "pdb_id", "ligand_id", "pose_key", "receptor_source", "coord_frame"):
        if column in frame.columns:
            features[column] = frame[column].to_numpy()
    return features, {
        "n_rows": int(len(features)),
        "n_unique_pose_keys": int(features["pose_key"].nunique()),
        "methods_available": sorted(features["method"].astype(str).unique().tolist()),
        "n_methods": int(features["method"].nunique()),
        "methods_note": "Vinardo is a second scoring function on the same Vina poses; method diversity requires Boltz-2/other pose generators",
        "features_partial": bool(partial),
        "missing_features": missing_features,
    }


def _group_cv(
    features: pd.DataFrame,
    groups: np.ndarray,
    n_splits: int = 5,
    weights: Optional[Dict[str, float]] = None,
) -> Dict[str, Any]:
    from sklearn.model_selection import GroupKFold

    cfg = default_config()
    if weights is not None:
        cfg = replace(cfg, scores=replace(cfg.scores, weights=dict(weights)))
    unique_groups = np.unique(groups)
    n_splits = int(min(n_splits, len(unique_groups)))
    if n_splits < 2:
        return {"predictions": np.full(len(features), np.nan), "n_splits": 0,
                "note": "too few groups for grouped CV"}
    predictions = np.full(len(features), np.nan)
    fold_rows = []
    gkf = GroupKFold(n_splits=n_splits)
    for fold, (train_idx, test_idx) in enumerate(gkf.split(features, groups=groups)):
        train = features.iloc[train_idx]
        test = features.iloc[test_idx]
        try:
            model = MAESModel(config=cfg).fit(train)
            predictions[test_idx] = model.predict(test)
            fold_rows.append({"fold": fold, "n_train": int(len(train)), "n_test": int(len(test)),
                              "fit_ok": True})
        except Exception as error:  # noqa: BLE001 - recorded in report
            fold_rows.append({"fold": fold, "n_train": int(len(train)), "n_test": int(len(test)),
                              "fit_ok": False, "error": str(error)[:200]})
    return {"predictions": predictions, "n_splits": n_splits, "folds": fold_rows}


def _cluster_bootstrap_rho(
    predictions: np.ndarray,
    rmsd: np.ndarray,
    clusters: np.ndarray,
    n_boot: int = 2000,
    seed: int = 42,
) -> Dict[str, Any]:
    valid = np.isfinite(predictions) & np.isfinite(rmsd)
    p, r, c = predictions[valid], rmsd[valid], clusters[valid]
    if len(p) < 10:
        return {"rho": float("nan"), "ci_lo": float("nan"), "ci_hi": float("nan"), "n_boot": 0}
    unique = np.unique(c)
    groups = [np.where(c == value)[0] for value in unique]
    rng = np.random.default_rng(int(seed))
    draws = []
    for _ in range(int(n_boot)):
        idx = np.concatenate([groups[j] for j in rng.integers(0, len(unique), len(unique))])
        value = spearman_rho(p[idx], r[idx])
        if np.isfinite(value):
            draws.append(value)
    if not draws:
        return {"rho": float(spearman_rho(p, r)), "ci_lo": float("nan"), "ci_hi": float("nan"), "n_boot": 0}
    return {
        "rho": float(spearman_rho(p, r)),
        "ci_lo": float(np.percentile(draws, 2.5)),
        "ci_hi": float(np.percentile(draws, 97.5)),
        "n_boot": int(len(draws)),
        "method": "cluster bootstrap over systems/CF",
    }


def _reliability_bins(predictions: np.ndarray, labels: np.ndarray, n_bins: int = 5) -> List[Dict[str, Any]]:
    valid = np.isfinite(predictions)
    pred, label = predictions[valid], labels[valid].astype(int)
    if len(pred) < n_bins:
        return []
    ranks = pd.Series(pred).rank(method="first")
    bin_ids = pd.qcut(ranks, q=min(n_bins, len(np.unique(ranks))), labels=False, duplicates="drop")
    rows = []
    for bin_id in sorted(pd.Series(bin_ids).dropna().unique()):
        mask = bin_ids == bin_id
        rows.append({
            "bin": int(bin_id),
            "mean_pred": float(np.mean(pred[mask])),
            "observed": float(np.mean(label[mask])),
            "n": int(mask.sum()),
        })
    return rows


def _stratified(features: pd.DataFrame, predictions: np.ndarray, label: np.ndarray, column: str) -> List[Dict[str, Any]]:
    rows = []
    frame = features.copy()
    frame["_pred"] = predictions
    frame["_label"] = label
    for value, group in frame.groupby(column, dropna=False):
        valid = np.isfinite(group["_pred"].to_numpy(dtype=float))
        y = group["_label"].to_numpy(dtype=int)[valid]
        s = group["_pred"].to_numpy(dtype=float)[valid]
        if len(y) < 4 or y.sum() == 0 or y.sum() == len(y):
            rows.append({"stratum": str(value), "n": int(len(y)), "auroc": float("nan"), "note": "too few/one-class"})
            continue
        rows.append({"stratum": str(value), "n": int(len(y)), "auroc": auroc(y, s),
                     "rho_rmsd": spearman_rho(s, group["rmsd"].to_numpy(dtype=float)[valid]) if "rmsd" in group else float("nan")})
    return rows


def _figure(report: Dict[str, Any], tag: str) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = PROJECT_ROOT / "figures" / f"fig6_w12_calibration_{tag}.png"
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8))
    axes[0].plot([0, 1], [0, 1], ls="--", color="0.6")
    bins = report.get("primary", {}).get("reliability_bins", [])
    if bins:
        axes[0].plot([row["mean_pred"] for row in bins], [row["observed"] for row in bins], "o-", color="#225ea8")
        for row in bins:
            axes[0].annotate(f"n={row['n']}", (row["mean_pred"], row["observed"]), fontsize=6, xytext=(2, 2), textcoords="offset points")
    axes[0].set_xlim(0, 1)
    axes[0].set_ylim(0, 1)
    axes[0].set_xlabel("Predicted MAES (primary V+M)")
    axes[0].set_ylabel("Observed acceptable-pose rate")
    axes[0].set_title(f"W12 calibration ({tag}; n={report.get('n')})")
    metrics = report.get("metrics", {})
    axes[1].axis("off")
    lines = [
        f"AUROC: {metrics.get('auroc_primary', float('nan')):.3f}",
        f"|rho|: {abs(metrics.get('rho_maes_rmsd', float('nan'))):.3f}",
        f"Brier: {metrics.get('brier_primary', float('nan')):.3f}",
        f"ECE: {metrics.get('ece_primary', float('nan')):.3f}",
        f"Calibrator: {report.get('calibrator', '?')}",
        f"Verdict: {report.get('verdict', '?')}",
    ]
    axes[1].text(0.02, 0.95, "\n".join(lines), va="top", fontsize=9, family="monospace")
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out


def run(tag: str = "dryrun", n_boot: int = 4999) -> Dict[str, Any]:
    poses = _read_pose_tables()
    poses = _merge_features(poses)
    features, readiness = _prepare(poses)
    if readiness["features_partial"] and tag != "dryrun":
        raise RuntimeError(
            "W12 features are partial (M/S neutral). Provide manifest_features.csv with V/M/S or "
            "z/angle/phi primitives, then rerun with tag='formal'."
        )
    manifest_path = W12 / "manifest_features.csv"
    manifest_sha256 = (
        hashlib.sha256(manifest_path.read_bytes()).hexdigest() if manifest_path.exists() else None
    )
    afdb_excluded_by_value = (
        features["excluded_predicted"].astype(str).value_counts().to_dict()
        if "excluded_predicted" in features.columns else {}
    )
    if tag == "formal":
        if "phi_source" not in features.columns:
            raise RuntimeError(
                "formal W12 requires row-level phi_source for every method x pose; "
                "ask B task-37 to add phi_source (reference_structure|pose)."
            )
        duplicates = int(features.duplicated(["pose_key", "method"]).sum())
        if duplicates:
            raise RuntimeError(
                f"formal W12 found {duplicates} duplicate (pose_key, method) rows; "
                "Boltz-2 diffusion samples need unique sample_id/pose_key suffixes."
            )
    features["system_id"] = features["system"].astype(str)
    primary_features = features.copy()
    # Primary (frozen): V+M only; S is neutral for reference-structure phi.
    if "S" in primary_features.columns:
        primary_features["S"] = 1.0
    if "s_neutral" in primary_features.columns:
        primary_features["s_neutral"] = True
    primary_weights = {"V": 0.4999995, "M": 0.4999995, "S": 1e-6}

    def analyze(frame: pd.DataFrame, weights: Optional[Dict[str, float]], label: str) -> Dict[str, Any]:
        groups_here = frame["system"].astype(str).to_numpy()
        cv_here = _group_cv(frame, groups_here, weights=weights)
        predictions_here = cv_here["predictions"]
        valid_here = np.isfinite(predictions_here)
        eval_label_here = build_evaluation_labels(frame, default_config().labels, which="primary")
        metrics_here: Dict[str, Any] = {}
        if valid_here.sum() >= 10 and eval_label_here[valid_here].sum() >= 2:
            metrics_here = {
                "auroc_primary": float(auroc(eval_label_here[valid_here], predictions_here[valid_here])),
                "brier_primary": float(brier_score(eval_label_here[valid_here], predictions_here[valid_here])),
                "ece_primary": float(ece(eval_label_here[valid_here], predictions_here[valid_here])),
                "rho_maes_rmsd": float(spearman_rho(
                    predictions_here[valid_here],
                    frame.loc[valid_here, "rmsd"].to_numpy(dtype=float),
                )),
                "n_eval": int(valid_here.sum()),
                "n_eval_positive": int(eval_label_here[valid_here].sum()),
            }
        cfg_here = default_config()
        if weights is not None:
            cfg_here = replace(cfg_here, scores=replace(cfg_here.scores, weights=dict(weights)))
        model_here = MAESModel(config=cfg_here).fit(frame)
        wild_here = wild_cluster_bootstrap_auc(
            eval_label_here[valid_here], predictions_here[valid_here],
            groups_here[valid_here], n_boot=n_boot,
        ) if valid_here.sum() >= 10 else {}
        rho_here = _cluster_bootstrap_rho(
            predictions_here,
            frame["rmsd"].to_numpy(dtype=float),
            groups_here,
            n_boot=min(2000, n_boot),
        ) if "rmsd" in frame.columns else {}
        if rho_here and np.isfinite(rho_here.get("ci_lo", float("nan"))):
            metrics_here["rho_maes_rmsd_ci95"] = [rho_here["ci_lo"], rho_here["ci_hi"]]
        reliability_here = _reliability_bins(predictions_here, eval_label_here)
        strata_here = {
            "method": _stratified(frame, predictions_here, eval_label_here, "method"),
            "target": _stratified(frame, predictions_here, eval_label_here, "target"),
            "receptor_source": _stratified(frame, predictions_here, eval_label_here, "receptor_source")
            if "receptor_source" in frame else [],
        }
        return {
            "label": label,
            "metrics": metrics_here,
            "wild_cluster_bootstrap": wild_here,
            "rho_cluster_bootstrap": rho_here,
            "reliability_bins": reliability_here,
            "strata": strata_here,
            "calibrator": dict(model_here.metadata).get("calibrator"),
            "calibration_metadata": dict(model_here.metadata),
            "cv": {"n_splits": cv_here.get("n_splits"), "folds": cv_here.get("folds", [])},
            "predictions": predictions_here,
        }

    primary = analyze(primary_features, primary_weights, "V+M (primary; S neutral for reference phi)")
    exploratory = analyze(features, None, "V+M+S (exploratory; S axis limited)")
    predictions = primary["predictions"]
    metrics = primary["metrics"]
    n_pose_observations = int(len(features))
    n_unique_pose_keys = int(features["pose_key"].nunique())
    n_unique_pose_key_method = int(features[["pose_key", "method"]].drop_duplicates().shape[0])
    n_duplicates = int(features.duplicated(["pose_key", "method"]).sum())
    phi_source_by_method = (
        features.groupby(["method", "phi_source"]).size().unstack(fill_value=0).to_dict("index")
        if "phi_source" in features.columns else {}
    )
    go = bool(
        tag == "formal"
        and readiness["features_partial"] is False
        and n_pose_observations >= 100
        and readiness["n_methods"] >= 3
        and metrics.get("auroc_primary", 0) >= 0.65
        and abs(metrics.get("rho_maes_rmsd", 0)) >= 0.40
    )
    verdict = "GO" if go else ("NO_GO" if tag == "formal" and not readiness["features_partial"] else "DRY_RUN_PARTIAL")
    report = {
        "task": "W12 calibration gate (task-36)",
        "tag": tag,
        "dry_run": tag != "formal",
        "analysis_plan": "results/maes/w12/W12_ANALYSIS_PLAN.md",
        "readiness": readiness,
        "manifest_sha256": manifest_sha256,
        "afdb_excluded_by_value": afdb_excluded_by_value,
        "afdb_usable_rows": int(afdb_excluded_by_value.get("False", 0)),
        "afdb_excluded_rows": int(afdb_excluded_by_value.get("True", 0)),
        "n": n_pose_observations,
        "n_pose_observations": n_pose_observations,
        "n_unique_pose_keys": n_unique_pose_keys,
        "n_unique_pose_key_method": n_unique_pose_key_method,
        "n_duplicates": n_duplicates,
        "n_methods": readiness["n_methods"],
        "methods": readiness["methods_available"],
        "phi_source_by_method": phi_source_by_method,
        "features_partial": readiness["features_partial"],
        "missing_features": readiness["missing_features"],
        "primary": {key: value for key, value in primary.items() if key != "predictions"},
        "exploratory": {key: value for key, value in exploratory.items() if key != "predictions"},
        "cv": primary["cv"],
        "metrics": primary["metrics"],
        "wild_cluster_bootstrap": primary["wild_cluster_bootstrap"],
        "strata": primary["strata"],
        "calibrator": primary["calibrator"],
        "calibration_metadata": primary["calibration_metadata"],
        "go_criteria": {
            "auroc_ge_0_65": primary["metrics"].get("auroc_primary", 0) >= 0.65,
            "abs_rho_ge_0_40": abs(primary["metrics"].get("rho_maes_rmsd", 0)) >= 0.40,
            "n_pose_observations_ge_100": n_pose_observations >= 100,
            "n_methods_ge_3": readiness["n_methods"] >= 3,
            "features_complete": not readiness["features_partial"],
        },
        "verdict": verdict,
        "go": go,
        "no_go_plan": (
            "If No-go: MAES is reported descriptively only (no Go claim); calibration/reliability stays in SI; "
            "MAES is not used for method ranking. Q2 fallback: keep the two-regime benchmark + per-target "
            "AUC/EF5/CI and failure-mode analyses as the primary contribution; MAES appears as an exploratory "
            "pose-quality score with explicit limitations."
        ),
        "caveats": [
            "dry-run/partial features: M/S neutral until B provides MAES primitives; do not cite these numbers",
            "primary analysis is V+M; reference-structure S is set neutral (s_neutral=True) and S is exploratory only",
            "Vinardo is a second scoring function on the same Vina poses; method diversity requires Boltz-2/other engines",
            "grouped CV by system; wild cluster bootstrap by system/CF",
        ] if tag != "formal" else [
            "primary analysis is V+M; reference-structure S is neutral and S is exploratory only",
            "Vinardo is a second scoring function on the same Vina poses; method diversity requires Boltz-2/other engines",
            "grouped CV by system; wild cluster bootstrap by system/CF",
        ],
        "generated_at": pd.Timestamp.now().isoformat(timespec="seconds"),
    }
    W12.mkdir(parents=True, exist_ok=True)
    report_path = W12 / f"w12_report_{tag}.json"
    metrics_path = W12 / f"w12_metrics_{tag}.csv"
    predictions_path = W12 / f"w12_predictions_{tag}.csv"
    pd.DataFrame([{
        "n": report["n"], **primary["metrics"],
        "verdict": verdict, "features_partial": readiness["features_partial"],
    }]).to_csv(metrics_path, index=False)
    pd.DataFrame({
        "pose_key": features["pose_key"].astype(str),
        "method": features["method"].astype(str),
        "system": features["system"].astype(str),
        "rmsd": features["rmsd"],
        "MAES_primary_VM": primary["predictions"],
        "MAES_exploratory_VMS": exploratory["predictions"],
    }).to_csv(predictions_path, index=False)
    report["report_json"] = str(report_path)
    report["metrics_csv"] = str(metrics_path)
    report["predictions_csv"] = str(predictions_path)
    report["figure_png"] = str(_figure(report, tag))
    report_path.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", default="dryrun", choices=["dryrun", "formal"])
    parser.add_argument("--n-boot", type=int, default=4999)
    args = parser.parse_args()
    out = run(tag=args.tag, n_boot=args.n_boot)
    print(json.dumps({k: out[k] for k in ("tag", "n", "n_unique_pose_keys", "n_methods", "features_partial", "verdict", "metrics", "report_json", "figure_png")}, indent=2, default=str))
