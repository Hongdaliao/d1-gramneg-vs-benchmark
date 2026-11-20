"""Track B prediction-table ingestion report (first real-data validation).

When B delivers ``data/trackB/vina_predictions.csv`` this module produces, in one
call: the input audit (``check_trackb_inputs``), failure classification, metric
distributions, M-term coverage (``m_neutral`` / ``m_frame_mismatch`` /
``coord_frame``), the first MAES calibration on real poses (AUROC / Brier / ECE /
rho + ablation), the ICC-corrected n_eff with the final
target x CF x sub-site cluster definition, and the beta-barrel plane-choice
sensitivity.  A 5-line summary is printed first (Lead request) before the full
JSON is written.

Usage::

    python -m maes ingest --predictions data/trackB/vina_predictions.csv
    python -c "from maes.ingest import run_ingest_report; run_ingest_report()"
"""

from __future__ import annotations

import copy
import gzip
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from .api import (
    attach_pose_ids,
    build_features_from_manifest,
    check_trackb_inputs,
    cluster_keys,
    load_pose_manifest,
    select_analysis_set,
)
from .calibration import MAESModel, build_evaluation_labels
from .combine import entropy_weights
from .config import default_config
from .evaluation import (
    ablation_study,
    design_effect_n_eff,
    evaluate_predictions,
    wild_cluster_bootstrap_auc,
)
from .membrane_map import MembranePlane
from .subscores import insertion_energy, membrane_consistency, principal_axis_angle_deg

__all__ = [
    "beta_barrel_m_sensitivity",
    "distribution_summary",
    "failure_classification",
    "read_structure_atoms",
    "run_ingest_report",
]

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = PROJECT_ROOT / "data" / "trackB" / "pose_manifest_curated.csv"
DEFAULT_PREDICTIONS = PROJECT_ROOT / "results" / "trackB" / "vina_predictions.csv"
DEFAULT_BETA_PLANES = PROJECT_ROOT / "data" / "trackB" / "beta_barrel_planes.csv"
DEFAULT_STRUCT_DIR = PROJECT_ROOT / "data" / "trackB" / "cache" / "struct"
DEFAULT_PB_CHECKS = PROJECT_ROOT / "results" / "trackB" / "pb_checks.csv"
FIGURES_DIR = PROJECT_ROOT / "figures"


# ------------------------------------------------------------------ structures
def _read_lines(path: Path) -> List[str]:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", errors="ignore") as fh:
            return fh.read().splitlines()
    return path.read_text(encoding="utf-8", errors="ignore").splitlines()


def read_structure_atoms(pdb_id: str, struct_dir: Path = DEFAULT_STRUCT_DIR) -> Optional[pd.DataFrame]:
    """Parse ATOM/HETATM records from a PDB or mmCIF file (gzip aware)."""
    candidates = [struct_dir / f"{pdb_id}{suffix}" for suffix in
                  (".pdb", ".pdb.gz", ".cif", ".cif.gz", ".mmcif", ".mmcif.gz", ".opm.pdb")]
    path = next((p for p in candidates if p.exists()), None)
    if path is None:
        return None
    lines = _read_lines(path)
    # NB: for gzip files Path.suffix is '.gz', so inspect all suffixes
    is_cif = any(s in {".cif", ".mmcif"} for s in (x.lower() for x in path.suffixes))
    if not is_cif:
        rows = []
        for line in lines:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            rows.append(
                {
                    "group": line[:6].strip(),
                    "atom": line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "chain": line[21:22].strip(),
                    "altloc": line[16:17].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "element": (line[76:78].strip() or line[12:14].strip()),
                }
            )
        return pd.DataFrame(rows)
    # mmCIF: locate the _atom_site loop and map columns
    col_names: List[str] = []
    rows = []
    in_loop = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("_atom_site."):
            col_names.append(stripped.split(".", 1)[1].strip())
            in_loop = True
            continue
        if in_loop:
            if not stripped or stripped.startswith("#") or stripped.startswith("_") or stripped == "loop_":
                in_loop = False
                continue
            tokens = stripped.split()
            if len(tokens) != len(col_names):
                continue
            rec = dict(zip(col_names, tokens))
            group = str(rec.get("group_PDB", "ATOM")).strip()
            if group not in {"ATOM", "HETATM"}:
                continue
            element = str(rec.get("type_symbol", "C")).strip()
            rows.append(
                {
                    "group": group,
                    "atom": str(rec.get("label_atom_id", "")).strip(),
                    "resname": str(rec.get("label_comp_id", "")).strip(),
                    "chain": str(rec.get("label_asym_id", rec.get("auth_asym_id", ""))).strip(),
                    "auth_chain": str(rec.get("auth_asym_id", "")).strip(),
                    "altloc": str(rec.get("label_alt_id", ".")).strip(),
                    "x": float(rec["Cartn_x"]),
                    "y": float(rec["Cartn_y"]),
                    "z": float(rec["Cartn_z"]),
                    "element": element,
                }
            )
    return pd.DataFrame(rows) if rows else None


def _extract_ligand_atoms(
    pdb_id: str,
    ligand_id: str,
    chain: Optional[str] = None,
    struct_dir: Path = DEFAULT_STRUCT_DIR,
):
    atoms = read_structure_atoms(pdb_id, struct_dir)
    if atoms is None or len(atoms) == 0:
        return None
    subset = atoms[(atoms["group"] == "HETATM") & (atoms["resname"] == str(ligand_id))]
    if "auth_chain" in subset.columns:
        chain_col = "auth_chain"
    else:
        chain_col = "chain"
    if chain is not None and str(chain).strip() and chain_col in subset.columns:
        subset = subset[subset[chain_col].astype(str) == str(chain).strip()]
    if len(subset) == 0 and chain is not None:
        # never pool multiple copies of the same ligand: fall back to the first copy
        all_lig = atoms[(atoms["group"] == "HETATM") & (atoms["resname"] == str(ligand_id))]
        if chain_col in all_lig.columns:
            first_chains = all_lig[chain_col].astype(str).dropna().unique()
            if len(first_chains):
                subset = all_lig[all_lig[chain_col].astype(str) == first_chains[0]]
        else:
            subset = all_lig
    if len(subset) < 3:
        return None
    coords = subset[["x", "y", "z"]].to_numpy(dtype=float)
    elements = [str(e).strip().capitalize() for e in subset["element"]]
    return coords, elements


# ------------------------------------------------------------------ beta barrel
def beta_barrel_m_sensitivity(
    planes_csv: str | Path = DEFAULT_BETA_PLANES,
    struct_dir: str | Path = DEFAULT_STRUCT_DIR,
    cfg=None,
) -> Dict[str, Any]:
    """Score the 3 beta-barrel reference poses with the geometry vs OPM plane.

    Both planes are expressed in the RCSB original frame; the same ligand
    coordinates are used, so the M difference is driven purely by the plane
    definition.  ``M_cross_*`` treats the other plane's pose attributes as the
    reference, i.e. how much M would change if the reference was defined with
    the other plane.
    """
    planes_csv = Path(planes_csv)
    if not planes_csv.exists():
        return {"available": False, "reason": "beta_barrel_planes.csv missing"}
    cfg = cfg or default_config()
    table = pd.read_csv(planes_csv)
    rows: List[Dict[str, Any]] = []
    for _, r in table.iterrows():
        extracted = _extract_ligand_atoms(str(r["pdb_id"]), str(r["ligand_id"]), r.get("chain"), Path(struct_dir))
        if extracted is None:
            rows.append({"pdb_id": r["pdb_id"], "ligand_id": r["ligand_id"], "status": "no coordinates"})
            continue
        coords, elements = extracted
        geo = MembranePlane(
            normal=[r["geometry_normal_x"], r["geometry_normal_y"], r["geometry_normal_z"]],
            center=[r["geometry_center_x"], r["geometry_center_y"], r["geometry_center_z"]],
            thickness=float(r["geometry_thickness_A"]), source="geometric",
            confidence="medium", coord_frame="rcsb",
        )
        if not bool(r.get("opm_transfer_available", False)):
            # geometry plane only: report attributes and flag the missing contrast
            zg = geo.project(coords)
            ag = principal_axis_angle_deg(coords, geo.normal_array)
            dGg = insertion_energy(
                zg, elements, geo.half_thickness,
                interface_width=cfg.scores.interface_width_angstrom,
            )
            rows.append(
                {
                    "target": r["target"], "pdb_id": r["pdb_id"], "ligand_id": r["ligand_id"],
                    "status": "geometry_only",
                    "mean_z_geometry": float(np.mean(zg)),
                    "angle_geometry_deg": float(ag),
                    "dG_geometry": float(dGg),
                    "note": "no OPM transfer for this entry; geometry plane only",
                }
            )
            continue
        opm = MembranePlane(
            normal=[r["opm_normal_x"], r["opm_normal_y"], r["opm_normal_z"]],
            center=[r["opm_center_x"], r["opm_center_y"], r["opm_center_z"]],
            thickness=float(r["opm_thickness_A"]), source="opm_transfer",
            confidence=str(r.get("opm_transfer_confidence", "low")), coord_frame="rcsb",
        )

        def attributes(plane):
            z = plane.project(coords)
            angle = principal_axis_angle_deg(coords, plane.normal_array)
            dg = insertion_energy(
                z, elements, plane.half_thickness,
                interface_width=cfg.scores.interface_width_angstrom,
            )
            return z, angle, dg, coords.mean(axis=0)

        zg, ag, dGg, cg = attributes(geo)
        zo, ao, dGo, co = attributes(opm)
        res_g = membrane_consistency(
            cfg.scores, z_lig=zg, z_ref=zo, plane=geo,
            angle_lig_deg=ag, angle_ref_deg=ao, dG_lig=dGg, dG_ref=dGo,
            centroid_lig=cg, centroid_ref=cg,
        )
        res_o = membrane_consistency(
            cfg.scores, z_lig=zo, z_ref=zg, plane=opm,
            angle_lig_deg=ao, angle_ref_deg=ag, dG_lig=dGo, dG_ref=dGg,
            centroid_lig=co, centroid_ref=co,
        )
        rows.append(
            {
                "target": r["target"], "pdb_id": r["pdb_id"], "ligand_id": r["ligand_id"],
                "status": "ok",
                "mean_z_geometry": float(np.mean(zg)), "mean_z_opm": float(np.mean(zo)),
                "delta_mean_z": float(abs(np.mean(zg) - np.mean(zo))),
                "angle_geometry_deg": float(ag), "angle_opm_deg": float(ao),
                "delta_angle_deg": float(abs(ag - ao)),
                "dG_geometry": float(dGg), "dG_opm": float(dGo), "delta_dG": float(abs(dGg - dGo)),
                "M_cross_geometry": float(res_g.M), "M_cross_opm": float(res_o.M),
                "M_cross_min": float(min(res_g.M, res_o.M)),
                "plane_normal_angle_deg": float(r.get("normal_angle_deg", np.nan)),
                "plane_center_distance_A": float(r.get("center_distance_A", np.nan)),
            }
        )
    frame = pd.DataFrame(rows)
    ok = frame[frame["status"] == "ok"] if len(frame) else frame
    geometry_only = frame[frame["status"] == "geometry_only"] if len(frame) else frame
    if len(ok) == 0:
        return {
            "available": True, "n_poses": 0, "n_with_both_planes": 0,
            "n_geometry_only": int(len(geometry_only)), "rows": rows,
            "recommendation": "no evaluable two-plane poses",
        }
    min_m = float(ok["M_cross_min"].min())
    recommendation = (
        "report the beta-barrel subset SEPARATELY (plane choice materially changes M; "
        "min cross-plane M = {:.3f})".format(min_m)
        if min_m < 0.90
        else "beta-barrel subset can stay pooled (plane choice has limited M effect)"
    )
    if len(geometry_only):
        recommendation += (
            "; {} pose(s) have geometry plane only (no OPM transfer) and must be "
            "reported with that caveat".format(len(geometry_only))
        )
    return {
        "available": True,
        "n_poses": int(len(ok)),
        "n_with_both_planes": int(len(ok)),
        "n_geometry_only": int(len(geometry_only)),
        "max_delta_mean_z_A": float(ok["delta_mean_z"].max()),
        "max_delta_angle_deg": float(ok["delta_angle_deg"].max()),
        "max_delta_dG": float(ok["delta_dG"].max()),
        "min_M_cross": min_m,
        "mean_M_cross": float(ok["M_cross_min"].mean()),
        "recommendation": recommendation,
        "rows": ok.round(4).to_dict(orient="records"),
        "geometry_only_rows": geometry_only.round(4).to_dict(orient="records") if len(geometry_only) else [],
    }


# ------------------------------------------------------------------- reporting
def normalise_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    """Coerce the B-line table's stringly-typed metric columns.

    The final Track B table encodes unavailable values as strings such as
    ``NA(no_reference_coords)`` / ``NA(posebusters_unavailable)``; pandas must
    see them as NaN, and ``posebusters_pass`` must become a real bool column
    (for now the strict bools are pending, so this can be all False).
    """
    out = predictions.copy()
    numeric_columns = (
        "rmsd", "plif_recovery", "pb_pass_rate", "score", "runtime_min",
        "z_lig_mean", "z_ref_mean", "angle_lig_deg", "angle_ref_deg",
        "dG_lig", "dG_ref", "strain_kcal",
    )
    for column in numeric_columns:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    if "posebusters_pass" in out.columns:
        mapping = {
            "true": True, "1": True, "1.0": True, "yes": True,
            "false": False, "0": False, "0.0": False, "no": False,
        }
        out["posebusters_pass"] = (
            out["posebusters_pass"].astype(str).str.strip().str.lower().map(mapping)
        )
    if "rmsd" in out.columns:
        raw_rmsd = predictions["rmsd"].astype(str)
        out["rmsd_reason"] = raw_rmsd.where(raw_rmsd.str.lower().str.startswith("na("), "")
    if "posebusters_pass" in out.columns:
        raw_pb = predictions["posebusters_pass"].astype(str)
        out["posebusters_reason"] = raw_pb.where(raw_pb.str.lower().str.startswith("na("), "")
    return out


def distribution_summary(values) -> Dict[str, Any]:
    """Quantile summary of a numeric column (NaNs are counted, not imputed)."""
    series = pd.to_numeric(pd.Series(values), errors="coerce")
    valid = series.dropna()
    out: Dict[str, Any] = {
        "n": int(len(series)),
        "n_missing": int(series.isna().sum()),
    }
    if len(valid):
        out.update(
            {
                "min": float(valid.min()),
                "q25": float(valid.quantile(0.25)),
                "median": float(valid.median()),
                "q75": float(valid.quantile(0.75)),
                "max": float(valid.max()),
                "mean": float(valid.mean()),
            }
        )
    return out


def failure_classification(predictions: pd.DataFrame) -> Dict[str, Any]:
    """Count run statuses and classify failure messages (OOM / timeout / parse / other)."""
    out: Dict[str, Any] = {"n_rows": int(len(predictions))}
    if "status" not in predictions.columns:
        out["note"] = "no status column; assume all rows attempted"
        out["status_counts"] = {"unknown": int(len(predictions))}
        return out
    status = predictions["status"].astype(str).str.strip().str.lower()
    out["status_counts"] = {str(k): int(v) for k, v in status.value_counts().items()}
    failed = predictions[status != "ok"]
    # metrics-only re-delivery: rows can be `ok` while rmsd/posebusters_pass are
    # still empty (B-line dock_pdb_parse stage not yet run) — report explicitly
    n_ok_with_missing_metrics: Dict[str, int] = {}
    for column in ("rmsd", "posebusters_pass", "plif_recovery"):
        if column in predictions.columns:
            n_ok_with_missing_metrics[column] = int(predictions.loc[status == "ok", column].isna().sum())
    out["n_ok_with_missing_metrics"] = n_ok_with_missing_metrics
    buckets: Dict[str, int] = {"oom": 0, "timeout": 0, "parse_error": 0, "missing_input": 0, "other": 0}
    for _, row in failed.iterrows():
        text = (str(row.get("error_message", "")) + " " + str(row.get("status", ""))).lower()
        if "memory" in text or "oom" in text or "out of memory" in text:
            buckets["oom"] += 1
        elif "timeout" in text or "timed out" in text or "time limit" in text:
            buckets["timeout"] += 1
        elif "parse" in text or "malformed" in text or "format" in text:
            buckets["parse_error"] += 1
        elif "missing" in text or "not found" in text or "no such" in text:
            buckets["missing_input"] += 1
        else:
            buckets["other"] += 1
    out["failure_buckets"] = buckets
    return out


# Final Track B table: 40 rows, 2 failed + 4 `rmsd=NA(no_reference_coords)` rows
# are retained for audit but excluded from calibration/evaluation.  Readiness
# therefore requires >= 30 successful rows WITH an RMSD (Lead 2026-11-16).
MIN_OK_WITH_RMSD = 30
MIN_ANALYSIS_ROWS = 20
# Amendment #6 (draft): operational PoseBusters criterion used if strict
# all-pass is not attainable on PDBQT-derived hydrogen conformations
PB_OPERATIONAL_THRESHOLD = 0.95
# Amendment #6 (v1.0 rewrite): PoseBusters redock includes `rmsd_≤_2å`, which is
# circular with the RMSD-based calibration label; exclude every RMSD-based check
# from the physics/chemistry PB definition.
RMSD_BASED_CHECK_PATTERNS = ("rmsd", "reference_distance", "distance_to_reference", "rmsd_to_ref")


def _to_bool_series(series: pd.Series) -> pd.Series:
    """Parse a boolean-ish column ('True'/'False'/1/0/yes/no/NA...) to object bool/NaN."""
    mapping = {
        "true": True, "1": True, "1.0": True, "yes": True, "pass": True, "t": True,
        "false": False, "0": False, "0.0": False, "no": False, "fail": False, "f": False,
    }
    return series.astype(str).str.strip().str.lower().map(mapping)


def physics_posebusters_pass(
    checks: pd.DataFrame,
    exclude_patterns: tuple = RMSD_BASED_CHECK_PATTERNS,
) -> Dict[str, Any]:
    """Strict all-pass over non-RMSD PoseBusters checks (per-check bool matrix).

    Drops every RMSD-based column (including the redock ``rmsd_≤_2å`` check)
    before computing the physics/chemistry strict bool, so the PB term cannot be
    circular with the ``RMSD < 2 Å`` component of the calibration label.
    """
    if checks is None or len(checks.columns) == 0:
        raise ValueError("empty PoseBusters check matrix")
    excluded = [
        column for column in checks.columns
        if any(pattern in str(column).lower() for pattern in exclude_patterns)
    ]
    kept = [column for column in checks.columns if column not in excluded]
    if not kept:
        raise ValueError("all PoseBusters checks were RMSD-based; no physics checks left")
    bools = checks[kept].apply(_to_bool_series)
    any_value = bools.notna().any(axis=1)
    all_true = bools.astype("boolean").fillna(False).astype(bool).all(axis=1)
    passed = all_true.astype(object).where(any_value, np.nan)
    n_true = int((passed == True).sum())  # noqa: E712
    n_false = int((passed == False).sum())  # noqa: E712
    n_missing = int(passed.isna().sum())
    non_missing = int(any_value.sum())
    return {
        "passed": passed,
        "n_checks": int(len(kept)),
        "excluded_checks": excluded,
        "physics_pass_rate": float(n_true / non_missing) if non_missing else float("nan"),
        "n_true": n_true,
        "n_false": n_false,
        "n_missing": n_missing,
        "n_failed": n_false,
    }


def _is_booleanish(series: pd.Series) -> bool:
    """True when every non-null value parses as a boolean (True/False/1/0/...)."""
    mapped = _to_bool_series(series)
    non_null = series.notna()
    if not bool(non_null.any()):
        return False
    return bool(mapped[non_null].notna().all())


PB_METADATA_COLUMNS = ("pose_key", "method", "status", "method_version", "note", "error_message")


def verify_pb_checks(predictions: pd.DataFrame, checks: pd.DataFrame) -> Dict[str, Any]:
    """Independently recompute the physics PB pass from B's per-check matrix.

    Drops metadata columns (``pose_key/method/status/...``) and every RMSD-based
    check (``rmsd_≤_2å``), computes the strict physics all-pass, then aligns it
    with B's ``posebusters_pass`` column by ``pose_key`` and reports
    agreement/mismatches (the matrix is the authority when they differ).
    """
    if "pose_key" not in checks.columns:
        raise KeyError("pb_checks.csv needs a 'pose_key' column")
    metadata_cols = [c for c in checks.columns if c in PB_METADATA_COLUMNS]
    candidate_cols = [c for c in checks.columns if c not in metadata_cols]
    check_cols = [c for c in candidate_cols if _is_booleanish(checks[c])]
    non_check_cols = [c for c in candidate_cols if c not in check_cols]
    if not check_cols:
        raise ValueError("pb_checks.csv has no boolean check columns")
    matrix = checks[check_cols]
    physics = physics_posebusters_pass(matrix)
    by_pose = pd.DataFrame(
        {
            "pose_key": checks["pose_key"].astype(str),
            "physics_pass": physics["passed"],
        }
    )
    by_pose = by_pose.drop_duplicates(subset=["pose_key"], keep="first")
    merged = predictions.merge(by_pose, on="pose_key", how="left")
    b_bool = (
        _to_bool_series(merged["posebusters_pass"])
        if "posebusters_pass" in merged.columns
        else pd.Series(np.nan, index=merged.index)
    )
    p_bool = merged["physics_pass"]
    both = b_bool.notna() & p_bool.notna()
    agreement = {
        "n_compared": int(both.sum()),
        "n_physics_true": physics["n_true"],
        "n_physics_false": physics["n_false"],
        "n_physics_missing": physics["n_missing"],
        "n_b_true": int((b_bool == True).sum()),  # noqa: E712
        "n_b_false": int((b_bool == False).sum()),  # noqa: E712
        "n_b_missing": int(b_bool.isna().sum()),
        "n_both_true": int(((b_bool == True) & (p_bool == True)).sum()),  # noqa: E712
        "n_both_false": int(((b_bool == False) & (p_bool == False)).sum()),  # noqa: E712
        "n_mismatch": int((both & (b_bool != p_bool)).sum()),
    }
    mismatch_rows = []
    for idx in np.where(both & (b_bool != p_bool))[0]:
        row = merged.iloc[idx]
        mismatch_rows.append(
            {
                "pose_key": str(row["pose_key"]),
                "b_posebusters_pass": bool(b_bool.loc[row.name]),
                "physics_pass": bool(p_bool.loc[row.name]),
                "rmsd": None if pd.isna(row.get("rmsd")) else float(row["rmsd"]),
            }
        )
    rmsd = pd.to_numeric(merged.get("rmsd", pd.Series(np.nan, index=merged.index)), errors="coerce")
    expected_positives = int(((rmsd < 2.0) & (p_bool == True)).sum())  # noqa: E712
    expected_analysis = int((rmsd.notna() & p_bool.notna()).sum())
    independent_ok = bool(
        agreement["n_mismatch"] == 0
        and agreement["n_physics_true"] == agreement["n_b_true"]
        and agreement["n_physics_false"] == agreement["n_b_false"]
        and agreement["n_physics_missing"] == agreement["n_b_missing"]
    )
    return {
        "pb_independent_ok": independent_ok,
        "n_checks_total": int(len(check_cols) + len(non_check_cols)),
        "n_checks_physics": physics["n_checks"],
        "excluded_checks": physics["excluded_checks"],
        "metadata_columns": metadata_cols,
        "non_boolean_columns": non_check_cols,
        "physics_pass_rate": physics["physics_pass_rate"],
        "agreement": agreement,
        "mismatch_rows": mismatch_rows,
        "expected_analysis_rows": expected_analysis,
        "expected_calibration_positives": expected_positives,
        "merged": merged,
    }


def exclusion_breakdown(predictions: pd.DataFrame, features: Optional[pd.DataFrame] = None) -> Dict[str, Any]:
    """Count retained-but-excluded rows and their reasons (audit, not deletion)."""
    out: Dict[str, Any] = {"n_rows": int(len(predictions)), "reasons": {}}
    text_cols = [
        c for c in ("note", "notes", "error_message", "rmsd_reason", "reason", "status")
        if c in predictions.columns
    ]
    if "status" in predictions.columns:
        ok = predictions["status"].astype(str).str.strip().str.lower().eq("ok")
    else:
        ok = pd.Series(True, index=predictions.index)
    rmsd_missing = (
        predictions["rmsd"].isna() if "rmsd" in predictions.columns else pd.Series(False, index=predictions.index)
    )
    pb_missing = (
        predictions["posebusters_pass"].isna() if "posebusters_pass" in predictions.columns else pd.Series(False, index=predictions.index)
    )

    def classify(mask) -> Dict[str, int]:
        reasons: Dict[str, int] = {}
        for _, row in predictions.loc[mask].iterrows():
            text = " ".join(str(row.get(c, "")) for c in text_cols).lower()
            if "no_reference_coords" in text or "no reference" in text or "no_ref" in text:
                reason = "no_reference_coords"
            elif "no_pred_sdf" in text:
                reason = "no_prediction_sdf"
            elif "no_results_col" in text:
                reason = "posebusters_pending"
            elif "parse" in text:
                reason = "dock_pdb_parse_fail"
            elif "failed" in text:
                reason = "run_failed"
            elif "nan" in text:
                reason = "metric_nan"
            else:
                reason = "other"
            reasons[reason] = reasons.get(reason, 0) + 1
        return reasons

    out["n_failed_rows"] = int((~ok).sum())
    out["n_ok_without_rmsd"] = int((ok & rmsd_missing).sum())
    out["n_ok_without_posebusters"] = int((ok & pb_missing).sum())
    out["failed_reasons"] = classify(~ok)
    out["no_rmsd_reasons"] = classify(ok & rmsd_missing)
    if features is not None and "rmsd" in features.columns:
        rmsd = pd.to_numeric(features["rmsd"], errors="coerce")
        pb = features["posebusters_pass"] if "posebusters_pass" in features.columns else pd.Series(np.nan, index=features.index)
        out["n_analysis_rows"] = int((rmsd.notna() & pd.Series(pb).notna()).sum())
    return out


def _n_eff_block(manifest: pd.DataFrame, features: pd.DataFrame, cfg=None, pb_definition: str = "strict") -> Dict[str, Any]:
    cfg = cfg or default_config()
    attach = attach_pose_ids(manifest)
    primary = select_analysis_set(attach, "primary_main")
    keys = cluster_keys(primary, "target_cf_subsite")
    counts = keys.value_counts()
    kish = float((counts.sum() ** 2) / np.square(counts.to_numpy(dtype=float)).sum()) if len(counts) else float("nan")
    pair_to_key = dict(zip(primary["pair_id"].astype(str), keys.astype(str)))
    if "pair_id" in features.columns:
        cluster_of_row = features["pair_id"].astype(str).map(pair_to_key)
    else:
        pdb_to_key = dict(zip(primary["pdb_id"].astype(str), keys.astype(str)))
        cluster_of_row = features["pdb_id"].astype(str).map(pdb_to_key)
    cluster_of_row = cluster_of_row.fillna("unknown").to_numpy()
    rmsd_ok = pd.to_numeric(features["rmsd"], errors="coerce").to_numpy(dtype=float) < 2.0
    if pb_definition == "operational" and "pb_pass_rate" in features.columns:
        pb_ok = (
            pd.to_numeric(features["pb_pass_rate"], errors="coerce").to_numpy(dtype=float)
            >= PB_OPERATIONAL_THRESHOLD
        )
    else:
        pb_ok = features["posebusters_pass"].astype("boolean").fillna(False).astype(bool).to_numpy()
    accepted = (rmsd_ok & pb_ok).astype(float)
    icc = design_effect_n_eff(accepted, cluster_of_row) if len(accepted) else {}
    return {
        "n_obs": int(len(features)),
        "n_clusters": int(keys.nunique()),
        "cluster_definition": "target x crystal form x sub-site",
        "pb_definition": pb_definition,
        "kish_n_eff": kish,
        "icc": float(icc.get("icc", float("nan"))),
        "design_effect": float(icc.get("design_effect", float("nan"))),
        "n_eff": float(icc.get("n_eff", float("nan"))),
        "note": (
            "ICC computed on the acceptable indicator (RMSD<2 A & PB) with "
            "target x CF x sub-site clusters; pb_definition='strict' unless "
            "Amendment #6 operational criterion is frozen"
        ),
    }


def five_line_summary(report: Dict[str, Any]) -> List[str]:
    """The 5-line summary requested by the Lead (sent before the full analysis)."""
    audit = report.get("audit", {})
    fc = report.get("failure_classification", {})
    neff = report.get("n_eff", {})
    feat = report.get("features", {})
    cal = report.get("calibration", {})
    bb = report.get("beta_barrel", {})
    if neff.get("n_obs") is None:
        neff_line = "2) n_eff: pending metrics (target x CF x sub-site; expect n_obs=40, n_clusters=15)"
    else:
        neff_line = "2) n_eff: n_obs={} n_clusters={} ICC={} n_eff={} (target x CF x sub-site)".format(
            neff.get("n_obs"), neff.get("n_clusters"),
            None if neff.get("icc") is None else round(float(neff.get("icc", float("nan"))), 3),
            None if neff.get("n_eff") is None else round(float(neff.get("n_eff", float("nan"))), 2),
        )
    if not feat:
        coverage_line = "3) coverage: pending metrics (m_neutral / m_frame_mismatch / coord_frame)"
        calib_line = "4) calibration: pending metrics (AUROC / Brier / ECE / rho)"
    else:
        excl = report.get("exclusions", {})
        coverage_line = (
            "3) coverage: analysis n={}/{} (failed={}, no_rmsd={}: {}); "
            "m_neutral={:.2%} m_frame_mismatch={:.2%} coord_frame={}".format(
                feat.get("n_analysis_rows", "?"), fc.get("n_rows", "?"),
                excl.get("n_failed_rows", "?"), excl.get("n_ok_without_rmsd", "?"),
                excl.get("no_rmsd_reasons", {}),
                float(feat.get("m_neutral_fraction", float("nan"))),
                float(feat.get("m_frame_mismatch_fraction", float("nan"))),
                feat.get("coord_frame_counts", {}),
            )
        )
        calib_line = "4) calibration[n={}, calib_pos={}, eval_pos={}]{}: AUROC={} Brier={} ECE={} rho(MAES,RMSD)={}".format(
            cal.get("n_analysis", report.get("features", {}).get("n_analysis_rows")),
            (report.get("model_metadata") or {}).get("n_calibration_positive"),
            cal.get("n_primary_positives"),
            " [PROVISIONAL: strict PB pending]" if report.get("provisional") else "",
            None if cal.get("auroc_primary") is None else round(float(cal["auroc_primary"]), 3),
            None if cal.get("brier_primary") is None else round(float(cal["brier_primary"]), 3),
            None if cal.get("ece_primary") is None else round(float(cal["ece_primary"]), 3),
            None if cal.get("rho_maes_rmsd") is None else round(float(cal["rho_maes_rmsd"]), 3),
        )
    lines = [
        "1) n_ok/rows: {} ok of {} attempts ({}); failures={}; PB true/false/missing={}; metrics_pending={}".format(
            fc.get("status_counts", {}).get("ok", "?"), fc.get("n_rows", "?"),
            fc.get("status_counts", {}), fc.get("failure_buckets", {}),
            report.get("pb_counts", {}), report.get("metrics_pending", False),
        ),
        neff_line,
        coverage_line,
        calib_line,
        "5) beta-barrel M sensitivity: {} (n={}, min M_cross={})".format(
            bb.get("recommendation", "n/a"), bb.get("n_poses"),
            None if bb.get("min_M_cross") is None else round(float(bb["min_M_cross"]), 3),
        ),
    ]
    return lines


def divergence_cases(
    analysis_frame: pd.DataFrame,
    maes: Sequence[float],
    low_maes: float = 0.2,
    high_maes: float = 0.5,
) -> Dict[str, Any]:
    """MAES-vs-RMSD divergence cases in both directions (SI table).

    ``low_MAES_good_RMSD``: MAES < 0.2 but RMSD < 1 Å (over-penalised good pose);
    ``high_MAES_bad_RMSD``: MAES >= 0.5 but RMSD > 2 Å (over-confident bad pose).
    """
    rmsd = pd.to_numeric(analysis_frame["rmsd"], errors="coerce").to_numpy(dtype=float)
    scores = np.asarray(maes, dtype=float)
    low_good = (scores < float(low_maes)) & (rmsd < 1.0)
    high_bad = (scores >= float(high_maes)) & (rmsd > 2.0)
    # rank discordance (SI requirement cannot always be met by the strict
    # thresholds on a 36-row set): |MAES rank - RMSD rank| >= 10 positions
    maes_rank = pd.Series(scores).rank(ascending=False, method="average").to_numpy(dtype=float)
    rmsd_rank = pd.Series(rmsd).rank(ascending=True, method="average").to_numpy(dtype=float)
    rank_gap = np.abs(maes_rank - rmsd_rank)
    rank_discordant = rank_gap >= 10
    cases = []
    for idx in np.where(low_good | high_bad | rank_discordant)[0]:
        row = analysis_frame.iloc[idx]
        if low_good[idx]:
            direction = "low_MAES_good_RMSD"
        elif high_bad[idx]:
            direction = "high_MAES_bad_RMSD"
        else:
            direction = "rank_discordant"
        cases.append(
            {
                "pose_key": str(row.get("pose_key", row.get("pair_id", ""))),
                "target": str(row.get("system_id", "")),
                "rmsd": round(float(rmsd[idx]), 3),
                "MAES": round(float(scores[idx]), 3),
                "direction": direction,
                "rank_gap": round(float(rank_gap[idx]), 1),
            }
        )
    cases = sorted(cases, key=lambda c: c["MAES"])
    return {
        "n_cases": int(len(cases)),
        "low_MAES_good_RMSD": int(low_good.sum()),
        "high_MAES_bad_RMSD": int(high_bad.sum()),
        "n_rank_discordant": int(rank_discordant.sum()),
        "rank_gap_threshold": 10,
        "thresholds": {"low_maes": float(low_maes), "high_maes": float(high_maes), "good_rmsd": 1.0, "bad_rmsd": 2.0},
        "cases": cases,
    }


def _cluster_labels_for_features(manifest: pd.DataFrame, features: pd.DataFrame) -> np.ndarray:
    """Target × CF × sub-site cluster label per feature row (for wild bootstrap)."""
    attach = attach_pose_ids(manifest)
    primary = select_analysis_set(attach, "primary_main")
    keys = cluster_keys(primary, "target_cf_subsite")
    mapping = dict(zip(primary["pair_id"].astype(str), keys.astype(str)))
    if "pair_id" in features.columns:
        labels = features["pair_id"].astype(str).map(mapping)
    else:  # pragma: no cover - defensive
        labels = pd.Series(["unknown"] * len(features), index=features.index)
    return labels.fillna("unknown").to_numpy()


def run_ingest_report(
    predictions_path: Optional[str | Path] = None,
    manifest_path: Optional[str | Path] = None,
    tag: Optional[str] = None,
    beta_planes_path: Optional[str | Path] = None,
    write: bool = True,
    verbose: bool = True,
    provisional_pb: bool = False,
    pb_definition: str = "strict",
    pb_checks_path: Optional[str | Path] = None,
) -> Dict[str, Any]:
    """Full ingestion + first-real-data validation report (see module docstring)."""
    cfg = default_config()
    preds_path = Path(predictions_path) if predictions_path is not None else DEFAULT_PREDICTIONS
    manifest = load_pose_manifest(Path(manifest_path) if manifest_path is not None else DEFAULT_MANIFEST)
    report: Dict[str, Any] = {
        "predictions_path": str(preds_path),
        "predictions_exists": preds_path.exists(),
        "manifest_path": str(manifest_path or DEFAULT_MANIFEST),
    }
    beta = beta_barrel_m_sensitivity(beta_planes_path or DEFAULT_BETA_PLANES, DEFAULT_STRUCT_DIR, cfg)
    report["beta_barrel"] = beta
    if not preds_path.exists():
        report["status"] = "waiting_for_predictions"
        report["five_line_summary"] = five_line_summary(report)
        if verbose:
            print("\n".join(report["five_line_summary"]))
        return report

    predictions = pd.read_csv(preds_path)
    predictions = normalise_predictions(predictions)
    # auditable independent physics PB recomputed from the per-check matrix
    # (`results/trackB/pb_checks.csv`); the matrix is the authority when it
    # disagrees with B's aggregate bool column.
    checks_path = (
        Path(pb_checks_path) if pb_checks_path is not None
        else (DEFAULT_PB_CHECKS if DEFAULT_PB_CHECKS.exists() else None)
    )
    if checks_path is not None and Path(checks_path).exists():
        checks_table = pd.read_csv(checks_path)
        verification = verify_pb_checks(predictions, checks_table)
        report["pb_independent_recheck"] = {
            key: value for key, value in verification.items() if key != "merged"
        }
        predictions = verification["merged"]
        predictions["posebusters_pass"] = predictions["physics_pass"]
    # operational sensitivity (Amendment #6): overwrite the PB term with
    # pb_pass_rate >= 0.95 before feature/label construction; strict run untouched
    if pb_definition == "operational" and "pb_pass_rate" in predictions.columns:
        predictions["posebusters_pass_strict"] = predictions["posebusters_pass"]
        predictions["posebusters_pass"] = (
            pd.to_numeric(predictions["pb_pass_rate"], errors="coerce") >= PB_OPERATIONAL_THRESHOLD
        )
    report["pb_definition"] = pb_definition
    # provisional marker: the final table can land before strict PoseBusters
    # bools exist (currently 0/36 True because heavy-atom-only checks fail and
    # AddHs fix is pending); V uses continuous pb_pass_rate meanwhile
    pb_true = int((predictions["posebusters_pass"] == True).sum()) if "posebusters_pass" in predictions.columns else 0
    pb_present = int(predictions["posebusters_pass"].notna().sum()) if "posebusters_pass" in predictions.columns else 0
    report["provisional"] = bool(provisional_pb or (pb_present > 0 and pb_true == 0))
    report["pb_label_source"] = (
        "posebusters_pass column (strict bool pending: 0 True so far; V uses pb_pass_rate "
        "continuous values where available)" if report["provisional"] else "posebusters_pass (strict bool)"
    )
    report["pb_counts"] = {
        "n_true": pb_true,
        "n_false": int((predictions["posebusters_pass"] == False).sum()) if "posebusters_pass" in predictions.columns else 0,
        "n_missing": int(predictions["posebusters_pass"].isna().sum()) if "posebusters_pass" in predictions.columns else None,
    }
    audit = check_trackb_inputs(manifest, predictions, min_ok_with_rmsd=MIN_OK_WITH_RMSD)
    report["audit"] = audit
    report["failure_classification"] = failure_classification(predictions)
    report["distributions"] = {
        "rmsd": distribution_summary(predictions["rmsd"]) if "rmsd" in predictions.columns else {},
        "plif_recovery": distribution_summary(predictions["plif_recovery"]) if "plif_recovery" in predictions.columns else {},
        "posebusters_pass_rate": (
            float(pd.Series(predictions["posebusters_pass"]).dropna().astype(bool).mean())
            if "posebusters_pass" in predictions.columns
            and pd.Series(predictions["posebusters_pass"]).notna().any()
            else None
        ),
        "pb_pass_rate_continuous": (
            distribution_summary(predictions["pb_pass_rate"])
            if "pb_pass_rate" in predictions.columns else {}
        ),
        "coord_frame": (
            {str(k): int(v) for k, v in predictions["coord_frame"].astype(str).value_counts().items()}
            if "coord_frame" in predictions.columns else {}
        ),
    }
    if not audit.get("ready"):
        missing_metric_cols = {"rmsd", "posebusters_pass"}
        fc = report.get("failure_classification", {})
        incomplete = set(audit.get("missing", [])) | set(audit.get("empty_required", []))
        report["metrics_pending"] = bool(
            missing_metric_cols.intersection(incomplete)
            and fc.get("n_ok_with_missing_metrics")
        )
        report["status"] = "awaiting_metrics" if report["metrics_pending"] else "not_ready"
        report["five_line_summary"] = five_line_summary(report)
        if verbose:
            print("\n".join(report["five_line_summary"]))
        return report

    # Features are computed only on rows that carry labels: failed runs and
    # rmsd=NA(no_reference_coords) rows have no valid V/label inputs, so they are
    # excluded BEFORE feature construction (their count/reasons are recorded in
    # report["exclusions"] from the full table).
    analysis_ready = (
        predictions["rmsd"].notna()
        & predictions["posebusters_pass"].notna()
    ) if "rmsd" in predictions.columns and "posebusters_pass" in predictions.columns else pd.Series(False, index=predictions.index)
    analysis_input = predictions.loc[analysis_ready].copy()
    report["exclusions"] = exclusion_breakdown(predictions)
    report["exclusions"]["n_analysis_rows"] = int(len(analysis_input))
    if len(analysis_input) < MIN_ANALYSIS_ROWS:
        report["status"] = "insufficient_labels"
        report["features"] = {"n_rows": int(len(analysis_input)), "n_analysis_rows": int(len(analysis_input))}
        report["five_line_summary"] = five_line_summary(report)
        if verbose:
            print("\n".join(report["five_line_summary"]))
        return report
    features = build_features_from_manifest(manifest, analysis_input, subset="primary_main")
    n_analysis = int(len(features))
    report["features"] = {
        "n_rows": n_analysis,
        "n_analysis_rows": n_analysis,
        "n_rows_excluded_from_calibration": int(len(predictions) - n_analysis),
        "n_pairs": int(features["pair_id"].nunique()),
        "m_neutral_fraction": float(features["m_neutral"].astype(bool).mean()),
        "m_frame_mismatch_fraction": float(features["m_frame_mismatch"].astype(bool).mean()),
        "m_confidence_counts": {str(k): int(v) for k, v in features["m_confidence"].value_counts().items()},
        "coord_frame_counts": (
            {str(k): int(v) for k, v in features["coord_frame"].astype(str).value_counts().items()}
            if "coord_frame" in features.columns else {}
        ),
    }
    analysis_frame = features
    try:
        model = MAESModel(cfg).fit(analysis_frame)
        report["calibration"] = evaluate_predictions(analysis_frame, model)
        report["model_metadata"] = model.metadata
        cal = report["calibration"]
        if isinstance(cal, dict) and "error" not in cal:
            # stable fill keys for the Results draft (notes/maes_r2_fillkeys.md)
            cal["brier"] = cal.get("brier_primary")
            cal["ece"] = cal.get("ece_primary")
            cal["rho_rmsd"] = cal.get("rho_maes_rmsd")
            cal["n_analysis"] = cal.get("n")
            cal["n_primary_positives"] = cal.get("n_eval_positive_primary")
            cal["secondary_auroc_plif"] = cal.get("auroc_secondary")
        report["divergence"] = divergence_cases(analysis_frame, model.predict(analysis_frame))
        ablation = ablation_study(analysis_frame, cfg).to_dict(orient="records")
        report["ablation"] = ablation
        report["ablation_by_variant"] = {
            str(row.get("variant")): {k: v for k, v in row.items() if k != "variant"}
            for row in ablation
        }
        # ------------------------------------------------ robustness (audit M3)
        cluster_labels = _cluster_labels_for_features(manifest, analysis_frame)
        y_eval = build_evaluation_labels(analysis_frame, cfg.labels, which="primary")
        wild = wild_cluster_bootstrap_auc(
            y_eval,
            model.predict(analysis_frame),
            cluster_labels,
            n_boot=2000,
            seed=cfg.scores.random_state,
        )
        report["robustness"] = {
            "auc_eval_primary": wild,
            "n_clusters": wild.get("n_clusters"),
            "note": (
                "15 clusters only: wild cluster bootstrap reported alongside the "
                "ligand-level bootstrap; CIs may be narrow (few clusters) — limitation"
            ),
        }
        cal["auc_wild_ci_lo"] = wild.get("ci_lo")
        cal["auc_wild_ci_hi"] = wild.get("ci_hi")
        # ------------------------------------------ weight-leakage audit (M4)
        entropy_w = entropy_weights(analysis_frame)
        # config.validate requires strictly positive weights; clip exact zeros
        # (constant sub-score columns can get zero entropy weight) and renormalise
        entropy_w = {key: max(float(value), 1e-6) for key, value in entropy_w.items()}
        entropy_total = sum(entropy_w.values())
        entropy_w = {key: value / entropy_total for key, value in entropy_w.items()}
        cfg_ent = copy.deepcopy(cfg)
        cfg_ent.scores.weights = dict(entropy_w)
        cfg_ent.validate()
        model_ent = MAESModel(config=cfg_ent).fit(analysis_frame)
        entropy_metrics = evaluate_predictions(analysis_frame, model_ent)
        report["weight_sensitivity"] = {
            "fixed": {
                "weights": dict(model.weights),
                "auroc": cal.get("auroc_primary"),
                "brier": cal.get("brier_primary"),
                "ece": cal.get("ece_primary"),
            },
            "entropy": {
                "weights": dict(entropy_w),
                "auroc": entropy_metrics.get("auroc_primary"),
                "brier": entropy_metrics.get("brier_primary"),
                "ece": entropy_metrics.get("ece_primary"),
            },
            "note": (
                "entropy weights estimated on the calibration rows only (in-sample "
                "sensitivity); the primary fit uses fixed a-priori V=M=S=1/3 weights "
                "frozen at inference — no data-derived weighting in fit (audit M4)"
            ),
        }
        # --------------------------------------- small-n calibration policy (M3)
        report["calibration_policy"] = {
            "primary": "platt",
            "isotonic_min_positives": (model.metadata.get("calibrator_policy") or {}).get(
                "isotonic_min_positives"
            ),
            "reliability_bins": 4,
            "irt": "descriptive only (36 observations x 15 clusters — no Rasch/2PL inference)",
        }
    except Exception as exc:  # single-class labels etc.
        report["calibration"] = {"error": f"{type(exc).__name__}: {exc}"}
        report["ablation"] = []
        report["ablation_by_variant"] = {}
        report["divergence"] = {"n_cases": None, "note": "not computed (calibration failed)"}
    report["n_eff"] = _n_eff_block(manifest, analysis_frame, cfg, pb_definition=pb_definition)
    report["irt"] = {
        "n_clusters": report["n_eff"].get("n_clusters"),
        "icc": report["n_eff"].get("icc"),
        "n_eff": report["n_eff"].get("n_eff"),
        "pb_definition": report["n_eff"].get("pb_definition"),
        "cluster_definition": report["n_eff"].get("cluster_definition"),
    }
    report["status"] = "ok"

    report["five_line_summary"] = five_line_summary(report)
    if verbose:
        print("\n".join(report["five_line_summary"]))
    if write:
        FIGURES_DIR.mkdir(parents=True, exist_ok=True)
        suffix = tag or pd.Timestamp.now().strftime("%Y%m%d_%H%M")
        out_path = FIGURES_DIR / f"ingest_report_{suffix}.json"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
        report["report_path"] = str(out_path)
    return report
