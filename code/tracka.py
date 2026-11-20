"""Track A enrichment analysis framework (W3).

Reads one or more ``<target>_scores.csv`` tables (B-line output; one row per
ligand, label 1 = active, 0 = decoy) and produces, per (target, method,
metal-condition):

* EF1% (reported only when n_active >= 150, per reporting rule R1),
* EF5%, BEDROC (alpha = 20), AUROC,
* percentile bootstrap CIs over ligands (2.5 / 97.5), plus
* BH-adjusted p-values for AUROC vs 0.5 (Hanley-McNeil SE), and
* metal / non-metal stratification (LpxC = Zn; GyrA/GyrB = Mg; DHFR/FabI/LepB/
  TEM-1 = none; ``reactive`` metadata from ``evaluation.TARGET_METAL_DEPENDENCE``).

The full Track A run (B-line ``results/trackA/full/<target>_scores.csv``) uses
the same entry point::

    python -m maes.tracka --scores-dir ../results/trackA/full --tag full
    python -m maes.tracka --scores-dir ../results/trackA/pilot_v2 --tag pilot_v2

Outputs: ``figures/trackA_enrichment_<tag>.json`` and ``.csv``.
"""

from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import norm

from .evaluation import (
    BEDROC_ALPHA,
    TARGET_METAL_DEPENDENCE,
    auroc,
    bedroc,
    ef1_allowed,
    enrichment_factor,
    enrichment_metrics,
    metal_stratum,
)
from .properties import (
    compute_rdkit_properties,
    largest_ecfp_cluster,
    largest_scaffold_group,
    load_ligand_properties,
    property_sensitivity,
)

__all__ = [
    "IncompleteTrackADataError",
    "bootstrap_enrichment_ci",
    "compute_tracka_metrics",
    "load_score_table",
    "run_tracka_framework",
]


class IncompleteTrackADataError(RuntimeError):
    """Score tables are shorter than the target_metadata ``n_rows`` (run in progress)."""

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TRACKA_DIR = PROJECT_ROOT / "results" / "trackA"
FIGURES_DIR = PROJECT_ROOT / "figures"

#: per-target docking metadata (Amendment #7): decoy subsample + exhaustiveness
METADATA_FIELDS = (
    "decoy_ratio",
    "sampling_seed",
    "decoy_subsample_sha256",
    "exhaustiveness",
    "panel_role",
    "box_size_A",
    "n_prep_missing",
    "n_rows",
    "n_act",
    "n_dec",
)
DEFAULT_METADATA_CANDIDATES = (
    DEFAULT_TRACKA_DIR / "full" / "target_metadata.csv",
    DEFAULT_TRACKA_DIR / "target_metadata.csv",
    PROJECT_ROOT / "data" / "trackA_target_metadata.csv",
)

#: per-ligand box-clipping audit (Amendment #10)
DEFAULT_BOX_CLIPPED_CANDIDATES = (
    DEFAULT_TRACKA_DIR / "box_clipped_ligands.csv",
    PROJECT_ROOT / "data" / "trackA" / "box_clipped_ligands.csv",
)
BOX_CLIPPED_FIELDS = ("box_clipped", "max_dim_A", "box_size_A")

DEFAULT_NOVEL_SPLITS_DIR = PROJECT_ROOT / "data" / "trackA" / "splits"
DEFAULT_DECOYS_DIR = PROJECT_ROOT / "data" / "decoys"
DEFAULT_PANEL_READY = DEFAULT_TRACKA_DIR / "PANEL_READY.json"


def _load_chemotype_evidence() -> Optional[List[Dict[str, Any]]]:
    path = PROJECT_ROOT / "data" / "trackA" / "chemotype_analysis.csv"
    if not path.exists():
        return None
    frame = pd.read_csv(path)
    keep = [
        "target", "n_actives", "n_scaffolds_unique", "max_scaffold_frac",
        "top5_scaffold_frac", "n_ecfp_clusters", "max_ecfp_cluster_frac",
        "top5_ecfp_cluster_frac", "internal_maxT_median", "pchembl_median",
        "n_documents",
    ]
    columns = [column for column in keep if column in frame.columns]
    return frame[columns].to_dict(orient="records")


def _load_decoy_overlap_evidence() -> Optional[List[Dict[str, Any]]]:
    path = PROJECT_ROOT / "data" / "decoys" / "decoy_overlap_lpxc.csv"
    if not path.exists():
        return None
    return pd.read_csv(path).to_dict(orient="records")


def load_panel_ready_signal(path: Optional[str | Path] = None) -> Optional[Dict[str, Any]]:
    """B-line completion signal: ``results/trackA/full/PANEL_READY.json``.

    The formal 9/9 run may be triggered by B's message/signal; this only records
    the signal in the report (panel completeness itself is still checked from
    the score tables / metadata).
    """
    target_path = Path(path) if path is not None else DEFAULT_PANEL_READY
    if not target_path.exists():
        return None
    try:
        payload = json.loads(target_path.read_text(encoding="utf-8"))
    except Exception:
        payload = {}
    if not isinstance(payload, dict):
        payload = {"payload": payload}
    payload.setdefault("signal_path", str(target_path))
    return payload

#: expected Track A main panel (9 targets); the formal JSON is written only when
#: every expected target has a completed score table (Lead 2026-11-17)
EXPECTED_TRACKA_PANEL = (
    "LpxC_Ec",
    "LpxC_Pa",
    "DHFR_Ec",
    "Gyrase_Ec",
    "GyrA_Ec",
    "GyrB_Ec",
    "FabI_Ec",
    "LepB_Ec",
    "TEM1_Ec",
)
_DECOY_CLUSTER_CACHE: Dict[str, Optional[Dict[str, Any]]] = {}


def _load_decoy_cluster_map(target: str) -> Optional[Dict[str, Any]]:
    """Map decoy ChEMBL ids to their matched active (inference cluster source).

    Reads ``data/decoys/<target>_decoys_s1.csv`` (A-line DUD-E-style table) and
    returns ``{decoy_norm_id -> {'inchikey':..., 'active_chembl':...}}`` plus
    ``active_chembl -> inchikey``.  Falls back to None when the pool is absent.
    """
    target = str(target)
    if target in _DECOY_CLUSTER_CACHE:
        return _DECOY_CLUSTER_CACHE[target]
    path = DEFAULT_DECOYS_DIR / f"{target}_decoys_s1.csv"
    if not path.exists():
        _DECOY_CLUSTER_CACHE[target] = None
        return None
    pool = pd.read_csv(
        path,
        usecols=lambda column: column in {
            "decoy_chembl_id", "decoy_inchikey",
            "matched_active_inchikey", "matched_active_compound_chembl_id",
        },
    )
    by_decoy: Dict[str, Dict[str, Any]] = {}
    by_decoy_inchikey: Dict[str, Dict[str, Any]] = {}
    active_to_inchikey: Dict[str, str] = {}
    for _, row in pool.iterrows():
        decoy_key = _normalize_ligand_key(row.get("decoy_chembl_id", ""))
        active_chembl = row.get("matched_active_compound_chembl_id")
        active_chembl = "" if pd.isna(active_chembl) else _normalize_ligand_key(active_chembl)
        inchikey = row.get("matched_active_inchikey")
        inchikey = "" if pd.isna(inchikey) else str(inchikey).strip()
        decoy_inchikey = row.get("decoy_inchikey")
        decoy_inchikey = "" if pd.isna(decoy_inchikey) else str(decoy_inchikey).strip()
        info = {"inchikey": inchikey, "active_chembl": active_chembl, "decoy_inchikey": decoy_inchikey}
        if decoy_key:
            by_decoy[decoy_key] = info
        if decoy_inchikey:
            by_decoy_inchikey[decoy_inchikey] = info
        if active_chembl and inchikey:
            active_to_inchikey[active_chembl] = inchikey
    result = {
        "by_decoy": by_decoy,
        "by_decoy_inchikey": by_decoy_inchikey,
        "active_to_inchikey": active_to_inchikey,
        "source": str(path),
    }
    _DECOY_CLUSTER_CACHE[target] = result
    return result


def cluster_bootstrap_enrichment_ci(
    labels: Sequence[int],
    scores: Sequence[float],
    clusters: Sequence[Any],
    n_boot: int = 2000,
    seed: int = 42,
    alpha: float = BEDROC_ALPHA,
) -> Dict[str, Any]:
    """Cluster bootstrap CI: resample **clusters** (matched active groups), not rows."""
    y = np.asarray(labels, dtype=int)
    s = np.asarray(scores, dtype=float)
    c = np.asarray([str(value) for value in clusters])
    valid = np.isfinite(s)
    y, s, c = y[valid], s[valid], c[valid]
    if len(y) == 0 or y.sum() == 0 or y.sum() == len(y):
        return {"n_clusters": 0, "n_boot": 0}
    unique, inverse = np.unique(c, return_inverse=True)
    groups = [np.where(inverse == i)[0] for i in range(len(unique))]
    rng = np.random.default_rng(int(seed))
    draws: Dict[str, List[float]] = {"auc": [], "ef5_percent": [], "bedroc20": []}
    if ef1_allowed(int(y.sum())):
        draws["ef1_percent"] = []
    for _ in range(int(n_boot)):
        chosen = rng.integers(0, len(unique), size=len(unique))
        idx = np.concatenate([groups[j] for j in chosen])
        yb, sb = y[idx], s[idx]
        if yb.sum() == 0 or yb.sum() == len(yb):
            continue
        for name, value in (
            ("auc", auroc(yb, sb)),
            ("ef5_percent", enrichment_factor(yb, sb, 0.05)),
            ("bedroc20", bedroc(yb, sb, alpha=alpha)),
        ):
            if np.isfinite(value):
                draws[name].append(float(value))
        if "ef1_percent" in draws:
            value = enrichment_factor(yb, sb, 0.01)
            if np.isfinite(value):
                draws["ef1_percent"].append(float(value))
    ci = {
        name: (float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5)))
        for name, values in draws.items()
        if len(values) >= 20
    }
    return {
        "ci": ci,
        "n_clusters": int(len(unique)),
        "n_boot": int(n_boot),
        "n_nonempty_resamples": int(len(draws.get("auc", []))),
    }


def _load_novel_scaffold_splits(splits_dir: Optional[str | Path] = None) -> pd.DataFrame:
    """Amendment #9 novel-scaffold test subset files (``*_temporal_novel.csv``)."""
    base = Path(splits_dir) if splits_dir is not None else DEFAULT_NOVEL_SPLITS_DIR
    files = sorted(base.glob("*_temporal_novel.csv"))
    if not files:
        return pd.DataFrame(columns=["target", "ligand_id_norm", "ligand_id", "max_tanimoto", "murcko_scaffold"])
    parts = []
    for path in files:
        frame = pd.read_csv(path)
        target = str(frame["target_key"].iloc[0]) if "target_key" in frame.columns and len(frame) else path.stem.split("_temporal")[0]
        ligand_col = "compound_chembl_id" if "compound_chembl_id" in frame.columns else (
            "ligand_id" if "ligand_id" in frame.columns else frame.columns[0]
        )
        rows = pd.DataFrame(
            {
                "target": target,
                "ligand_id": frame[ligand_col].astype(str),
                "ligand_id_norm": frame[ligand_col].map(_normalize_ligand_key),
                "max_tanimoto": pd.to_numeric(frame.get("max_tanimoto_to_train"), errors="coerce"),
                "murcko_scaffold": frame.get("murcko_scaffold"),
                "source_file": path.name,
            }
        )
        parts.append(rows)
    return pd.concat(parts, ignore_index=True).drop_duplicates(subset=["target", "ligand_id_norm"])


def _to_bool_value(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = str(value).strip().lower()
    if "clip" in text:
        return True
    if text in {"false", "0", "0.0", "no", "ok", "none", "nan", ""}:
        return False
    return text in {"true", "1", "1.0", "yes", "y"}


def _normalize_ligand_key(value: Any) -> str:
    """Normalise ligand ids so ``act_CHEMBL123`` matches ``CHEMBL123``."""
    text = str(value).strip()
    match = re.search(r"CHEMBL\d+", text, flags=re.IGNORECASE)
    if match:
        return match.group(0).upper()
    return text


def _normalise_status(value: Any) -> str:
    """Map B-line statuses to the analysis status; ``ok(reuse)`` counts as ok."""
    text = str(value).strip().lower()
    if text in {"ok", "ok(reuse)", "ok (reuse)", "ok_reuse"}:
        return "ok"
    return text


def load_box_clipped_ligands(path: Optional[str | Path] = None) -> pd.DataFrame:
    """Load A-line per-ligand box-clipping flags (Amendment #10).

    Accepts the A-line schema ``target, compound_id, inchikey, max_dim_A,
    box_size_A, flag`` (flag = ``box_clipped``) and tolerant variants:
    ``compound_id`` is mapped to ``ligand_id``; ``flag`` is mapped to
    ``box_clipped``; when only ``max_dim_A``/``box_size_A`` are present the flag
    is derived as ``max_dim_A > box_size_A``.  Ligand ids are normalised
    (``act_CHEMBL123`` -> ``CHEMBL123``) for joining against score tables.
    Missing file -> empty frame (no exclusion).
    """
    candidates: List[Path] = []
    if path is not None:
        candidates = [Path(path)]
    else:
        candidates = [p for p in DEFAULT_BOX_CLIPPED_CANDIDATES if p.exists()]
    if not candidates or not candidates[0].exists():
        return pd.DataFrame(columns=["target", "ligand_id", "ligand_id_norm", *BOX_CLIPPED_FIELDS])
    frame = pd.read_csv(candidates[0])
    if "target" not in frame.columns:
        raise KeyError("box_clipped_ligands.csv needs a target column")
    if "ligand_id" not in frame.columns:
        for alias in ("compound_id", "compound_chembl_id", "ligand", "inchikey"):
            if alias in frame.columns:
                frame = frame.rename(columns={alias: "ligand_id"})
                break
    if "ligand_id" not in frame.columns:
        raise KeyError("box_clipped_ligands.csv needs ligand_id / compound_id")
    if "box_clipped" not in frame.columns:
        if "flag" in frame.columns:
            frame["box_clipped"] = frame["flag"].map(_to_bool_value)
        elif "max_dim_A" in frame.columns and "box_size_A" in frame.columns:
            frame["box_clipped"] = (
                pd.to_numeric(frame["max_dim_A"], errors="coerce")
                > pd.to_numeric(frame["box_size_A"], errors="coerce")
            )
        else:
            raise KeyError(
                "box_clipped_ligands.csv needs 'flag'/'box_clipped' or "
                "max_dim_A + box_size_A columns"
            )
    frame["target"] = frame["target"].astype(str)
    frame["ligand_id"] = frame["ligand_id"].astype(str)
    frame["ligand_id_norm"] = frame["ligand_id"].map(_normalize_ligand_key)
    frame["box_clipped"] = frame["box_clipped"].map(_to_bool_value)
    keep = ["target", "ligand_id", "ligand_id_norm"] + [c for c in BOX_CLIPPED_FIELDS if c in frame.columns]
    return frame[keep].drop_duplicates(subset=["target", "ligand_id_norm"], keep="first")


def load_target_metadata(path: Optional[str | Path] = None) -> Dict[str, Dict[str, Any]]:
    """Load per-target docking metadata (decoy ratio / seed / SHA256 / exh / role).

    Accepts a CSV/JSON with a ``target``/``target_id`` column; when no path is
    given, looks for the standard Track A metadata files.  Missing files return
    an empty mapping (the framework then falls back to score-file columns or the
    observed decoy ratio).
    """
    candidates: List[Path] = []
    if path is not None:
        candidates = [Path(path)]
    else:
        candidates = [p for p in DEFAULT_METADATA_CANDIDATES if p.exists()]
    if not candidates or not candidates[0].exists():
        return {}
    source = candidates[0]
    if source.suffix.lower() == ".json":
        payload = json.loads(source.read_text(encoding="utf-8"))
        if isinstance(payload, dict) and "targets" in payload:
            payload = payload["targets"]
        records = payload if isinstance(payload, list) else [
            {"target": key, **value} for key, value in payload.items() if isinstance(value, dict)
        ]
        frame = pd.DataFrame(records)
    else:
        frame = pd.read_csv(source)
    key_col = "target" if "target" in frame.columns else ("target_id" if "target_id" in frame.columns else None)
    if key_col is None:
        return {}
    out: Dict[str, Dict[str, Any]] = {}
    for _, row in frame.iterrows():
        target = str(row[key_col])
        values: Dict[str, Any] = {}
        for field in METADATA_FIELDS:
            if field not in frame.columns:
                continue
            value = row[field]
            if pd.isna(value):
                continue
            values[field] = value.item() if hasattr(value, "item") else value
        # tolerate the short alias `box_size` from older metadata exports
        if "box_size_A" not in values and "box_size" in frame.columns:
            alias_value = row["box_size"]
            if not pd.isna(alias_value):
                values["box_size_A"] = (
                    alias_value.item() if hasattr(alias_value, "item") else alias_value
                )
        out[target] = values
    return out

LABEL_COLUMNS = ("label", "active", "is_active", "y", "outcome")
SCORE_COLUMN_RE = re.compile(r"^(?P<method>.+?)_score$", flags=re.IGNORECASE)
ID_LIKE = {
    "ligand_id", "ligand", "compound_id", "compound_chembl_id", "inchikey", "smiles",
    "target", "method", "status", "n_modes", "returncode", "elapsed_s", "seed",
    "exhaustiveness", "label", "active", "is_active", "y", "outcome",
}
#: scoring functions whose raw values are better when LOWER (negated for ranking)
LOWER_BETTER_HINTS = ("vina", "drugclip", "unidock", "glide", "score")
HIGHER_BETTER_HINTS = ("boltz", "af3", "confidence", "prob", "rank_confidence")

#: filename stem -> Track A target id (used by the pilot files and C-line naming)
TARGET_ALIASES = {    "lpxc": "LpxC_Ec",
    "lpxc_ec": "LpxC_Ec",
    "lpxc_pa": "LpxC_Pa",
    "lpxc_zn": "LpxC_Ec",
    "lpxc_nozn": "LpxC_Ec",
    "pilot_scores": "LpxC_Ec",
    "pilot_lpxc": "LpxC_Ec",
    "docking_scores": "LpxC_Ec",
    "dhfr": "DHFR_Ec",
    "dhfr_ec": "DHFR_Ec",
    "gyra": "GyrA_Ec",
    "gyra_ec": "GyrA_Ec",
    "gyrase": "Gyrase_Ec",
    "gyrase_ec": "Gyrase_Ec",
    "gyrase_ec_scores": "Gyrase_Ec",
    "gyrb": "GyrB_Ec",
    "gyrb_ec": "GyrB_Ec",
    "fabi": "FabI_Ec",
    "fabi_ec": "FabI_Ec",
    "lepb": "LepB_Ec",
    "lepb_ec": "LepB_Ec",
    "tem1": "TEM1_Ec",
    "tem1_ec": "TEM1_Ec",
    "tem-1": "TEM1_Ec",
}


# --------------------------------------------------------------------- loading
def _label_column(df: pd.DataFrame) -> Optional[str]:
    for column in LABEL_COLUMNS:
        if column in df.columns:
            return column
    return None


def _score_columns(df: pd.DataFrame) -> Dict[str, str]:
    """Map method name -> score column (``score``/``vina_score``/``<method>_score``)."""
    found: Dict[str, str] = {}
    for column in df.columns:
        lowered = str(column).lower()
        if lowered in ID_LIKE or lowered in LABEL_COLUMNS:
            continue
        if lowered == "score":
            found.setdefault("vina", column)
            continue
        match = SCORE_COLUMN_RE.match(str(column))
        if match:
            method = match.group("method").strip().lower()
            if method in {"vina", "vina_zn", "vina_nozn"}:
                found.setdefault("vina", column) if method == "vina" else found.setdefault(method, column)
            else:
                found.setdefault(method, column)
    if not found:  # fallback: single numeric non-id column
        numeric = [
            c for c in df.columns
            if str(c).lower() not in ID_LIKE and str(c).lower() not in LABEL_COLUMNS
            and pd.api.types.is_numeric_dtype(df[c])
        ]
        if numeric:
            found["vina"] = numeric[0]
    return found


def _resolve_target(stem: str, target_map: Optional[Dict[str, str]]) -> str:
    key = stem.strip().lower()
    if target_map:
        for alias, target in target_map.items():
            if str(alias).lower() == key:
                return str(target)
    if key in TARGET_ALIASES:
        return TARGET_ALIASES[key]
    key = key.removesuffix("_scores").removesuffix("_score")
    if key in TARGET_ALIASES:
        return TARGET_ALIASES[key]
    # `<target>_<condition>` fallback, e.g. `gyrb_pa` -> `GyrB_Pa`
    parts = key.split("_")
    if len(parts) >= 2:
        base = "_".join(parts[:1] + parts[-1:])
        if base in TARGET_ALIASES:
            return TARGET_ALIASES[base]
    return stem


def _condition_from_stem(stem: str) -> str:
    lowered = stem.lower()
    if "nozn" in lowered or "no_zn" in lowered:
        return "no_Zn"
    if "_zn" in lowered or lowered.endswith("zn"):
        return "with_Zn"
    return "native"


def _score_direction(method: str, column: str) -> str:
    """Ranking direction for a raw score column (higher/lower is better)."""
    text = f"{method} {column}".lower()
    if any(hint in text for hint in HIGHER_BETTER_HINTS):
        return "higher_better"
    if any(hint in text for hint in LOWER_BETTER_HINTS):
        return "lower_better"
    return "higher_better"


def load_score_table(
    path: str | Path,
    target_map: Optional[Dict[str, str]] = None,
    condition: Optional[str] = None,
    score_direction: str = "auto",
) -> pd.DataFrame:
    """Load a ``<target>_scores.csv`` table into long format.

    Columns: ``target, method, condition, ligand_id, label, score, raw_score,
    score_direction``.  ``score`` is oriented for ranking (higher = better):
    raw Vina/DrugCLIP scores are negated by default because lower is better;
    pass ``score_direction='higher_better'`` to disable the auto-convention.
    Rows with a failed/NA score are dropped.
    """
    path = Path(path)
    frame = pd.read_csv(path)
    label_col = _label_column(frame) or "label"
    if label_col not in frame.columns:
        raise KeyError(f"{path.name}: no label column (expected one of {LABEL_COLUMNS})")
    stem = path.stem
    target = _resolve_target(stem, target_map)
    condition = condition or _condition_from_stem(stem)
    score_cols = _score_columns(frame)
    if not score_cols:
        raise KeyError(f"{path.name}: no score column found")
    ligand_col = "ligand_id" if "ligand_id" in frame.columns else (
        "compound_id" if "compound_id" in frame.columns else (
            "compound_chembl_id" if "compound_chembl_id" in frame.columns else (
                "ligand" if "ligand" in frame.columns else frame.columns[0]
            )
        )
    )
    label = pd.to_numeric(frame[label_col], errors="coerce").fillna(0).astype(int)
    status_raw = frame["status"].astype(str).str.strip().str.lower() if "status" in frame.columns else None
    status = status_raw.map(_normalise_status) if status_raw is not None else None
    reuse_raw = (
        status_raw.str.contains("reuse", na=False) if status_raw is not None
        else pd.Series(False, index=frame.index)
    )
    file_meta: Dict[str, Any] = {}
    for field in METADATA_FIELDS:  # metadata may be embedded as constant columns
        column = field if field in frame.columns else ("box_size" if field == "box_size_A" and "box_size" in frame.columns else None)
        if column is not None:
            values = frame[column].dropna()
            if len(values):
                value = values.iloc[0]
                file_meta[field] = value.item() if hasattr(value, "item") else value
    # A-line score tables now embed split/box context columns; carry them through
    per_row_context = [c for c in ("temporal_split", "document_year", "box_clipped") if c in frame.columns]
    box_mask = (
        frame["box_clipped"].map(_to_bool_value)
        if "box_clipped" in frame.columns else pd.Series(False, index=frame.index)
    )
    rows: List[pd.DataFrame] = []
    for method, column in score_cols.items():
        raw = pd.to_numeric(frame[column], errors="coerce")
        mask = raw.notna()
        if status is not None:
            mask &= status == "ok"
        # box-clipped ligands may have no docked score; keep them for the
        # Amendment #10 audit (they are removed before enrichment metrics)
        if "box_clipped" in frame.columns:
            mask = mask | box_mask
        direction = (
            score_direction if score_direction in {"higher_better", "lower_better"}
            else _score_direction(str(method), str(column))
        )
        ranked = raw if direction == "higher_better" else -raw
        # prep_missing / other non-ok rows are reported for the exclusion audit
        status_series = status if status is not None else pd.Series("ok", index=frame.index)
        prep_mask = status_series == "prep_missing"
        other_mask = ~status_series.isin({"ok", "prep_missing"})
        other_active_mask = other_mask & (label == 1) & ~box_mask
        other_decoy_mask = other_mask & (label == 0) & ~box_mask
        reuse_mask = reuse_raw & mask  # accepted rows that arrived as ok(reuse)
        payload = {
            "target": target,
            "method": method,
            "condition": condition,
            "ligand_id": frame.loc[mask, ligand_col].astype(str),
            "ligand_id_norm": frame.loc[mask, ligand_col].map(_normalize_ligand_key),
            "label": label[mask].to_numpy(dtype=int),
            "score": ranked[mask].to_numpy(dtype=float),
            "raw_score": raw[mask].to_numpy(dtype=float),
            "score_direction": direction,
            "source_file": path.name,
            "n_ok_reuse": int(reuse_mask.sum()),
            "n_prep_missing": int(prep_mask.sum()),
            "n_prep_missing_actives": int((prep_mask & (label == 1)).sum()),
            "n_prep_missing_decoys": int((prep_mask & (label == 0)).sum()),
            "n_non_ok_other": int(other_mask.sum()),
            "n_non_ok_other_actives": int(other_active_mask.sum()),
            "n_non_ok_other_decoys": int(other_decoy_mask.sum()),
            "n_actives_table_total": int(label.sum()),
        }
        for extra in per_row_context:
            payload[extra] = frame.loc[mask, extra].to_numpy()
        payload.update(
            {field: file_meta[field] for field in METADATA_FIELDS if field in file_meta and field not in payload}
        )
        rows.append(pd.DataFrame(payload))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=[
            "target", "method", "condition", "ligand_id", "label",
            "score", "raw_score", "score_direction", "source_file",
        ]
    )


# --------------------------------------------------------------------- metrics
def compute_tracka_metrics(labels: Sequence[int], scores: Sequence[float]) -> Dict[str, Any]:
    """R1-compliant point estimates for one target/method entry."""
    y = np.asarray(labels, dtype=int)
    s = np.asarray(scores, dtype=float)
    valid = np.isfinite(s)
    y, s = y[valid], s[valid]
    n_active = int(y.sum())
    metrics = enrichment_metrics(y, s, n_active=n_active, alpha=BEDROC_ALPHA)
    ef5 = enrichment_factor(y, s, 0.05)
    top10 = (
        float(y[np.argsort(-s, kind="mergesort")[: min(10, len(y))]].sum()) if len(y) else float("nan")
    )
    return {
        "n": int(len(y)),
        "n_active": n_active,
        "n_decoy": int(len(y) - n_active),
        "ef1_percent": metrics["ef1"],  # None when n_active < 150 (R1)
        "ef1_allowed": bool(metrics["ef1_allowed"]),
        "ef1_suppressed": bool(metrics["ef1_suppressed"]),
        "primary_metric": metrics["primary_metric"],
        "ef5_percent": float(ef5),
        "bedroc20": float(metrics["bedroc"]),
        "auc": float(metrics["auroc"]),
        "top10_actives": top10,
    }


def bootstrap_enrichment_ci(
    labels: Sequence[int],
    scores: Sequence[float],
    n_boot: int = 2000,
    seed: int = 42,
    alpha: float = BEDROC_ALPHA,
) -> Dict[str, Tuple[float, float]]:
    """Percentile bootstrap CI over ligands (AUC / EF5 / BEDROC / EF1 when allowed)."""
    y = np.asarray(labels, dtype=int)
    s = np.asarray(scores, dtype=float)
    valid = np.isfinite(s)
    y, s = y[valid], s[valid]
    if len(y) == 0 or y.sum() == 0 or y.sum() == len(y):
        return {}
    rng = np.random.default_rng(int(seed))
    draws: Dict[str, List[float]] = {"auc": [], "ef5_percent": [], "bedroc20": []}
    if ef1_allowed(int(y.sum())):
        draws["ef1_percent"] = []
    n = len(y)
    for _ in range(int(n_boot)):
        idx = rng.integers(0, n, size=n)
        yb, sb = y[idx], s[idx]
        if yb.sum() == 0 or yb.sum() == len(yb):
            continue
        a = auroc(yb, sb)
        ef5 = enrichment_factor(yb, sb, 0.05)
        b = bedroc(yb, sb, alpha=alpha)
        for name, value in (("auc", a), ("ef5_percent", ef5), ("bedroc20", b)):
            if np.isfinite(value):
                draws[name].append(float(value))
        if "ef1_percent" in draws:
            ef1 = enrichment_factor(yb, sb, 0.01)
            if np.isfinite(ef1):
                draws["ef1_percent"].append(float(ef1))
    return {
        name: (float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5)))
        for name, values in draws.items()
        if len(values) >= 20
    }


def _hanley_mcneil_p(auc_value: float, n_active: int, n_decoy: int) -> float:
    """Two-sided p for AUROC vs 0.5 using the Hanley-McNeil standard error."""
    if not np.isfinite(auc_value) or n_active < 2 or n_decoy < 2:
        return float("nan")
    if auc_value <= 0.0 or auc_value >= 1.0:
        return 0.0
    q1 = auc_value / (2.0 - auc_value)
    q2 = 2.0 * auc_value * auc_value / (1.0 + auc_value)
    var = (
        auc_value * (1.0 - auc_value)
        + (n_active - 1) * (q1 - auc_value ** 2)
        + (n_decoy - 1) * (q2 - auc_value ** 2)
    ) / (n_active * n_decoy)
    if var <= 0:
        return float("nan")
    z = (auc_value - 0.5) / np.sqrt(var)
    return float(2.0 * norm.sf(abs(z)))


def _bh_adjust(pvalues: Sequence[float]) -> List[float]:
    values = np.asarray(pvalues, dtype=float)
    adjusted = np.full(values.shape, np.nan)
    finite = np.where(np.isfinite(values))[0]
    if len(finite) == 0:
        return adjusted.tolist()
    ordered = finite[np.argsort(values[finite])]
    m = len(ordered)
    running = 1.0
    for rank, idx in enumerate(ordered[::-1], start=1):
        position = m - rank + 1
        running = min(running, values[idx] * m / position)
        adjusted[idx] = min(1.0, running)
    return adjusted.tolist()


# ---------------------------------------------------------------------- driver
def _load_sampling_meta(scores_dir: Path) -> Dict[str, Dict[str, Any]]:
    """Per-target ``sampling_meta.json`` files under ``<scores_dir>/*/``.

    B's full-run layout writes ``<target>/sampling_meta.json`` with decoy ratio,
    sampling seed and subsample SHA256; these fill ``target_metadata`` fields when
    no explicit metadata manifest is present.
    """
    out: Dict[str, Dict[str, Any]] = {}
    for path in sorted(scores_dir.glob("*/sampling_meta.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except Exception:  # pragma: no cover - defensive
            continue
        target = str(payload.get("target") or path.parent.name)
        values: Dict[str, Any] = {}
        for field in METADATA_FIELDS:
            if field in payload and payload[field] is not None:
                values[field] = payload[field]
        out[target] = values
    return out


def run_tracka_framework(
    scores_dir: str | Path = DEFAULT_TRACKA_DIR,
    pattern: str = "*_scores.csv",
    target_map: Optional[Dict[str, str]] = None,
    n_boot: int = 2000,
    seed: int = 42,
    tag: str = "pilot",
    write: bool = True,
    metadata_path: Optional[str | Path] = None,
    box_clipped_path: Optional[str | Path] = None,
    allow_incomplete: bool = False,
    dry_run: bool = False,
    allow_partial_provisional: bool = False,
) -> Dict[str, Any]:
    """Run the Track A enrichment framework over a directory of score tables."""
    scores_dir = Path(scores_dir)
    files = sorted(scores_dir.glob(pattern))
    if not files:  # B full-run layout: <target>/<target>_scores.csv
        files = sorted(scores_dir.glob(f"*/{pattern}"))
    if not files:
        raise FileNotFoundError(f"no files matching {pattern!r} in {scores_dir}")
    metadata = load_target_metadata(metadata_path)
    sampling_meta: Dict[str, Dict[str, Any]] = {}
    if not metadata:
        sampling_meta = _load_sampling_meta(scores_dir)
        if sampling_meta:
            metadata = dict(sampling_meta)
    box_clipped = load_box_clipped_ligands(box_clipped_path)
    long_parts = [load_score_table(path, target_map) for path in files]
    frame = pd.concat(long_parts, ignore_index=True)
    # score tables may already carry a per-ligand box_clipped flag (B line); keep
    # it under a separate name so the A-line flag merge below does not collide
    if "box_clipped" in frame.columns:
        frame = frame.rename(columns={"box_clipped": "box_clipped_score"})
    # Amendment #10: exclude box-clipped ligands from the primary enrichment
    box_clipped_source = None
    box_clipped_match_mode = None
    box_a_flag_matched = 0
    if len(box_clipped):
        box_clipped_source = (
            str(box_clipped_path) if box_clipped_path is not None
            else next((str(p) for p in DEFAULT_BOX_CLIPPED_CANDIDATES if p.exists()), None)
        )
        merged_flags = frame.merge(
            box_clipped[["target", "ligand_id_norm", "box_clipped"]],
            on=["target", "ligand_id_norm"],
            how="left",
        )
        box_clipped_match_mode = "target+normalized_ligand_id"
        if merged_flags["box_clipped"].notna().sum() == 0:  # raw-id fallback
            merged_flags = frame.merge(
                box_clipped[["target", "ligand_id", "box_clipped"]].rename(
                    columns={"ligand_id": "ligand_id_box"}
                ),
                left_on=["target", "ligand_id"],
                right_on=["target", "ligand_id_box"],
                how="left",
            )
            box_clipped_match_mode = "target+raw_ligand_id"
        if merged_flags["box_clipped"].notna().sum() == 0:  # normalized-id-only fallback
            merged_flags = frame.merge(
                box_clipped[["ligand_id_norm", "box_clipped"]].drop_duplicates("ligand_id_norm"),
                on="ligand_id_norm",
                how="left",
            )
            box_clipped_match_mode = "normalized_ligand_id_only"
        frame = merged_flags
        score_flag = (
            frame["box_clipped_score"].map(_to_bool_value)
            if "box_clipped_score" in frame.columns
            else pd.Series(False, index=frame.index)
        )
        ai_flag = frame["box_clipped"].astype("boolean")
        frame["box_clipped"] = ai_flag.fillna(score_flag).astype(bool)
        box_a_flag_matched = int(ai_flag.notna().sum())
        # consistency between A-line audit flags and the score-table column
        both = ai_flag.notna() & score_flag.notna()
        frame["box_clipped_flag_mismatch"] = bool(both.any()) and bool((ai_flag[both] != score_flag[both]).any())
    else:
        if "box_clipped_score" in frame.columns:
            frame["box_clipped"] = frame["box_clipped_score"].map(_to_bool_value).fillna(False).astype(bool)
            box_clipped_source = "score_file"
        else:
            frame["box_clipped"] = False
    if "box_clipped_score" in frame.columns:
        frame = frame.drop(columns=["box_clipped_score"])
    if "box_clipped_flag_mismatch" not in frame.columns:
        frame["box_clipped_flag_mismatch"] = False
    # completeness guard: B writes scores at dock-stage completion, so a
    # mismatch against target_metadata.n_rows means the run is still in progress
    incomplete_targets: List[Dict[str, Any]] = []
    for target_name, target_meta in metadata.items():
        if "n_rows" not in target_meta:
            continue
        try:
            expected_rows = int(float(target_meta["n_rows"]))
        except (TypeError, ValueError):
            continue
        sub = frame[frame["target"].astype(str) == str(target_name)]
        if len(sub) == 0:
            continue
        prep_missing = (
            int(sub["n_prep_missing"].iloc[0])
            if "n_prep_missing" in sub.columns and pd.notna(sub["n_prep_missing"].iloc[0])
            else 0
        )
        other_non_ok = (
            int(sub["n_non_ok_other"].iloc[0])
            if "n_non_ok_other" in sub.columns and pd.notna(sub["n_non_ok_other"].iloc[0])
            else 0
        )
        rows_seen = int(len(sub) + prep_missing + other_non_ok)
        if rows_seen < expected_rows:
            incomplete_targets.append(
                {
                    "target": str(target_name),
                    "rows_seen": rows_seen,
                    "n_rows_expected": expected_rows,
                    "fraction_seen": round(float(rows_seen) / expected_rows, 4),
                }
            )
    if incomplete_targets and not allow_incomplete:
        detail = ", ".join(
            f"{item['target']} {item['rows_seen']}/{item['n_rows_expected']}"
            for item in incomplete_targets
        )
        raise IncompleteTrackADataError(
            "Track A score tables are incomplete (dock stage still running): "
            + detail
            + "; rerun with --allow-incomplete for a dry-run"
        )
    frame_all = frame.copy()
    frame = frame_all[~frame_all["box_clipped"]].reset_index(drop=True)
    # ------------------------------------------------ inference clusters (M-fix)
    # decoys are grouped by their matched active (matched_active_inchikey); each
    # active is its own cluster.  This is the primary CI/p inference unit.
    cluster_source_map: Dict[str, str] = {}
    frame["inference_cluster"] = "cluster:unknown"
    for target_name, sub in frame.groupby("target"):
        pool = _load_decoy_cluster_map(str(target_name))
        if pool is None:
            keys = [
                ("active:" + str(norm)) if int(label) == 1 else ("decoy:" + str(norm))
                for norm, label in zip(sub["ligand_id_norm"], sub["label"])
            ]
            cluster_source_map[str(target_name)] = "ligand_fallback"
        else:
            active_to_inchikey = pool["active_to_inchikey"]
            by_decoy = pool["by_decoy"]
            by_decoy_inchikey = pool.get("by_decoy_inchikey", {})
            keys = []
            for raw_id, norm, label in zip(sub["ligand_id"], sub["ligand_id_norm"], sub["label"]):
                norm = str(norm)
                if int(label) == 1:
                    inchikey = active_to_inchikey.get(_normalize_ligand_key(norm), "")
                    keys.append("active:" + str(inchikey) if inchikey else "active:" + norm)
                else:
                    # B score decoy ids are `dec_<inchikey>`; fall back to ChEMBL id
                    decoy_inchikey = str(raw_id).strip()
                    decoy_inchikey = decoy_inchikey[4:] if decoy_inchikey.lower().startswith("dec_") else decoy_inchikey
                    info = by_decoy_inchikey.get(decoy_inchikey)
                    if info is None:
                        info = by_decoy.get(norm)
                    inchikey = info.get("inchikey", "") if info else ""
                    if not inchikey and info:
                        inchikey = active_to_inchikey.get(info.get("active_chembl", ""), "")
                    keys.append("active:" + str(inchikey) if inchikey else "decoy:" + norm)
            cluster_source_map[str(target_name)] = pool["source"]
        frame.loc[sub.index, "inference_cluster"] = keys
    # ligand-property descriptors (MW / cLogP / HBD) from B manifests, for the
    # property-stratified and nearest-decoy sensitivity requested by Lead
    property_source_map: Dict[str, str] = {}
    property_parts: List[pd.DataFrame] = []
    for target_name in sorted(frame["target"].astype(str).unique()):
        candidates = [
            scores_dir / str(target_name) / "manifest.csv",
            scores_dir / f"{target_name}_manifest.csv",
        ]
        manifest_path = next((path for path in candidates if path.exists()), None)
        if manifest_path is None:
            continue
        try:
            props = load_ligand_properties(str(manifest_path))
        except Exception:
            continue
        props["target"] = str(target_name)
        property_source_map[str(target_name)] = str(manifest_path)
        property_parts.append(props)
    if property_parts:
        props_all = pd.concat(property_parts, ignore_index=True)
        frame = frame.merge(
            props_all[["target", "ligand_id", "mw", "logp", "hbd"]],
            on=["target", "ligand_id"],
            how="left",
        )
        if "mw" not in frame_all.columns:
            frame_all = frame_all.merge(
                props_all[["target", "ligand_id", "mw", "logp", "hbd"]],
                on=["target", "ligand_id"],
                how="left",
            )
    metadata_origin_map = (
        {target: "sampling_meta" for target in sampling_meta}
        if sampling_meta
        else {target: "manifest" for target in metadata}
    )
    # boron rescue surrogates (B -> C atom type; B-line rescues): main analysis
    # includes them; surrogate_sensitivity below excludes them for a sensitivity run
    surrogate_sources: Dict[str, str] = {}
    frame["surrogate"] = False
    frame_all["surrogate"] = False
    for target_name in sorted(frame["target"].astype(str).unique()):
        rescue_path = scores_dir / target_name / "rescue_boron_surrogate.csv"
        if not rescue_path.exists():
            continue
        try:
            rescue = pd.read_csv(rescue_path)
        except Exception:
            continue
        if "ligand_id" not in rescue.columns:
            continue
        rescue_ids = set(rescue["ligand_id"].astype(str).map(_normalize_ligand_key))
        mask = (frame["target"].astype(str) == target_name) & frame["ligand_id_norm"].astype(str).isin(rescue_ids)
        frame.loc[mask, "surrogate"] = True
        mask_all = (frame_all["target"].astype(str) == target_name) & frame_all["ligand_id_norm"].astype(str).isin(rescue_ids)
        frame_all.loc[mask_all, "surrogate"] = True
        surrogate_sources[target_name] = str(rescue_path)
    entries: List[Dict[str, Any]] = []
    for (target, method, condition), group in frame.groupby(["target", "method", "condition"]):
        metrics = compute_tracka_metrics(group["label"], group["score"])
        ligand_ci = bootstrap_enrichment_ci(group["label"], group["score"], n_boot=n_boot, seed=seed)
        cluster_bootstrap = cluster_bootstrap_enrichment_ci(
            group["label"], group["score"], group["inference_cluster"], n_boot=n_boot, seed=seed
        )
        cluster_ci = cluster_bootstrap.get("ci", {})
        cluster_values = group["inference_cluster"].astype(str)
        n_active_clusters = int(cluster_values[cluster_values.str.startswith("active:")].nunique())
        n_unmatched_decoy_clusters = int(cluster_values[cluster_values.str.startswith("decoy:")].nunique())
        metal = TARGET_METAL_DEPENDENCE.get(str(target), "unknown")
        target_meta = dict(metadata.get(str(target), {}))
        metadata_source = metadata_origin_map.get(str(target), "") if target_meta else ""
        for field in METADATA_FIELDS:
            if field in target_meta:
                continue
            if field in group.columns:
                values = group[field].dropna()
                if len(values):
                    value = values.iloc[0]
                    target_meta[field] = value.item() if hasattr(value, "item") else value
                    if not metadata_source:
                        metadata_source = "score_file"
        all_group = frame_all[
            (frame_all["target"].astype(str) == str(target))
            & (frame_all["method"].astype(str) == str(method))
            & (frame_all["condition"].astype(str) == str(condition))
        ]
        n_clipped = int(all_group["box_clipped"].sum())
        clipped_actives = int(((all_group["box_clipped"]) & (all_group["label"] == 1)).sum())
        clipped_decoys = int(((all_group["box_clipped"]) & (all_group["label"] == 0)).sum())
        n_act_before = int(all_group["label"].sum())
        # prep_missing audit (SMILES -> 3D failures are stored with status)
        prep_total = int(group["n_prep_missing"].iloc[0]) if "n_prep_missing" in group.columns else 0
        prep_actives = int(group["n_prep_missing_actives"].iloc[0]) if "n_prep_missing_actives" in group.columns else 0
        prep_decoys = int(group["n_prep_missing_decoys"].iloc[0]) if "n_prep_missing_decoys" in group.columns else 0
        other_non_ok = int(group["n_non_ok_other"].iloc[0]) if "n_non_ok_other" in group.columns else 0
        other_non_ok_actives = (
            int(group["n_non_ok_other_actives"].iloc[0])
            if "n_non_ok_other_actives" in group.columns else 0
        )
        other_non_ok_decoys = (
            int(group["n_non_ok_other_decoys"].iloc[0])
            if "n_non_ok_other_decoys" in group.columns else 0
        )
        n_actives_table_total = (
            int(group["n_actives_table_total"].iloc[0]) if "n_actives_table_total" in group.columns else n_act_before
        )
        prep_active_fraction = (
            float(prep_actives) / float(n_actives_table_total) if n_actives_table_total > 0 else float("nan")
        )
        n_ok_reuse = int(group["n_ok_reuse"].iloc[0]) if "n_ok_reuse" in group.columns else 0
        n_surrogate_actives = int((group["surrogate"].astype(bool) & (group["label"] == 1)).sum()) if "surrogate" in group.columns else 0
        n_surrogate_decoys = int((group["surrogate"].astype(bool) & (group["label"] == 0)).sum()) if "surrogate" in group.columns else 0
        panel_role = target_meta.get("panel_role")
        auxiliary = bool(panel_role is not None and str(panel_role).lower() != "main")
        property_result = None
        if all(column in group.columns for column in ("mw", "logp", "hbd")):
            props = group[["mw", "logp", "hbd"]]
            if props.notna().any().any():
                try:
                    property_result = property_sensitivity(
                        group["label"], group["score"], props, float(metrics["auc"])
                    )
                except Exception:
                    property_result = None
        n_temporal_missing = (
            int(group["temporal_split"].isna().sum()) if "temporal_split" in group.columns else None
        )
        n_temporal_missing_actives = (
            int((group["temporal_split"].isna() & (group["label"] == 1)).sum())
            if "temporal_split" in group.columns else None
        )
        observed_ratio = (
            float(metrics["n_decoy"]) / float(metrics["n_active"]) if metrics["n_active"] > 0 else float("nan")
        )
        meta_n_act = target_meta.get("n_act")
        meta_n_dec = target_meta.get("n_dec")
        planned_ratio = (
            float(meta_n_dec) / float(meta_n_act)
            if meta_n_act is not None and meta_n_dec is not None
            and pd.notna(meta_n_act) and pd.notna(meta_n_dec) and float(meta_n_act) > 0
            else float("nan")
        )
        entries.append(
            {
                "target": str(target),
                "method": str(method),
                "condition": str(condition),
                "metal_dependent": metal,
                "metal_stratum": metal_stratum(metal),
                "panel_role": panel_role,
                "auxiliary": auxiliary,
                **metrics,
                "n_box_clipped": n_clipped,
                "n_box_clipped_actives": clipped_actives,
                "n_box_clipped_decoys": clipped_decoys,
                "n_act_before_box_clipped_exclusion": n_act_before,
                "n_act_after_box_clipped_exclusion": int(metrics["n_active"]),
                "box_clipped": bool(n_clipped > 0),
                "box_clipped_source": box_clipped_source,
                "n_prep_missing": prep_total,
                "n_prep_missing_actives": prep_actives,
                "n_prep_missing_decoys": prep_decoys,
                "n_non_ok_other": other_non_ok,
                "n_non_ok_other_actives": other_non_ok_actives,
                "n_non_ok_other_decoys": other_non_ok_decoys,
                "non_ok_other_actives_fraction": (
                    float(other_non_ok_actives) / float(n_actives_table_total)
                    if n_actives_table_total > 0 else float("nan")
                ),
                "non_ok_other_actives_warning": bool(
                    n_actives_table_total > 0 and other_non_ok_actives / n_actives_table_total > 0.05
                ),
                "n_ok_reuse": n_ok_reuse,
                "n_surrogate_actives": n_surrogate_actives,
                "n_surrogate_decoys": n_surrogate_decoys,
                "surrogate_source": surrogate_sources.get(str(target)),
                "n_temporal_split_missing": n_temporal_missing,
                "n_temporal_split_missing_actives": n_temporal_missing_actives,
                "n_actives_table_total": n_actives_table_total,
                "n_act_after_prep_missing": int(n_actives_table_total - prep_actives),
                "prep_missing_active_fraction": prep_active_fraction,
                "prep_missing_actives_warning": bool(
                    np.isfinite(prep_active_fraction) and prep_active_fraction > 0.05
                ),
                "decoy_ratio": target_meta.get("decoy_ratio"),
                "decoy_ratio_observed": round(observed_ratio, 3) if np.isfinite(observed_ratio) else float("nan"),
                "decoy_ratio_planned": round(planned_ratio, 3) if np.isfinite(planned_ratio) else float("nan"),
                "n_act_metadata": meta_n_act,
                "n_dec_metadata": meta_n_dec,
                "sampling_seed": target_meta.get("sampling_seed"),
                "decoy_subsample_sha256": target_meta.get("decoy_subsample_sha256"),
                "exhaustiveness": target_meta.get("exhaustiveness"),
                "panel_role": target_meta.get("panel_role"),
                "box_size_A": target_meta.get("box_size_A"),
                "metadata_source": metadata_source or "observed",
                "inference_unit": (
                    "active-cluster bootstrap (decoys grouped by matched active)"
                    if cluster_source_map.get(str(target), "ligand_fallback") != "ligand_fallback"
                    else "ligand fallback (decoy matching map unavailable)"
                ),
                "cluster_source": cluster_source_map.get(str(target)),
                "n_inference_clusters": int(cluster_bootstrap.get("n_clusters", 0)),
                "n_active_clusters": n_active_clusters,
                "n_unmatched_decoy_clusters": n_unmatched_decoy_clusters,
                "auc_ci_lo": cluster_ci.get("auc", (float("nan"), float("nan")))[0],
                "auc_ci_hi": cluster_ci.get("auc", (float("nan"), float("nan")))[1],
                "ef5_ci_lo": cluster_ci.get("ef5_percent", (float("nan"), float("nan")))[0],
                "ef5_ci_hi": cluster_ci.get("ef5_percent", (float("nan"), float("nan")))[1],
                "bedroc_ci_lo": cluster_ci.get("bedroc20", (float("nan"), float("nan")))[0],
                "bedroc_ci_hi": cluster_ci.get("bedroc20", (float("nan"), float("nan")))[1],
                "ef1_ci_lo": cluster_ci.get("ef1_percent", (float("nan"), float("nan")))[0],
                "ef1_ci_hi": cluster_ci.get("ef1_percent", (float("nan"), float("nan")))[1],
                # ligand-level bootstrap retained as sensitivity only
                "auc_ci_ligand_lo": ligand_ci.get("auc", (float("nan"), float("nan")))[0],
                "auc_ci_ligand_hi": ligand_ci.get("auc", (float("nan"), float("nan")))[1],
                "ef5_ci_ligand_lo": ligand_ci.get("ef5_percent", (float("nan"), float("nan")))[0],
                "ef5_ci_ligand_hi": ligand_ci.get("ef5_percent", (float("nan"), float("nan")))[1],
                "bedroc_ci_ligand_lo": ligand_ci.get("bedroc20", (float("nan"), float("nan")))[0],
                "bedroc_ci_ligand_hi": ligand_ci.get("bedroc20", (float("nan"), float("nan")))[1],
                "p_auc": _hanley_mcneil_p(metrics["auc"], metrics["n_active"], metrics["n_decoy"]),
                "p_auc_interpretation": (
                    "descriptive only — ligand-level Hanley-McNeil; invalid under "
                    "matched-decoy clustering (do not use for inference)"
                ),
                "auc_property_stratified": (
                    property_result["auc_property_stratified"] if property_result else None
                ),
                "auc_property_by_property": (
                    property_result["by_property"] if property_result else None
                ),
                "auc_property_matched": (
                    property_result["auc_property_matched"] if property_result else None
                ),
                "property_delta_vs_raw": (
                    property_result["delta_vs_raw"] if property_result else None
                ),
                "property_sensitivity_method": (
                    property_result["method"] if property_result else None
                ),
                "property_source": property_source_map.get(str(target)),
                "n_input_rows": int(len(group)),
                "source_files": sorted({str(g["source_file"].iloc[0]) for _, g in group.groupby("source_file")}),
            }
        )
    for method in {entry["method"] for entry in entries}:
        mask = [i for i, entry in enumerate(entries) if entry["method"] == method]
        adjusted = _bh_adjust([entries[i]["p_auc"] for i in mask])
        for i, p_adj in zip(mask, adjusted):
            entries[i]["p_auc_bh"] = p_adj
    summary = frame.shape[0]
    result: Dict[str, Any] = {
        "tag": tag,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "scores_dir": str(scores_dir),
        "n_files": len(files),
        "n_rows_long": int(summary),
        "incomplete": bool(incomplete_targets),
        "incomplete_targets": incomplete_targets,
        "dry_run": bool(dry_run or (allow_incomplete and bool(incomplete_targets))),
        "entries": entries,
        "by_metal_stratum": {},
        "notes": [
            "EF1% is reported only when n_active >= 150 (R1); otherwise ef1_percent=None and primary metric is EF5%.",
            "EF5% is the enrichment factor at the top 5% (ratio; 2.0 = 2x).",
            "BEDROC uses alpha=20; AUROC bootstrap is a ligand-level percentile bootstrap (2.5/97.5).",
            "p_auc is Hanley-McNeil (AUROC vs 0.5), BH-adjusted within method across target entries.",
            "Metal strata: LpxC=Zn, GyrA/GyrB=Mg, DHFR/FabI/LepB/TEM-1=none (evaluation.TARGET_METAL_DEPENDENCE).",
        ],
    }
    for stratum in ("metal", "non-metal", "unknown"):
        subset = [entry for entry in entries if entry["metal_stratum"] == stratum]
        if not subset:
            continue
        result["by_metal_stratum"][stratum] = {
            "n_entries": len(subset),
            "targets": sorted({entry["target"] for entry in subset}),
            "mean_auc": float(np.nanmean([entry["auc"] for entry in subset])),
            "mean_ef5": float(np.nanmean([entry["ef5_percent"] for entry in subset])),
        }
    # Amendment #10 box-clipping audit (per target: pre/post exclusion counts)
    clip_rows: List[Dict[str, Any]] = []
    for target, sub in frame_all.groupby("target"):
        clipped = sub["box_clipped"].astype(bool)
        clip_rows.append(
            {
                "target": str(target),
                "n_ligands": int(len(sub)),
                "n_box_clipped": int(clipped.sum()),
                "n_box_clipped_actives": int((clipped & (sub["label"] == 1)).sum()),
                "n_box_clipped_decoys": int((clipped & (sub["label"] == 0)).sum()),
                "n_act_before": int(sub["label"].sum()),
                "n_act_after": int(((~clipped) & (sub["label"] == 1)).sum()),
                "box_size_A": (metadata.get(str(target), {}) or {}).get("box_size_A"),
            }
        )
    # boron-surrogate sensitivity: main analysis includes B->C rescues; this
    # block recomputes metrics excluding the rescued ligands and reports deltas
    surrogate_sensitivity: Dict[str, Any] = {}
    for target_name, source in surrogate_sources.items():
        sub = frame[frame["target"].astype(str) == target_name]
        entry = next((entry for entry in entries if entry["target"] == str(target_name)), None)
        if entry is None or len(sub) == 0:
            continue
        sub_no = sub[~sub["surrogate"].astype(bool)]
        if len(sub_no) == 0:
            continue
        metrics_no = compute_tracka_metrics(sub_no["label"], sub_no["score"])
        cluster_no = cluster_bootstrap_enrichment_ci(
            sub_no["label"], sub_no["score"], sub_no["inference_cluster"], n_boot=n_boot, seed=seed
        )
        ci_no = cluster_no.get("ci", {})
        surrogate_sensitivity[str(target_name)] = {
            "source": source,
            "n_surrogate_actives": entry["n_surrogate_actives"],
            "n_surrogate_decoys": entry["n_surrogate_decoys"],
            "with_surrogate": {
                "auc": entry["auc"],
                "ef5_percent": entry["ef5_percent"],
                "bedroc20": entry["bedroc20"],
                "n_active": entry["n_active"],
                "n_decoy": entry["n_decoy"],
            },
            "without_surrogate": {
                "auc": metrics_no["auc"],
                "auc_ci_cluster": [
                    ci_no.get("auc", (float("nan"), float("nan")))[0],
                    ci_no.get("auc", (float("nan"), float("nan")))[1],
                ],
                "ef5_percent": metrics_no["ef5_percent"],
                "ef5_ci_cluster": [
                    ci_no.get("ef5_percent", (float("nan"), float("nan")))[0],
                    ci_no.get("ef5_percent", (float("nan"), float("nan")))[1],
                ],
                "bedroc20": metrics_no["bedroc20"],
                "n_active": metrics_no["n_active"],
                "n_decoy": metrics_no["n_decoy"],
            },
            "delta_auc": float(metrics_no["auc"] - entry["auc"]),
            "delta_ef5": float(metrics_no["ef5_percent"] - entry["ef5_percent"]),
            "n_act_main": entry["n_active"],
            "n_act_sensitivity": metrics_no["n_active"],
            "pre_rescue_formal": {
                "auc": metrics_no["auc"],
                "ef5_percent": metrics_no["ef5_percent"],
                "n_active": metrics_no["n_active"],
                "note": "pre-rescue run had only non-surrogate actives; this equals the sensitivity arm",
            },
            "delta_auc_vs_pre_rescue": float(metrics_no["auc"] - entry["auc"]),
            "delta_ef5_vs_pre_rescue": float(metrics_no["ef5_percent"] - entry["ef5_percent"]),
            "note": (
                "main analysis includes B->C isosteric surrogates; this sensitivity excludes "
                "rescued ligands; difference reported for SI"
            ),
        }
    result["surrogate_sensitivity"] = surrogate_sensitivity
    result["surrogate_note"] = (
        f"{sum(len(pd.read_csv(path)) if path.endswith('.csv') else 0 for path in surrogate_sources.values())} "
        "boron-containing ligands were rescued with a B->C isosteric surrogate (Vina atom-type "
        "limitation); main analysis includes them and surrogate_sensitivity reports the "
        "exclude-surrogate sensitivity."
        if surrogate_sources else
        "no boron-surrogate ligands in the current score tables"
    )
    result["box_clipped_audit"] = {
        "source": box_clipped_source,
        "match_mode": box_clipped_match_mode,
        "n_flagged_rows": int(frame_all["box_clipped"].sum()),
        "n_a_file_rows": int(len(box_clipped)),
        "n_a_flags_matched_to_scores": int(box_a_flag_matched),
        "flag_match_rate": (
            float(box_a_flag_matched) / len(box_clipped) if len(box_clipped) else float("nan")
        ),
        "unmatched_flag_ligands": [
            f"{row['target']}|{row['ligand_id']}"
            for _, row in box_clipped.iterrows()
            if (
                str(row["target"]), str(row["ligand_id_norm"])
            ) not in set(
                zip(frame_all["target"].astype(str), frame_all["ligand_id_norm"].astype(str))
            )
        ][:20],
        "rule": "max_dim_A > box_size_A -> box_clipped=True -> excluded from primary enrichment (Amendment #10)",
        "by_target": clip_rows,
        "note": (
            "clipped ligands are not docked and cannot enter the primary enrichment denominators; "
            "SI reports n_box_clipped / n_box_clipped_actives / n_act_after"
        ),
    }
    # prep_missing audit (data-processing rule; Lead 2026-11-17)
    prep_rows: List[Dict[str, Any]] = []
    for target, sub in frame_all.groupby("target"):
        def _first(column: str) -> int:
            if column not in sub.columns or len(sub) == 0:
                return 0
            value = sub[column].iloc[0]
            return int(value) if pd.notna(value) else 0

        prep_total = _first("n_prep_missing")
        prep_actives = _first("n_prep_missing_actives")
        prep_decoys = _first("n_prep_missing_decoys")
        table_actives = _first("n_actives_table_total")
        fraction = float(prep_actives) / float(table_actives) if table_actives > 0 else float("nan")
        meta_value = (metadata.get(str(target), {}) or {}).get("n_prep_missing")
        meta_rows = (metadata.get(str(target), {}) or {}).get("n_rows")
        rows_seen = int(len(sub) + prep_total + _first("n_non_ok_other"))
        table_incomplete = bool(
            meta_rows is not None and pd.notna(meta_rows) and rows_seen < int(meta_rows)
        )
        mismatched = bool(
            meta_value is not None
            and pd.notna(meta_value)
            and not table_incomplete
            and int(meta_value) != prep_total
        )
        prep_rows.append(
            {
                "target": str(target),
                "n_prep_missing": prep_total,
                "n_prep_missing_actives": prep_actives,
                "n_prep_missing_decoys": prep_decoys,
                "n_actives_table_total": table_actives,
                "n_act_after_prep_missing": int(table_actives - prep_actives),
                "prep_missing_active_fraction": fraction,
                "prep_missing_actives_warning": bool(np.isfinite(fraction) and fraction > 0.05),
                "metadata_n_prep_missing": meta_value,
                "metadata_n_rows": meta_rows,
                "rows_seen_in_score_table": rows_seen,
                "table_incomplete": table_incomplete,
                "metadata_mismatch": mismatched,
                "metadata_mismatch_pending_table_completion": bool(
                    meta_value is not None and pd.notna(meta_value) and table_incomplete
                ),
            }
        )
    result["prep_missing_audit"] = {
        "rule": (
            "SMILES->3D failures (status=prep_missing) are excluded from enrichment "
            "denominators and reported separately; no result-based filtering"
        ),
        "by_target": prep_rows,
        "warning_targets": [row["target"] for row in prep_rows if row["prep_missing_actives_warning"]],
        "note": "flag raised to Lead when prep_missing actives exceed 5% of a target's actives",
    }
    # other non-ok statuses (docking fails, box_clipped_not_docked, ...) — data quality
    dq_rows: List[Dict[str, Any]] = []
    for target, sub in frame_all.groupby("target"):
        def _dq_first(column: str) -> int:
            if column not in sub.columns or len(sub) == 0:
                return 0
            value = sub[column].iloc[0]
            return int(value) if pd.notna(value) else 0

        other_act = _dq_first("n_non_ok_other_actives")
        other_dec = _dq_first("n_non_ok_other_decoys")
        table_act = _dq_first("n_actives_table_total")
        dq_rows.append(
            {
                "target": str(target),
                "n_non_ok_other_actives": other_act,
                "n_non_ok_other_decoys": other_dec,
                "n_actives_table_total": table_act,
                "non_ok_other_actives_fraction": (
                    float(other_act) / table_act if table_act > 0 else float("nan")
                ),
                "warning": bool(table_act > 0 and other_act / table_act > 0.05),
            }
        )
    result["data_quality_audit"] = {
        "by_target": dq_rows,
        "warning_targets": [row["target"] for row in dq_rows if row["warning"]],
        "note": (
            "non-ok, non-prep_missing statuses (e.g. 'fail', 'box_clipped_not_docked') are "
            "excluded from enrichment denominators; targets with >5% actives affected are "
            "flagged for Lead review"
        ),
    }
    # temporal split coverage audit (Amendment #9 prep; B merge may leave actives empty)
    temporal_rows: List[Dict[str, Any]] = []
    for target, sub in frame_all.groupby("target"):
        if "temporal_split" not in sub.columns:
            continue
        missing = sub["temporal_split"].isna()
        temporal_rows.append(
            {
                "target": str(target),
                "n_rows": int(len(sub)),
                "n_missing": int(missing.sum()),
                "n_missing_actives": int((missing & (sub["label"] == 1)).sum()),
                "n_missing_decoys": int((missing & (sub["label"] == 0)).sum()),
                "known_values": sorted(
                    {str(value) for value in sub["temporal_split"].dropna().unique()}
                )[:5],
            }
        )
    result["temporal_split_audit"] = {
        "by_target": temporal_rows,
        "warning": bool(any(row["n_missing_actives"] > 0 for row in temporal_rows)),
        "note": (
            "Amendment #9 novel-scaffold analysis needs temporal_split for actives; "
            "if actives are missing, B must re-merge by normalised ChEMBL id"
        ),
    }
    # observed decoy ratio (matched decoys can fall short of the 20:1 request)
    ratios = [entry["decoy_ratio_observed"] for entry in entries if np.isfinite(entry["decoy_ratio_observed"])]
    if ratios:
        planned = [
            entry["decoy_ratio_planned"]
            for entry in entries
            if np.isfinite(entry.get("decoy_ratio_planned", float("nan")))
        ]
        result["decoy_ratio_summary"] = {
            "min": float(min(ratios)),
            "max": float(max(ratios)),
            "median": float(np.median(ratios)),
            "planned_min": float(min(planned)) if planned else None,
            "planned_max": float(max(planned)) if planned else None,
            "by_entry": {
                f"{entry['target']}|{entry['condition']}|{entry['method']}": {
                    "observed": entry["decoy_ratio_observed"],
                    "planned": entry.get("decoy_ratio_planned"),
                }
                for entry in entries
            },
            "note": (
                f"up to 20 matched decoys per active (deterministic seed 42; observed "
                f"{min(ratios):.2f}-{max(ratios):.2f}, median {np.median(ratios):.2f}); "
                "per-target decoy_ratio_observed retained; planned values include prep-failed/"
                "box-clipped actives and are not the analysis denominator; EF/AUC computed on "
                "actual N"
            ),
        }
        result["decoy_ratio_definition"] = (
            f"up to 20 matched decoys per active (deterministic seed 42; observed "
            f"{min(ratios):.2f}-{max(ratios):.2f}, median {np.median(ratios):.2f})"
        )
    primary_entries = [entry for entry in entries if not entry.get("auxiliary")]
    auxiliary_entries = [entry for entry in entries if entry.get("auxiliary")]
    metal_primary = [entry for entry in primary_entries if entry["metal_stratum"] == "metal"]
    nonmetal_primary = [entry for entry in primary_entries if entry["metal_stratum"] == "non-metal"]
    inverted_metal = [
        entry for entry in metal_primary
        if np.isfinite(entry["auc"]) and entry["auc"] < 0.5
    ]
    strong_metal = [
        entry for entry in metal_primary
        if np.isfinite(entry["auc"]) and entry["auc"] >= 0.6
    ]
    chemotype_evidence = _load_chemotype_evidence()
    decoy_overlap_evidence = _load_decoy_overlap_evidence()
    chemotype_by_target = {
        str(row["target"]): row for row in (chemotype_evidence or [])
    }
    result["interpretation_notes"] = []
    if primary_entries:
        inverted_text = ", ".join(f"{entry['target']} ({entry['auc']:.3f})" for entry in inverted_metal) or "none"
        strong_text = ", ".join(f"{entry['target']} ({entry['auc']:.3f})" for entry in strong_metal) or "none"
        nonmetal_text = ", ".join(
            f"{entry['target']} ({entry['auc']:.3f})" for entry in nonmetal_primary
        ) or "none yet"
        result["interpretation_notes"].append(
            "Scoring performance is TARGET-SPECIFIC, not uniformly metal-determined: "
            f"metal-dependent targets split into strong ({strong_text}) and inverted ({inverted_text}) "
            f"outcomes; non-metal control(s): {nonmetal_text}. Candidate explanations: metal "
            "coordination, ligand property/chemotype concentration, activity strength, assay "
            "heterogeneity — no causal conclusion is drawn."
        )
    lpxc_entries = [entry for entry in primary_entries if entry["target"] == "LpxC_Ec"]
    if lpxc_entries:
        result["interpretation_notes"].append(
            "LpxC_Ec ligand-property contrast (A-line audit): actives HBD 4 vs decoys 3; "
            "median cLogP 1.64 vs 2.05 — ligand polarity differences may contribute to the "
            "inverted enrichment, alongside metal coordination."
        )
        if all(column in frame_all.columns for column in ("mw", "logp", "hbd")):
            lpxc_frame = frame_all[frame_all["target"].astype(str) == "LpxC_Ec"]
            act = lpxc_frame[lpxc_frame["label"] == 1]
            dec = lpxc_frame[lpxc_frame["label"] == 0]
            if len(act) and len(dec):
                result["interpretation_notes"].append(
                    "LpxC_Ec property medians recomputed from manifests (actives vs decoys): "
                    f"HBD {act['hbd'].median():.1f} vs {dec['hbd'].median():.1f}, "
                    f"cLogP {act['logp'].median():.2f} vs {dec['logp'].median():.2f}, "
                    f"MW {act['mw'].median():.1f} vs {dec['mw'].median():.1f}."
                )
    gyra_rows = [entry for entry in auxiliary_entries if entry["target"] == "GyrA_Ec"]
    if gyra_rows:
        gyra = gyra_rows[0]
        gyra_chem = chemotype_by_target.get("GyrA_Ec", {})
        result["interpretation_notes"].append(
            "GyrA_Ec is low_confidence_auxiliary (chemotype concentration: max BM scaffold "
            f"{gyra_chem.get('max_scaffold_frac')}, n_documents {gyra_chem.get('n_documents')}, "
            f"pChEMBL median {gyra_chem.get('pchembl_median')}, observed ratio "
            f"{gyra['decoy_ratio_observed']}); its AUC {gyra['auc']:.3f} is a chemotype/data "
            "artefact candidate and is excluded from all primary/stratum aggregates (SI only)."
        )
    result["interpretation_notes"].append(
        "Stratum aggregates are heterogeneous and descriptive only (primary targets: "
        + ", ".join(entry["target"] for entry in primary_entries)
        + "); per-target reporting is primary. GyrA_Ec/LepB_Ec auxiliary targets are excluded "
          "from aggregates."
    )
    # LpxC_Pa vs LpxC_Ec narrative is added after the series-exclusion sensitivity
    # is computed below (it requires the exact delta values).
    if nonmetal_primary and any(entry["auc"] < 0.5 for entry in nonmetal_primary if np.isfinite(entry["auc"])):
        result["interpretation_notes"].append(
            "A non-metal control also shows AUROC < 0.5 — flag to Lead; this weakens any "
            "metal-specific interpretation and should be treated as a pipeline/assay check."
        )
    result["chemotype_evidence"] = chemotype_evidence
    result["decoy_overlap_evidence"] = decoy_overlap_evidence
    # -------- supplemental sensitivity analyses (Lead 2026-11-18)
    supplementary: Dict[str, Any] = {}
    supplementary["property_sensitivity"] = {
        entry["target"]: {
            "raw_auc": entry["auc"],
            "auc_property_stratified": entry.get("auc_property_stratified"),
            "auc_property_matched": entry.get("auc_property_matched"),
            "delta_vs_raw": entry.get("property_delta_vs_raw"),
            "panel_role": entry.get("panel_role"),
        }
        for entry in entries
        if entry.get("auc_property_stratified") is not None
    }
    # (b) LpxC_Pa series exclusion: remove largest BM scaffold / ECFP4 cluster
    pa_entry = next((entry for entry in entries if entry["target"] == "LpxC_Pa"), None)
    if pa_entry is not None and property_source_map.get("LpxC_Pa"):
        pa_manifest = pd.read_csv(property_source_map["LpxC_Pa"])
        pa_actives = pa_manifest[pa_manifest["label"] == 1].reset_index(drop=True)
        pa_frame = frame_all[frame_all["target"].astype(str) == "LpxC_Pa"]
        pa_act = pa_frame[pa_frame["label"] == 1]
        pa_dec = pa_frame[pa_frame["label"] == 0]
        smiles = pa_actives["smiles"].astype(str).tolist()
        scaffold = largest_scaffold_group(smiles)
        cluster = largest_ecfp_cluster(smiles)
        def _pa_excl_metrics(excluded: set[str]) -> Dict[str, Any]:
            keep = ~pa_act["ligand_id"].astype(str).isin(excluded)
            y = np.concatenate([np.ones(int(keep.sum())), np.zeros(len(pa_dec))])
            s = np.concatenate([pa_act.loc[keep, "score"].to_numpy(), pa_dec["score"].to_numpy()])
            return compute_tracka_metrics(y, s)
        scaf_ids = set(pa_actives.iloc[scaffold["indices"]]["ligand_id"].astype(str)) if scaffold["indices"] else set()
        clust_ids = set(pa_actives.iloc[cluster["indices"]]["ligand_id"].astype(str)) if cluster["indices"] else set()
        scaf_metrics = _pa_excl_metrics(scaf_ids) if scaf_ids else None
        clust_metrics = _pa_excl_metrics(clust_ids) if clust_ids else None
        supplementary["pa_series_exclusion"] = {
            "raw": {"auc": pa_entry["auc"], "ef5_percent": pa_entry["ef5_percent"], "n_active": pa_entry["n_active"]},
            "largest_bm_scaffold": {
                "scaffold": scaffold.get("scaffold"),
                "n_excluded": int(len(scaf_ids)),
                "auc_after": scaf_metrics["auc"] if scaf_metrics else None,
                "ef5_after": scaf_metrics["ef5_percent"] if scaf_metrics else None,
                "delta_auc_vs_raw": (scaf_metrics["auc"] - pa_entry["auc"]) if scaf_metrics else None,
            },
            "largest_ecfp_cluster": {
                "n_excluded": int(len(clust_ids)),
                "auc_after": clust_metrics["auc"] if clust_metrics else None,
                "ef5_after": clust_metrics["ef5_percent"] if clust_metrics else None,
                "delta_auc_vs_raw": (clust_metrics["auc"] - pa_entry["auc"]) if clust_metrics else None,
            },
            "note": (
                "Exclusion sensitivity only — the main panel keeps all actives; this checks "
                "whether Pa enrichment is driven by its dominant chemotype series."
            ),
        }
        pa_block = supplementary["pa_series_exclusion"]
        pa_interpretation = None
        if scaf_metrics is not None and clust_metrics is not None:
            pa_interpretation = (
                "LpxC_Pa series exclusion: removing the largest BM scaffold "
                f"({pa_block['largest_bm_scaffold']['n_excluded']} actives) gives AUC "
                f"{scaf_metrics['auc']:.3f} (Δ{scaf_metrics['auc'] - pa_entry['auc']:+.3f}) and "
                "removing the largest ECFP4 cluster "
                f"({pa_block['largest_ecfp_cluster']['n_excluded']} actives) gives AUC "
                f"{clust_metrics['auc']:.3f} (Δ{clust_metrics['auc'] - pa_entry['auc']:+.3f}); "
                f"EF5 does not decrease (raw {pa_entry['ef5_percent']:.3f} → "
                f"{scaf_metrics['ef5_percent']:.3f}/{clust_metrics['ef5_percent']:.3f}). "
                "This does NOT support a chemotype-concentration explanation for Pa's high "
                "AUC; the LpxC_Ec vs LpxC_Pa difference is retained as a target/receptor/"
                "activity-space difference. Caveat: the two LpxC share 4,739 s1 decoys "
                "(60.3% of the Ec pool) and 48 shared actives; the cross-species axis uses "
                "the paired comparison on 951 shared scored decoys."
            )
            result["interpretation_notes"].append(pa_interpretation)
        pa_block["interpretation"] = pa_interpretation
    # (c) LpxC_Ec vs LpxC_Pa shared-decoy paired analysis
    ec_frame = frame_all[frame_all["target"].astype(str) == "LpxC_Ec"]
    pa_frame_all = frame_all[frame_all["target"].astype(str) == "LpxC_Pa"]
    if len(ec_frame) and len(pa_frame_all):
        ec_decoys = ec_frame[ec_frame["label"] == 0]
        pa_decoys = pa_frame_all[pa_frame_all["label"] == 0]
        shared_ids = sorted(set(ec_decoys["ligand_id_norm"].astype(str)) & set(pa_decoys["ligand_id_norm"].astype(str)))
        if shared_ids:
            ec_shared = ec_decoys[ec_decoys["ligand_id_norm"].astype(str).isin(shared_ids)]
            pa_shared = pa_decoys[pa_decoys["ligand_id_norm"].astype(str).isin(shared_ids)]
            ec_actives = ec_frame[ec_frame["label"] == 1]
            pa_actives_frame = pa_frame_all[pa_frame_all["label"] == 1]
            auc_ec = auroc(
                np.concatenate([np.ones(len(ec_actives)), np.zeros(len(ec_shared))]),
                np.concatenate([ec_actives["score"], ec_shared["score"]]),
            )
            auc_pa = auroc(
                np.concatenate([np.ones(len(pa_actives_frame)), np.zeros(len(pa_shared))]),
                np.concatenate([pa_actives_frame["score"], pa_shared["score"]]),
            )
            supplementary["shared_decoy_paired"] = {
                "n_shared_scored_decoys": int(len(shared_ids)),
                "n_ec_scored_decoys": int(len(ec_decoys)),
                "n_pa_scored_decoys": int(len(pa_decoys)),
                "auc_ec_shared": float(auc_ec),
                "auc_pa_shared": float(auc_pa),
                "delta_pa_minus_ec": float(auc_pa - auc_ec),
                "note": (
                    "AUROC on the same scored decoy IDs for both targets; decoy sets are not "
                    "identical to the full panel and this is a sensitivity check, not a replacement."
                ),
            }
    result["supplementary_sensitivities"] = supplementary
    # Lead-requested alias + data-quality note for the boron rescues
    supplementary["boron_surrogate_sensitivity"] = surrogate_sensitivity
    if surrogate_sensitivity:
        dq = result.setdefault("data_quality_audit", {})
        dq["boron_surrogate"] = {
            "by_target": surrogate_sensitivity,
            "note": "main analysis includes B->C isosteric surrogates; exclude-surrogate arm is the SI sensitivity",
        }
        pieces = []
        for target_name, info in surrogate_sensitivity.items():
            pieces.append(
                f"{target_name}: n_act {info['n_act_main']}→{info['n_act_sensitivity']} when "
                f"boron surrogates excluded; ΔAUC(main−exclude)={info['delta_auc']:+.4f}"
            )
        dq["notes"] = (
            "TEM1/LepB 硼酸类 actives 因 Vina 1.2.5 无硼原子类型失败，经 B→C isosteric "
            "surrogate 补救；主分析纳入，敏感性排除后结果变化：" + "；".join(pieces) + "。"
        )
    # Methods evidence for the SI/Methods additions requested by Lead
    shared_intersection = None
    for row in (decoy_overlap_evidence or []):
        if str(row.get("pair")) == "LpxC_Ec_vs_LpxC_Pa":
            shared_intersection = int(row.get("n_intersection"))
    result["methods_evidence"] = {
        "chemotype_internal_maxT_median": {
            str(row["target"]): row.get("internal_maxT_median") for row in (chemotype_evidence or [])
        },
        "top5_scaffold_frac": {
            str(row["target"]): row.get("top5_scaffold_frac") for row in (chemotype_evidence or [])
        },
        "top5_ecfp_cluster_frac": {
            str(row["target"]): row.get("top5_ecfp_cluster_frac") for row in (chemotype_evidence or [])
        },
        "shared_decoy_intersection_lpxc": shared_intersection,
        "si_chemotype_table": "data/trackA/si_chemotype_table.csv",
        "pa_exclude_top_series": "data/trackA/LpxC_Pa_sensitivity_exclude_top.csv",
        "note": "Methods/SI additions: internal max Tanimoto median, top-5 scaffold/ECFP cluster shares, shared decoy intersection",
    }
    if surrogate_sensitivity:
        result["methods_evidence"]["boron_surrogate_note"] = (
            "boron-containing ligands (TEM1 23, LepB 2) failed Vina 1.2.5 atom-type B and were "
            "rescued with a B->C isosteric surrogate; main analysis includes them and "
            "supplementary_sensitivities.boron_surrogate_sensitivity reports the exclude-surrogate arm."
        )
    novel = _load_novel_scaffold_splits()
    if len(novel):
        novel_rows: List[Dict[str, Any]] = []
        for target in sorted(novel["target"].unique()):
            novel_target = novel[novel["target"] == str(target)]
            sub = frame_all[(frame_all["target"].astype(str) == str(target)) & (frame_all["label"] == 1)]
            test_actives = sub[sub["temporal_split"].astype(str).str.lower() == "test"] if "temporal_split" in sub.columns else sub.iloc[0:0]
            matched = int(test_actives["ligand_id_norm"].isin(novel_target["ligand_id_norm"]).sum())
            novel_rows.append(
                {
                    "target": str(target),
                    "n_novel_split_rows": int(len(novel_target)),
                    "n_test_actives_in_score_table": int(len(test_actives)),
                    "n_matched_to_score_table": matched,
                    "linkage_fraction": (float(matched) / len(test_actives)) if len(test_actives) else float("nan"),
                }
            )
        result["novel_scaffold_audit"] = {
            "by_target": novel_rows,
            "n_novel_total_in_splits": int(len(novel)),
            "n_matched_total": int(sum(row["n_matched_to_score_table"] for row in novel_rows)),
            "note": (
                "Amendment #9 strict-transfer subset; matching is by normalised ChEMBL id "
                "against score-table actives with temporal_split='test'"
            ),
        }
    done_targets = sorted({str(entry["target"]) for entry in entries})
    missing_targets = [target for target in EXPECTED_TRACKA_PANEL if target not in done_targets]
    panel_complete = len(missing_targets) == 0
    result["panel"] = {
        "expected_targets": list(EXPECTED_TRACKA_PANEL),
        "n_expected_targets": int(len(EXPECTED_TRACKA_PANEL)),
        "done_targets": done_targets,
        "n_done_targets": int(len(done_targets)),
        "missing_targets": missing_targets,
        "panel_complete": bool(panel_complete),
        "note": (
            "formal trackA_enrichment_full.json is written only when panel_complete=true; "
            "otherwise the output is forced to a _dryrun file"
        ),
    }
    if not panel_complete:
        result["dry_run"] = True
        result["forced_dry_run"] = True
        result["forced_dry_run_reason"] = "panel incomplete: missing " + ", ".join(missing_targets)
    result["panel_complete"] = bool(panel_complete)
    result["missing_targets"] = list(missing_targets)
    # explicit scope: 9-target panel vs 7 primary (+2 auxiliary), and active counts
    panel_actives = 0
    primary_actives = 0
    auxiliary_targets = {str(entry["target"]) for entry in auxiliary_entries}
    for target_name in EXPECTED_TRACKA_PANEL:
        n_act = (metadata.get(str(target_name), {}) or {}).get("n_act")
        if n_act is None or pd.isna(n_act):
            continue
        n_act = int(float(n_act))
        panel_actives += n_act
        if str(target_name) not in auxiliary_targets:
            primary_actives += n_act
    result["panel_scope"] = {
        "n_panel_targets": int(len(EXPECTED_TRACKA_PANEL)),
        "n_primary_targets": int(len(primary_entries)),
        "n_auxiliary": int(len(auxiliary_entries)),
        "n_actives_panel": int(panel_actives),
        "n_actives_primary": int(primary_actives),
        "primary_targets": [entry["target"] for entry in primary_entries],
        "auxiliary_targets": [entry["target"] for entry in auxiliary_entries],
        "note": (
            "main tables/figures/aggregates use the 7 primary targets only; GyrA_Ec and "
            "LepB_Ec are auxiliary and reported separately"
        ),
    }
    result["primary_caliber"] = "includes_boron_surrogate"
    panel_ready_signal = load_panel_ready_signal()
    if panel_ready_signal is not None:
        result["panel_ready_signal"] = panel_ready_signal
        result["panel_ready_signal_present"] = True
    else:
        result["panel_ready_signal_present"] = False
    # small-print inference summary
    result["inference"] = {
        "primary_unit": "active-level cluster bootstrap (decoys grouped by matched active; each active its own cluster)",
        "sensitivity_unit": "ligand-level bootstrap",
        "cluster_sources": cluster_source_map,
        "p_value_note": (
            "p_auc / p_auc_bh are ligand-level Hanley-McNeil descriptive only — "
            "invalid under matched-decoy clustering; do not use for inference"
        ),
    }
    # metal vs non-metal strata (used for the internal provisional read)
    # GyrA_Ec / LepB_Ec are auxiliary (low-confidence) and excluded from any
    # stratum aggregate; per-target values remain in their own blocks.
    primary_entries = [entry for entry in entries if not entry.get("auxiliary")]
    auxiliary_entries = [entry for entry in entries if entry.get("auxiliary")]
    strata: Dict[str, Any] = {}
    for stratum in ("metal", "non-metal", "unknown"):
        subset = [entry for entry in primary_entries if entry["metal_stratum"] == stratum]
        if not subset:
            continue
        strata[stratum] = {
            "n_targets": int(len({entry["target"] for entry in subset})),
            "targets": sorted({entry["target"] for entry in subset}),
            "mean_auc": float(np.nanmean([entry["auc"] for entry in subset])),
            "median_auc": float(np.nanmedian([entry["auc"] for entry in subset])),
            "mean_ef5": float(np.nanmean([entry["ef5_percent"] for entry in subset])),
            "per_target": [
                {
                    "target": entry["target"],
                    "auc": entry["auc"],
                    "auc_ci_cluster": [entry["auc_ci_lo"], entry["auc_ci_hi"]],
                    "ef5_percent": entry["ef5_percent"],
                    "ef5_ci_cluster": [entry["ef5_ci_lo"], entry["ef5_ci_hi"]],
                }
                for entry in subset
            ],
        }
    result["metal_vs_nonmetal"] = strata
    result["metal_vs_nonmetal"]["primary_caliber"] = "includes_boron_surrogate"
    result["aggregation_scope"] = {
        "primary_targets_only": [entry["target"] for entry in primary_entries],
        "auxiliary_excluded": [entry["target"] for entry in auxiliary_entries],
        "caveat": (
            "heterogeneous; per-target reporting is primary; stratum aggregates are "
            "descriptive only and must not be used as evidence of a uniform metal effect. "
            "Candidate confounders: metal coordination, chemotype concentration (not supported "
            "for the LpxC Ec/Pa difference), activity strength, assay heterogeneity."
        ),
    }
    result["primary_per_target"] = [
        {
            "target": entry["target"],
            "metal_dependent": entry["metal_dependent"],
            "metal_stratum": entry["metal_stratum"],
            "panel_role": entry.get("panel_role"),
            "n": entry["n"],
            "n_active": entry["n_active"],
            "auc": entry["auc"],
            "auc_ci_cluster": [entry["auc_ci_lo"], entry["auc_ci_hi"]],
            "ef5_percent": entry["ef5_percent"],
            "ef5_ci_cluster": [entry["ef5_ci_lo"], entry["ef5_ci_hi"]],
            "auc_property_stratified": entry.get("auc_property_stratified"),
            "auc_property_matched": entry.get("auc_property_matched"),
        }
        for entry in primary_entries
    ]
    if auxiliary_entries:
        result["auxiliary_per_target"] = [
            {
                "target": entry["target"],
                "panel_role": entry.get("panel_role"),
                "n": entry["n"],
                "n_active": entry["n_active"],
                "auc": entry["auc"],
                "ef5_percent": entry["ef5_percent"],
                "note": "low-confidence/auxiliary target; excluded from primary aggregates, SI only",
            }
            for entry in auxiliary_entries
        ]
    dhfr_entries = [entry for entry in entries if entry["target"] == "DHFR_Ec"]
    if dhfr_entries:
        entry = dhfr_entries[0]
        result["dhfr_summary"] = {
            "target": "DHFR_Ec",
            "metal_stratum": entry["metal_stratum"],
            "n": entry["n"],
            "n_active": entry["n_active"],
            "auc": entry["auc"],
            "auc_ci_cluster": [entry["auc_ci_lo"], entry["auc_ci_hi"]],
            "ef5_percent": entry["ef5_percent"],
            "ef5_ci_cluster": [entry["ef5_ci_lo"], entry["ef5_ci_hi"]],
            "bedroc20": entry["bedroc20"],
        }
    # internal provisional report: panel incomplete but >=5 targets incl DHFR
    if not panel_complete and allow_partial_provisional:
        nonmetal_done = any(entry["metal_stratum"] == "non-metal" for entry in entries)
        provisional_eligible = bool(
            "DHFR_Ec" in done_targets and nonmetal_done and len(done_targets) >= 5
        )
        result["provisional_eligible"] = provisional_eligible
        if provisional_eligible:
            result["provisional"] = True
            result["panel_complete"] = False
            result["missing_targets"] = list(missing_targets)
            result["provisional_note"] = (
                "INTERNAL provisional report — panel incomplete; not for paper/abstract; "
                "formal full JSON only at 9/9"
            )
            result["dry_run"] = False
            result["forced_dry_run"] = False
        else:
            result["provisional"] = False
            result["provisional_reason"] = (
                "requires >=5 done targets including DHFR_Ec and a non-metal control; "
                f"done={len(done_targets)}, missing={missing_targets}"
            )
    elif not panel_complete:
        result["provisional_eligible"] = False
    manifest_path = None
    if metadata:
        if metadata_path is not None:
            manifest_path = str(metadata_path)
        else:
            existing = [str(p) for p in DEFAULT_METADATA_CANDIDATES if p.exists()]
            manifest_path = existing[0] if existing else None
    result["compute_settings"] = {
        "amendment": (
            "Amendments #7-#10 (FROZEN): main panel = up to 20 matched decoys per active "
            "(deterministic seed 42, decoy_subsample_sha256 recorded; observed 9.82-19.98, "
            "median 15.15) + Vina exhaustiveness 4; sensitivity arms = (a) LpxC_Ec/DHFR "
            "1,000-ligand stratified subsets exh4 vs exh8, (b) decoy-ratio bootstrap "
            "subsampling 10:1/15:1/20:1 (no extra docking); non-EF1 targets exh2+20:1; "
            "box 26 A for LpxC_Ec/GyrB_Ec/Gyrase_Ec/LepB_Ec and 22 A otherwise (#8); "
            "box-clipped ligands excluded (#10); original 50:1 arm deleted; H1-H3 unchanged"
        ),
        "decoy_ratio_values": sorted(
            {entry["decoy_ratio"] for entry in entries if entry.get("decoy_ratio") is not None}, key=str
        ),
        "observed_decoy_ratio_by_entry": {
            f"{entry['target']}|{entry['condition']}|{entry['method']}": entry["decoy_ratio_observed"]
            for entry in entries
        },
        "sampling_seed_values": sorted(
            {entry["sampling_seed"] for entry in entries if entry.get("sampling_seed") is not None}, key=str
        ),
        "exhaustiveness_values": sorted(
            {entry["exhaustiveness"] for entry in entries if entry.get("exhaustiveness") is not None}, key=str
        ),
        "panel_role_values": sorted(
            {entry["panel_role"] for entry in entries if entry.get("panel_role") is not None}, key=str
        ),
        "box_size_A_values": sorted(
            {entry.get("box_size_A") for entry in entries if entry.get("box_size_A") is not None}, key=str
        ),
        "box_clipped_rule": "max_dim_A > box_size_A -> excluded from primary enrichment (Amendment #10)",
        "prep_missing_rule": (
            "status=prep_missing -> excluded from enrichment denominators; "
            "counts reported as n_prep_missing / actives / decoys; warning if actives >5%"
        ),
        "ef_definition": (
            "EF(f) = (hits in top k / n_act) / (k / N), with k = max(1, round(f*N)) "
            "(D-line loader implementation); EF1 additionally gated on n_act>=150"
        ),
        "dry_run": bool(result.get("dry_run")),
        "incomplete": bool(incomplete_targets),
        "metadata_manifest": manifest_path,
        "sampling_meta_targets": sorted(sampling_meta.keys()) if sampling_meta else [],
        "box_clipped_match_mode": box_clipped_match_mode,
        "metadata_source_counts": {
            source: int(sum(1 for entry in entries if entry["metadata_source"] == source))
            for source in {entry["metadata_source"] for entry in entries}
        },
    }
    if write and result.get("provisional") and not result.get("panel_complete"):
        FIGURES_DIR.mkdir(parents=True, exist_ok=True)
        json_path = FIGURES_DIR / "trackA_enrichment_partial_provisional.json"
        csv_path = FIGURES_DIR / "trackA_enrichment_partial_provisional.csv"
        result["provisional"] = True
        result["panel_complete"] = False
        result["missing_targets"] = list(result.get("panel", {}).get("missing_targets", []))
        result["write_mode"] = "partial_provisional"
        json_path.write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
        pd.DataFrame(entries).to_csv(csv_path, index=False)
        result["report_json"] = str(json_path)
        result["report_csv"] = str(csv_path)
    elif write:
        FIGURES_DIR.mkdir(parents=True, exist_ok=True)
        suffix = tag + ("_dryrun" if result.get("dry_run") and not tag.endswith("dryrun") else "")
        json_path = FIGURES_DIR / f"trackA_enrichment_{suffix}.json"
        csv_path = FIGURES_DIR / f"trackA_enrichment_{suffix}.csv"
        result["write_mode"] = (
            "full" if result.get("panel_complete") and not result.get("dry_run")
            else ("dry_run" if result.get("dry_run") else "report")
        )
        json_path.write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
        pd.DataFrame(entries).to_csv(csv_path, index=False)
        result["report_json"] = str(json_path)
        result["report_csv"] = str(csv_path)
    return result


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Track A enrichment framework")
    parser.add_argument("--scores-dir", default=str(DEFAULT_TRACKA_DIR / "pilot_v2"))
    parser.add_argument("--pattern", default="*_scores.csv")
    parser.add_argument("--tag", default="pilot_v2")
    parser.add_argument("--n-boot", type=int, default=2000)
    parser.add_argument("--metadata", default=None,
                        help="per-target metadata CSV/JSON (decoy_ratio / sampling_seed / sha256 / exhaustiveness / panel_role / box_size_A)")
    parser.add_argument("--box-clipped", default=None,
                        help="A-line box_clipped_ligands.csv (Amendment #10; default: auto-detect)")
    parser.add_argument("--allow-incomplete", action="store_true",
                        help="allow a dry-run on an in-progress (partial) score table; never cite those numbers")
    parser.add_argument("--allow-partial-provisional", action="store_true",
                        help="write an internal provisional JSON when panel incomplete but >=5 targets incl DHFR")
    args = parser.parse_args()
    out = run_tracka_framework(
        args.scores_dir,
        args.pattern,
        tag=args.tag,
        n_boot=args.n_boot,
        metadata_path=args.metadata,
        box_clipped_path=args.box_clipped,
        allow_incomplete=bool(args.allow_incomplete),
        dry_run=bool(args.allow_incomplete),
        allow_partial_provisional=bool(args.allow_partial_provisional),
    )
    print(json.dumps(
        {
            "tag": out["tag"],
            "n_files": out["n_files"],
            "entries": [
                {k: entry[k] for k in ("target", "method", "condition", "n", "n_active", "auc", "ef5_percent", "bedroc20", "metal_stratum")}
                for entry in out["entries"]
            ],
            "report_json": out.get("report_json"),
            "report_csv": out.get("report_csv"),
        },
        indent=2,
        default=str,
    ))
