"""Ligand-property sensitivity helpers for Track A (D line).

Implements the property-stratified / nearest-neighbour sensitivity requested by
Lead (2026-11-18): bin actives+decoys by MW / cLogP / HBD, compute within-bin
AUROC/EF5 and weight by bin size; also match each active to its nearest decoy in
standardised (MW, cLogP, HBD) space and recompute AUROC.

The functions are deliberately dependency-light at import time (RDKit is imported
lazily) so the rest of ``maes.tracka`` can run without generating descriptors.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .evaluation import auroc, enrichment_factor

PROPERTY_COLUMNS = ("mw", "logp", "hbd")


def compute_rdkit_properties(smiles: Sequence[str]) -> pd.DataFrame:
    """MW, cLogP and HBD count for each SMILES (NaN when RDKit fails)."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors

    RDLogger.DisableLog("rdApp.*")
    rows: List[Dict[str, Any]] = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(str(smi)) if isinstance(smi, str) else None
        if mol is None:
            rows.append({"mw": np.nan, "logp": np.nan, "hbd": np.nan})
            continue
        rows.append(
            {
                "mw": float(Descriptors.MolWt(mol)),
                "logp": float(Descriptors.MolLogP(mol)),
                "hbd": int(Descriptors.NumHDonors(mol)),
            }
        )
    return pd.DataFrame(rows)


def load_ligand_properties(manifest_path: str) -> pd.DataFrame:
    """Load ``ligand_id -> (mw, logp, hbd)`` from a B-line manifest."""
    manifest = pd.read_csv(manifest_path)
    if "ligand_id" not in manifest.columns or "smiles" not in manifest.columns:
        raise KeyError(f"{manifest_path} needs ligand_id / smiles columns")
    props = compute_rdkit_properties(manifest["smiles"].tolist())
    props.insert(0, "ligand_id", manifest["ligand_id"].astype(str).values)
    return props


def _bin_sizes(labels: np.ndarray, n_bins: int, seed: int = 42) -> np.ndarray:
    order = np.argsort(labels, kind="mergesort")
    bins = np.empty(len(labels), dtype=int)
    for rank, idx in enumerate(order):
        bins[idx] = min(n_bins - 1, int(rank * n_bins / len(labels)))
    return bins


def _weighted_bin_auc(values: np.ndarray, labels: np.ndarray, scores: np.ndarray, n_bins: int) -> Dict[str, Any]:
    finite = np.isfinite(values)
    if finite.sum() == 0:
        return {"auc": float("nan"), "coverage": 0.0, "n_bins_used": 0}
    v, y, s = values[finite], labels[finite], scores[finite]
    bins = _bin_sizes(v, n_bins)
    aucs, weights = [], []
    for b in range(n_bins):
        mask = bins == b
        if mask.sum() == 0:
            continue
        yb, sb = y[mask], s[mask]
        if yb.sum() == 0 or yb.sum() == len(yb):
            continue
        aucs.append(float(auroc(yb, sb)))
        weights.append(float(mask.sum()))
    if not aucs:
        return {"auc": float("nan"), "coverage": 0.0, "n_bins_used": 0}
    coverage = float(np.sum(weights) / len(y))
    return {
        "auc": float(np.average(aucs, weights=weights)),
        "coverage": coverage,
        "n_bins_used": int(len(aucs)),
    }


def _nearest_decoy_auc(
    labels: np.ndarray,
    scores: np.ndarray,
    props: np.ndarray,
    k: int = 1,
) -> Dict[str, Any]:
    from sklearn.neighbors import NearestNeighbors

    finite = np.isfinite(scores) & np.all(np.isfinite(props), axis=1)
    y, s, p = labels[finite], scores[finite], props[finite]
    active_idx = np.where(y == 1)[0]
    decoy_idx = np.where(y == 0)[0]
    if len(active_idx) == 0 or len(decoy_idx) < 2:
        return {"auc": float("nan"), "n_matched_decoys": 0}
    mean = p.mean(axis=0)
    std = p.std(axis=0)
    std[std == 0] = 1.0
    z = (p - mean) / std
    nn = NearestNeighbors(n_neighbors=min(int(k), len(decoy_idx)))
    nn.fit(z[decoy_idx])
    _, neighbour_idx = nn.kneighbors(z[active_idx])
    matched = decoy_idx[np.unique(neighbour_idx.ravel())]
    if len(matched) < 2:
        return {"auc": float("nan"), "n_matched_decoys": int(len(matched))}
    y_m = np.concatenate([np.ones(len(active_idx)), np.zeros(len(matched))])
    s_m = np.concatenate([s[active_idx], s[matched]])
    return {
        "auc": float(auroc(y_m, s_m)),
        "n_matched_decoys": int(len(matched)),
        "k": int(k),
    }


def property_sensitivity(
    labels: Sequence[int],
    scores: Sequence[float],
    properties: pd.DataFrame,
    raw_auc: float,
    n_bins: int = 3,
) -> Dict[str, Any]:
    """Property-stratified and nearest-neighbour sensitivity for one entry.

    Returns per-property weighted AUROC, the mean across properties, the
    nearest-decoy-matched AUROC, and deltas versus the raw AUROC.
    """
    y = np.asarray(labels, dtype=int)
    s = np.asarray(scores, dtype=float)
    out: Dict[str, Any] = {
        "method": "MW/logP/HBD quantile bins (3) weighted by bin size; nearest-decoy match in standardised property space",
        "n_bins": int(n_bins),
        "raw_auc": float(raw_auc),
        "by_property": {},
        "delta_vs_raw": {},
    }
    by_property: Dict[str, float] = {}
    for column in PROPERTY_COLUMNS:
        if column not in properties.columns:
            continue
        values = pd.to_numeric(properties[column], errors="coerce").to_numpy(dtype=float)
        result = _weighted_bin_auc(values, y, s, n_bins)
        out["by_property"][column] = result
        if np.isfinite(result["auc"]):
            by_property[column] = result["auc"]
            out["delta_vs_raw"][column] = float(result["auc"] - raw_auc)
    out["auc_property_stratified"] = (
        float(np.mean(list(by_property.values()))) if by_property else float("nan")
    )
    out["delta_vs_raw"]["stratified_mean"] = (
        float(out["auc_property_stratified"] - raw_auc) if np.isfinite(out["auc_property_stratified"]) else float("nan")
    )
    prop_matrix = properties[list(PROPERTY_COLUMNS)].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float) if all(
        column in properties.columns for column in PROPERTY_COLUMNS
    ) else np.full((len(properties), len(PROPERTY_COLUMNS)), np.nan)
    matched = _nearest_decoy_auc(y, s, prop_matrix)
    out["auc_property_matched"] = matched["auc"]
    out["n_matched_decoys"] = matched["n_matched_decoys"]
    out["delta_vs_raw"]["matched"] = (
        float(matched["auc"] - raw_auc) if np.isfinite(matched["auc"]) else float("nan")
    )
    return out


def largest_scaffold_group(smiles: Sequence[str]) -> Dict[str, Any]:
    """Largest Bemis-Murcko scaffold group among the supplied SMILES."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold

    RDLogger.DisableLog("rdApp.*")
    scaffolds: List[str] = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(str(smi)) if isinstance(smi, str) else None
        if mol is None:
            scaffolds.append("")
            continue
        try:
            scaffolds.append(MurckoScaffold.MurckoScaffoldSmiles(mol=mol))
        except Exception:
            scaffolds.append("")
    series = pd.Series(scaffolds)
    counts = series.value_counts()
    if len(counts) == 0:
        return {"indices": [], "scaffold": None, "count": 0}
    top = str(counts.index[0])
    return {"indices": list(np.where(series == top)[0]), "scaffold": top, "count": int(counts.iloc[0])}


def largest_ecfp_cluster(smiles: Sequence[str], cutoff: float = 0.4) -> Dict[str, Any]:
    """Largest ECFP4 Butina cluster (Tanimoto distance cutoff, default 0.4)."""
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdMolDescriptors
    from rdkit.ML.Cluster import Butina

    RDLogger.DisableLog("rdApp.*")
    fps = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(str(smi)) if isinstance(smi, str) else None
        fps.append(
            rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            if mol is not None else None
        )
    valid_idx = [i for i, fp in enumerate(fps) if fp is not None]
    if len(valid_idx) < 2:
        return {"indices": [], "count": 0, "n_clusters": 0}
    dists: List[float] = []
    for i in range(1, len(valid_idx)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[valid_idx[i]], [fps[j] for j in valid_idx[:i]])
        dists.extend([1.0 - sim for sim in sims])
    clusters = Butina.ClusterData(dists, len(valid_idx), cutoff, isDistData=True)
    largest = max(clusters, key=len)
    return {
        "indices": [valid_idx[i] for i in largest],
        "count": int(len(largest)),
        "n_clusters": int(len(clusters)),
    }

