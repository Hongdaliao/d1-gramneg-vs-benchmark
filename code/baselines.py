"""Naive 2-D baselines for the D1 benchmark (pure CPU, no GPU, no internet).

The paper needs at least two cheap, *non-structure* baselines so that a
reviewer cannot say the comparison only involves expensive methods
(reporting checklist B3):

* ``ECFP4OneNNBaseline``   — max ECFP4 Tanimoto similarity to the known actives
                             (1-nearest-neighbour, ligand-side only)
* ``ECFP4RandomForestBaseline`` — Random Forest on ECFP4 fingerprints

Both expose a uniform interface that can be written into the same long-format
score matrix as Vina / DrugCLIP / Boltz-2::

    target_id, compound_id, method, score

so they can be joined directly onto ``results/trackA_scores.csv``.

Usage::

    base = ECFP4OneNNBaseline().fit(actives["smiles"])
    scores = base.predict(candidates["smiles"])

    rf = ECFP4RandomForestBaseline(random_state=42).fit(actives["smiles"], y)
    prob = rf.predict_proba(candidates["smiles"])[:, 1]
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

__all__ = [
    "SCORE_MATRIX_COLUMNS",
    "ECFP4RandomForestBaseline",
    "ECFP4OneNNBaseline",
    "build_score_matrix",
    "make_ecfp4",
    "tanimoto_matrix",
]

#: Long-format scoring matrix shared by every method (Vina / DrugCLIP / baselines).
SCORE_MATRIX_COLUMNS = ("target_id", "compound_id", "method", "score")

_DEFAULT_RADIUS = 2          # ECFP4 == radius 2
_DEFAULT_N_BITS = 2048


def _rdkit_available() -> bool:
    try:  # pragma: no cover - environment dependent
        import rdkit  # noqa: F401

        return True
    except Exception:
        return False


def make_ecfp4(
    smiles: Sequence[str],
    radius: int = _DEFAULT_RADIUS,
    n_bits: int = _DEFAULT_N_BITS,
    use_features: bool = False,
) -> List[Any]:
    """Morgan/ECFP4 fingerprints; invalid SMILES become ``None``."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    fps: List[Any] = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(str(smi)) if smi is not None else None
        if mol is None:
            fps.append(None)
            continue
        if use_features:
            fps.append(AllChem.GetFeatureFingerprintAsBitVect(mol, int(radius), nBits=int(n_bits)))
        else:
            fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, int(radius), nBits=int(n_bits)))
    return fps


def tanimoto_matrix(
    query_fps: Sequence[Any], reference_fps: Sequence[Any]
) -> np.ndarray:
    """Dense Tanimoto similarity matrix (invalid fingerprints -> 0 rows)."""
    from rdkit import DataStructs

    n_q, n_r = len(query_fps), len(reference_fps)
    out = np.zeros((n_q, n_r), dtype=np.float32)
    if n_q == 0 or n_r == 0:
        return out
    for i, fp_q in enumerate(query_fps):
        if fp_q is None:
            continue
        valid_refs = [(j, r) for j, r in enumerate(reference_fps) if r is not None]
        if not valid_refs:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(fp_q, [r for _, r in valid_refs])
        for (j, _), sim in zip(valid_refs, sims):
            out[i, j] = float(sim)
    return out


@dataclass
class ECFP4OneNNBaseline:
    """1-NN on ECFP4: score = max Tanimoto similarity to the training actives."""

    radius: int = _DEFAULT_RADIUS
    n_bits: int = _DEFAULT_N_BITS
    actives: List[Any] = field(default_factory=list)
    name: str = "ecfp4_1nn"
    fitted: bool = False

    def fit(self, actives: Sequence[str]) -> "ECFP4OneNNBaseline":
        self.actives = make_ecfp4(actives, radius=self.radius, n_bits=self.n_bits)
        if not any(fp is not None for fp in self.actives):
            raise ValueError("no valid active SMILES supplied to ECFP4OneNNBaseline.fit")
        self.fitted = True
        return self

    def predict(self, smiles: Sequence[str]) -> np.ndarray:
        if not self.fitted:
            raise RuntimeError("ECFP4OneNNBaseline is not fitted")
        query = make_ecfp4(smiles, radius=self.radius, n_bits=self.n_bits)
        sims = tanimoto_matrix(query, self.actives)
        return np.clip(sims.max(axis=1) if sims.size else np.zeros(len(smiles)), 0.0, 1.0)


@dataclass
class ECFP4RandomForestBaseline:
    """2-D classifier: Random Forest on ECFP4 bits (probability of being active)."""

    n_estimators: int = 500
    max_depth: Optional[int] = None
    min_samples_leaf: int = 1
    random_state: int = 42
    # keep n_jobs=1 by default: -1 uses multiprocessing, which some sandboxes
    # (and the D1 offline environment) block at the named-pipe layer
    n_jobs: int = 1
    radius: int = _DEFAULT_RADIUS
    n_bits: int = _DEFAULT_N_BITS
    name: str = "rf_ecfp4"
    _model: Any = None
    _base_rate: float = 0.5
    fitted: bool = False

    def _matrix(self, smiles: Sequence[str]) -> np.ndarray:
        from rdkit import DataStructs

        fps = make_ecfp4(smiles, radius=self.radius, n_bits=self.n_bits)
        n_bits = int(self.n_bits)
        out = np.zeros((len(fps), n_bits), dtype=np.uint8)
        for i, fp in enumerate(fps):
            if fp is None:
                continue
            DataStructs.ConvertToNumpyArray(fp, out[i])
        return out

    def fit(self, smiles: Sequence[str], y: Sequence[int]) -> "ECFP4RandomForestBaseline":
        from sklearn.ensemble import RandomForestClassifier

        y_arr = np.asarray(y, dtype=int)
        if y_arr.size != len(smiles):
            raise ValueError("smiles and y differ in length")
        if len(set(y_arr.tolist())) < 2:
            raise ValueError("ECFP4RandomForestBaseline needs both classes")
        x = self._matrix(smiles)
        self._base_rate = float(y_arr.mean())
        self._model = RandomForestClassifier(
            n_estimators=int(self.n_estimators),
            max_depth=self.max_depth,
            min_samples_leaf=int(self.min_samples_leaf),
            random_state=int(self.random_state),
            n_jobs=int(self.n_jobs),
        ).fit(x, y_arr)
        self.fitted = True
        return self

    def predict_proba(self, smiles: Sequence[str]) -> np.ndarray:
        if not self.fitted:
            raise RuntimeError("ECFP4RandomForestBaseline is not fitted")
        x = self._matrix(smiles)
        proba = np.asarray(self._model.predict_proba(x), dtype=float)
        if proba.shape[1] == 1:  # degenerate fold
            return np.full((x.shape[0], 2), self._base_rate)
        return proba

    def predict(self, smiles: Sequence[str]) -> np.ndarray:
        return self.predict_proba(smiles)[:, 1]


def align_baseline_scores(
    scores: pd.DataFrame,
    poses: pd.DataFrame,
    subset: str = "primary",
    target_col: str = "target",
    compound_col: str = "ligand_id",
) -> pd.DataFrame:
    """Restrict baseline scores to the W2 analysis set (inhibitor subset by default).

    ``scores`` is the long-format matrix produced by :func:`build_score_matrix`
    (``target_id, compound_id, method, score``); ``poses`` is the Track B
    manifest.  The result carries a ``subset`` column so the main (n=40
    inhibitor poses) and sensitivity (n=68) analyses cannot be mixed by accident.
    """
    from .api import attach_pose_ids, select_analysis_set

    frame = attach_pose_ids(poses)
    frame = select_analysis_set(frame, subset)
    keep = {
        (str(target), str(compound))
        for target, compound in zip(frame[target_col], frame[compound_col])
    }
    missing = [c for c in ("target_id", "compound_id") if c not in scores.columns]
    if missing:
        raise KeyError(f"baseline score matrix lacks {missing}")
    pairs = list(zip(scores["target_id"].astype(str), scores["compound_id"].astype(str)))
    aligned = scores.loc[[pair in keep for pair in pairs]].copy()
    aligned["subset"] = subset
    return aligned


def build_score_matrix(
    target_id: str,
    candidate_ids: Sequence[str],
    smiles: Sequence[str],
    actives: Sequence[str] = (),
    decoys: Sequence[str] = (),
    baselines: Optional[Sequence[str]] = None,
    seed: int = 42,
) -> pd.DataFrame:
    """Run the requested baselines on one target and return a long score matrix.

    Parameters
    ----------
    target_id : str
    candidate_ids, smiles : candidate compounds (same length)
    actives : known actives — the 1-NN reference set and the positives of the
        2-D classifier.
    decoys : known inactives used to train the 2-D Random Forest (the RF needs
        both classes; in the D1 pipeline these are the real ChEMBL inactives).
    baselines : subset of ``("ecfp4_1nn", "rf_ecfp4")``; defaults to both.
    """
    baselines = tuple(baselines or ("ecfp4_1nn", "rf_ecfp4"))
    rows: List[Dict[str, Any]] = []
    if "ecfp4_1nn" in baselines:
        model = ECFP4OneNNBaseline().fit(actives)
        scores = model.predict(smiles)
        rows.extend(
            {"target_id": target_id, "compound_id": cid, "method": model.name, "score": float(s)}
            for cid, s in zip(candidate_ids, scores)
        )
    if "rf_ecfp4" in baselines:
        if not decoys:
            raise ValueError("rf_ecfp4 needs decoys (known inactives) to train the 2-D classifier")
        model = ECFP4RandomForestBaseline(random_state=seed).fit(
            list(actives) + list(decoys), [1] * len(actives) + [0] * len(decoys)
        )
        scores = model.predict(smiles)
        rows.extend(
            {"target_id": target_id, "compound_id": cid, "method": model.name, "score": float(s)}
            for cid, s in zip(candidate_ids, scores)
        )
    frame = pd.DataFrame(rows, columns=list(SCORE_MATRIX_COLUMNS))
    return frame
