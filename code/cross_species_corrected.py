import json
from pathlib import Path

import numpy as np
import pandas as pd

from maes.evaluation import auroc
from maes.tracka import load_score_table

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/trackA/axes"
TRACKA = ROOT / "results/trackA/full"


def accepted(target):
    frame = load_score_table(TRACKA / target / f"{target}_scores.csv")
    clipped = frame["box_clipped"].map(lambda v: str(v).strip().lower() in {"true", "1", "1.0", "yes"})
    frame = frame[~clipped]
    return frame[np.isfinite(frame["score"].astype(float))]


ec = accepted("LpxC_Ec")
pa = accepted("LpxC_Pa")
ec_dec = ec[ec["label"] == 0].copy().set_index("ligand_id_norm")
pa_dec = pa[pa["label"] == 0].copy().set_index("ligand_id_norm")
shared = sorted(set(ec_dec.index.astype(str)) & set(pa_dec.index.astype(str)))
ec_dec = ec_dec.loc[shared]
pa_dec = pa_dec.loc[shared]
ec_act = ec[ec["label"] == 1].copy()
pa_act = pa[pa["label"] == 1].copy()


def auc_pair(ec_act_scores, pa_act_scores, dec_idx):
    dec_ec = ec_dec["score"].to_numpy()[dec_idx]
    dec_pa = pa_dec["score"].to_numpy()[dec_idx]
    y_ec = np.r_[np.ones(len(ec_act_scores)), np.zeros(len(dec_ec))]
    y_pa = np.r_[np.ones(len(pa_act_scores)), np.zeros(len(dec_pa))]
    return (float(auroc(y_ec, np.r_[ec_act_scores, dec_ec])),
            float(auroc(y_pa, np.r_[pa_act_scores, dec_pa])))


ec_act_scores = ec_act["score"].to_numpy()
pa_act_scores = pa_act["score"].to_numpy()
point_ec, point_pa = auc_pair(ec_act_scores, pa_act_scores, np.arange(len(shared)))
rng = np.random.default_rng(42)

# Correct decoy-level bootstrap: resample shared decoys with replacement; keep actives fixed.
dec_draw_ec, dec_draw_pa, dec_delta = [], [], []
for _ in range(10000):
    dec_idx = rng.integers(0, len(shared), len(shared))
    a, b = auc_pair(ec_act_scores, pa_act_scores, dec_idx)
    dec_draw_ec.append(a); dec_draw_pa.append(b); dec_delta.append(b - a)

# Active-level bootstrap: resample each target's actives with replacement; keep decoys fixed.
act_draw_ec, act_draw_pa, act_delta = [], [], []
dec_idx_full = np.arange(len(shared))
for _ in range(10000):
    ec_idx = rng.integers(0, len(ec_act_scores), len(ec_act_scores))
    pa_idx = rng.integers(0, len(pa_act_scores), len(pa_act_scores))
    a, b = auc_pair(ec_act_scores[ec_idx], pa_act_scores[pa_idx], dec_idx_full)
    act_draw_ec.append(a); act_draw_pa.append(b); act_delta.append(b - a)


def ci(values):
    return [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]


out = {
    "point": {"auc_ec": point_ec, "auc_pa": point_pa, "delta_auc_pa_minus_ec": point_pa - point_ec,
              "n_shared_decoys": len(shared), "n_ec_actives": len(ec_act), "n_pa_actives": len(pa_act)},
    "decoy_level_bootstrap": {
        "method": "shared decoys resampled with replacement; actives fixed; 10,000 resamples; seed 42",
        "auc_ec_ci95": ci(dec_draw_ec), "auc_pa_ci95": ci(dec_draw_pa),
        "delta_auc_ci95": ci(dec_delta),
    },
    "active_level_bootstrap": {
        "method": "actives resampled with replacement independently per target; decoys fixed; 10,000 resamples; seed 42",
        "auc_ec_ci95": ci(act_draw_ec), "auc_pa_ci95": ci(act_draw_pa),
        "delta_auc_ci95": ci(act_delta),
    },
    "old_cross_species_csv_ci": [0.513489, 0.543268],
    "old_ci_issue": "The old cross_species.csv bootstrap used isin() on resampled IDs, collapsing with-replacement multiplicity to a no-replacement subset; this corrected version keeps multiplicity.",
}
(OUT / "cross_species_bootstrap_corrected.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print(json.dumps(out, indent=2))
