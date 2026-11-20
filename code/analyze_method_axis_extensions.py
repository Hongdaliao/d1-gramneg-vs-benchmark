"""Amendment #12/#13/#14 extension: DrugCLIP and Uni-Dock on the 806 subset.

Run only after both aggregate score tables are ready (>=780 ok each).  Does not
touch v3; writes S29-S31, six-method pooled JSON and Fig 7 A-F draft.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(r"D:\harness\projects\d1_gramneg_vs")
BASE = ROOT / "results/trackA/method_axis"
OUT = BASE / "analysis"
SUBSET = BASE / "subsets_diffdock.csv"
DRUGCLIP = BASE / "drugclip/drugclip_scores_6_folds.csv"
UNIDOCK = BASE / "unidock/scores/unidock_scores_806.csv"
TARGETS = ["LpxC_Ec", "LpxC_Pa", "DHFR_Ec", "Gyrase_Ec", "GyrA_Ec", "GyrB_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]
EVAL = [t for t in TARGETS if t != "LepB_Ec"]
N_BOOT = 2000
SEED = 42


def auc(y, s):
    y = np.asarray(y, int); s = np.asarray(s, float)
    if y.sum() == 0 or y.sum() == len(y):
        return float("nan")
    try:
        return float(roc_auc_score(y, s))
    except ValueError:
        return float("nan")


def boot_ci(y, s, seed):
    y = np.asarray(y, int); s = np.asarray(s, float); pos = np.where(y == 1)[0]; neg = np.where(y == 0)[0]
    rng = np.random.default_rng(seed); vals = []
    for _ in range(N_BOOT):
        idx = np.r_[rng.choice(pos, len(pos), replace=True), rng.choice(neg, len(neg), replace=True)]
        v = auc(y[idx], s[idx])
        if np.isfinite(v): vals.append(v)
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))) if vals else (float("nan"), float("nan"))


def top5(y, s):
    y = np.asarray(y, int); s = np.asarray(s, float); k = min(5, len(y))
    return float(y[np.argsort(-s)[:k]].sum() / min(5, max(1, y.sum())))


def rho(x, y):
    x = pd.Series(np.asarray(x, float)); y = pd.Series(np.asarray(y, float)); ok = x.notna() & y.notna()
    return float(x[ok].rank().corr(y[ok].rank())) if ok.sum() >= 3 else float("nan")


def load_scores(path: Path, method: str) -> pd.DataFrame:
    d = pd.read_csv(path)
    if method == "DrugCLIP":
        d = d.rename(columns={"score": "method_score"})
        d["method_score"] = pd.to_numeric(d["method_score"], errors="coerce")  # higher better
    elif method == "Uni-Dock":
        d = d.rename(columns={"affinity": "method_score"})
        d["method_score"] = -pd.to_numeric(d["method_score"], errors="coerce")  # lower affinity better
    else:
        raise ValueError(method)
    d["method"] = method
    return d[["target", "ligand_id", "label", "method_score", "status", "method"]]


def coverage(d: pd.DataFrame, method: str) -> pd.DataFrame:
    c = d.copy(); c["failed_at"] = ""; bad = c["status"].astype(str).str.lower() != "ok"
    c.loc[bad, "failed_at"] = method.lower() + "_failure"
    c.to_csv(OUT / f"coverage_{method.lower().replace('-', '')}.csv", index=False)
    return c


def main() -> Dict[str, Any]:
    for p in (DRUGCLIP, UNIDOCK):
        if not p.exists():
            raise SystemExit(f"missing {p}")
    sub = pd.read_csv(SUBSET)
    vina = sub[["target", "subset_id", "ligand_id", "label", "vina_score"]].copy()
    vina["vina_score"] = pd.to_numeric(vina["vina_score"], errors="coerce")

    all_metrics, all_paired, summaries = [], [], {}
    for method, path in (("DrugCLIP", DRUGCLIP), ("Uni-Dock", UNIDOCK)):
        d = load_scores(path, method)
        m = vina.merge(d[["target", "ligand_id", "method_score", "status"]], on=["target", "ligand_id"], how="left")
        cov = coverage(m[["target", "ligand_id", "label", "method_score", "status"]], method)
        rows = []
        for t in TARGETS:
            s = m[(m["target"] == t) & (m["status"].astype(str).str.lower() == "ok") & m["method_score"].notna() & m["vina_score"].notna()].copy()
            if s.empty:
                continue
            y = s["label"].to_numpy(int); sc = s["method_score"].to_numpy(float)
            lo, hi = boot_ci(y, sc, SEED + TARGETS.index(t) + (0 if method == "DrugCLIP" else 100))
            if t == "LepB_Ec":
                rows.append({"target": t, "method": method, "n_ligands": len(m[m.target == t]), "n_analysis": len(s),
                             "n_active": int(y.sum()), "n_decoy": int((y == 0).sum()),
                             "n_na_excluded": int(len(m[m.target == t]) - len(s)), "auc": float("nan"),
                             "auc_ci_lo": float("nan"), "auc_ci_hi": float("nan"), "top5_recall": float("nan"),
                             "rho_confidence_vina": rho(sc, s["vina_score"]), "status": "descriptive only (n=6)"})
                continue
            rows.append({"target": t, "method": method, "n_ligands": len(m[m.target == t]), "n_analysis": len(s),
                         "n_active": int(y.sum()), "n_decoy": int((y == 0).sum()), "n_na_excluded": int(len(m[m.target == t]) - len(s)),
                         "auc": auc(y, sc), "auc_ci_lo": lo, "auc_ci_hi": hi, "top5_recall": top5(y, sc),
                         "rho_confidence_vina": rho(sc, s["vina_score"]), "status": "ok"})
            pairs = [g for _, g in s.groupby("subset_id") if len(g) == 2 and set(g["label"].astype(int)) == {0, 1}]
            if pairs:
                rng = np.random.default_rng(SEED + 400 + TARGETS.index(t)); deltas = []
                for _ in range(N_BOOT):
                    sample = pd.concat([pairs[i] for i in rng.integers(0, len(pairs), len(pairs))], ignore_index=True)
                    yy = sample["label"].to_numpy(int)
                    a = auc(yy, sample["method_score"].to_numpy(float))
                    b = auc(yy, -sample["vina_score"].to_numpy(float))
                    if np.isfinite(a) and np.isfinite(b):
                        deltas.append(a - b)
                all_paired.append({"method": method, "target": t, "n_pairs": len(pairs),
                                   "delta_auc_method_minus_vina": auc(y, sc) - auc(y, -s["vina_score"].to_numpy(float)),
                                   "delta_ci_lo": float(np.percentile(deltas, 2.5)) if deltas else float("nan"),
                                   "delta_ci_hi": float(np.percentile(deltas, 97.5)) if deltas else float("nan"),
                                   "n_boot": len(deltas)})
        frame = pd.DataFrame(rows)
        all_metrics.append(frame)
        targ = frame[(frame["method"] == method) & frame["target"].isin(EVAL) & frame["auc"].notna()]
        if len(targ):
            arr = targ["auc"].to_numpy(float); rng = np.random.default_rng(SEED + (10 if method == "DrugCLIP" else 20))
            draws = [float(np.mean(arr[rng.integers(0, len(arr), len(arr))])) for _ in range(N_BOOT)]
            summaries[method] = {"mean_auc": float(np.mean(arr)), "ci95": [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))],
                                 "n_ok": int((frame["status"] == "ok").sum()), "n_rows": int(len(m)), "targets": targ.set_index("target")["auc"].to_dict()}
        else:
            summaries[method] = {"status": "pending"}
    metrics = pd.concat(all_metrics, ignore_index=True) if all_metrics else pd.DataFrame()
    metrics.to_csv(OUT / "metrics_drugclip_unidock.csv", index=False)
    metrics[metrics["method"] == "DrugCLIP"].to_csv(OUT / "table_S29_drugclip_per_target.csv", index=False)
    metrics[metrics["method"] == "Uni-Dock"].to_csv(OUT / "table_S30_unidock_per_target.csv", index=False)
    paired = pd.DataFrame(all_paired)
    paired.to_csv(OUT / "table_S31b_drugclip_unidock_paired.csv", index=False)
    (OUT / "metrics_pooled_suite.json").write_text(json.dumps({"methods": summaries, "paired": all_paired}, indent=2), encoding="utf-8")
    # six-method pooled table (include existing DD/Vina/Boltz + new methods; ECFP4 is full-panel reference not subset)
    existing = json.loads((OUT / "metrics_pooled.json").read_text(encoding="utf-8"))["pooled"]
    six = {
        "Vina": existing.get("vina_auc_mean"),
        "DiffDock": existing.get("diffdock_auc_mean"),
        "Boltz-2": existing.get("methods", {}).get("Boltz-2", {}).get("mean_auc"),
        "DrugCLIP": summaries.get("DrugCLIP", {}).get("mean_auc"),
        "Uni-Dock": summaries.get("Uni-Dock", {}).get("mean_auc"),
        "ECFP4": "full-panel test-only reference only (median 0.775; subset reference withdrawn, Addendum 12a)",
    }
    six_rows = [{"method": k, "pooled_mean_auc": v, "note": "full-panel 2D reference; no subset comparison" if k == "ECFP4" else "subset-balanced; not comparable to Table 1"} for k, v in six.items()]
    pd.DataFrame(six_rows).to_csv(OUT / "table_S31_six_method_pooled.csv", index=False)
    (OUT / "six_method_pooled.json").write_text(json.dumps({"pooled_mean_auc": six, "coverage": {k: {"n_ok": summaries.get(k, {}).get("n_ok"), "n_rows": summaries.get(k, {}).get("n_rows")} for k in ("DrugCLIP", "Uni-Dock")}}, indent=2), encoding="utf-8")
    return {"metrics": str(OUT / "metrics_drugclip_unidock.csv"), "summaries": summaries}


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, default=str))
