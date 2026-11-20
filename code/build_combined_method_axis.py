#!/usr/bin/env python
"""Build combined method-axis long/wide tables for D-line analysis.

Methods (subset 806, Addendum 12a: ECFP4 subset reference withdrawn):
  vina (lower better), unidock (lower), diffdock_main (higher), diffdock_singlechain (higher,
  sensitivity arm only for GyrA/FabI/LepB), boltz2 (higher, ligand_iptm), drugclip (higher).

Outputs (results/trackA/method_axis/):
  combined_method_axis_scores.csv   long: method,target,ligand_id,label,score,status,higher_is_better,arm
  combined_method_axis_wide.csv     wide: one row per (target,ligand_id) + label + per-method scores
  combined_method_axis_coverage.csv per method/target attempted/ok/NA
  combined_method_axis_manifest.json sha256 + preview metrics
"""
import csv
import hashlib
import json
import os

import numpy as np

ROOT = r"D:\harness\projects\d1_gramneg_vs"
B = os.path.join(ROOT, "results", "trackA", "method_axis")
TARGETS = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec",
           "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]
MULTICHAIN = ["GyrA_Ec", "FabI_Ec", "LepB_Ec"]
SUBSETS = os.path.join(B, "subsets_diffdock.csv")


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(p):
    return list(csv.DictReader(open(p, encoding="utf-8-sig"))) if os.path.exists(p) else []


def write_csv(p, rows, fields):
    tmp = p + ".tmp"
    with open(tmp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, p)


def auc(labels, scores):
    lab = np.asarray(labels, dtype=int)
    sc = np.asarray(scores, dtype=float)
    pos = lab == 1
    npos, nneg = int(pos.sum()), int((~pos).sum())
    if npos == 0 or nneg == 0:
        return None
    order = np.argsort(sc, kind="mergesort")
    ranks = np.empty(len(sc))
    ss = sc[order]
    i = 0
    while i < len(sc):
        j = i
        while j + 1 < len(sc) and ss[j + 1] == ss[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    return float((ranks[pos].sum() - npos * (npos + 1) / 2) / (npos * nneg))


def main():
    subs = read_csv(SUBSETS)
    labels = {(r["target"], r["ligand_id"]): int(r["label"]) for r in subs}
    per_target_ligs = {t: [r["ligand_id"] for r in subs if r["target"] == t] for t in TARGETS}

    # --- Vina (lower better) ---
    vina = {}
    for t in TARGETS:
        for r in read_csv(os.path.join(ROOT, "results", "trackA", "full", t, f"{t}_scores.csv")):
            if r.get("status") == "ok" and r.get("score") not in ("", None):
                vina[(t, r["ligand_id"])] = float(r["score"])
    # --- Uni-Dock (lower) ---
    unidock = {}
    for r in read_csv(os.path.join(B, "unidock", "scores", "unidock_scores_806.csv")):
        if r.get("status") == "ok" and r.get("affinity") not in ("", None):
            unidock[(r["target"], r["ligand_id"])] = float(r["affinity"])
    # --- DiffDock main/singlechain (higher) ---
    dd_main, dd_single = {}, {}
    for r in read_csv(os.path.join(B, "diffdock_main_scores.csv")):
        if r.get("status") == "ok" and r.get("confidence_rank1") not in ("", None):
            dd_main[(r["target"], r["ligand_id"])] = float(r["confidence_rank1"])
    for r in read_csv(os.path.join(B, "diffdock_singlechain_scores.csv")):
        if r.get("status") == "ok" and r.get("confidence_rank1") not in ("", None):
            dd_single[(r["target"], r["ligand_id"])] = float(r["confidence_rank1"])
    # --- Boltz-2 (higher) ---
    boltz = {}
    boltz_ligs = {(r["target"], r["ligand_id"]) for r in read_csv(os.path.join(B, "subsets_boltz2.csv"))}
    for r in read_csv(os.path.join(B, "boltz2_scores.csv")):
        if r.get("status") == "ok" and r.get("ligand_iptm") not in ("", None):
            boltz[(r["target"], r["ligand_id"])] = float(r["ligand_iptm"])
    # --- DrugCLIP (higher) ---
    drugclip = {}
    for r in read_csv(os.path.join(B, "drugclip", "drugclip_scores_6_folds.csv")):
        if r.get("status") == "ok" and r.get("score") not in ("", None):
            drugclip[(r["target"], r["ligand_id"])] = float(r["score"])

    methods = [
        ("vina", vina, "lower", TARGETS, "all"),
        ("unidock", unidock, "lower", TARGETS, "all"),
        ("diffdock_main", dd_main, "higher", TARGETS, "all"),
        ("diffdock_singlechain", dd_single, "higher", MULTICHAIN, "all"),
        ("boltz2", boltz, "higher", TARGETS, "nested"),
        ("drugclip", drugclip, "higher", TARGETS, "all"),
    ]
    long_rows, wide = [], {}
    coverage = []
    for method, table, direction, mtargets, universe in methods:
        for t in mtargets:
            n_ok, n_na = 0, 0
            lids_t = ([lid for lid in per_target_ligs[t] if (t, lid) in boltz_ligs]
                      if universe == "nested" else per_target_ligs[t])
            for lid in lids_t:
                key = (t, lid)
                lab = labels[key]
                if key in table:
                    n_ok += 1
                    long_rows.append({"method": method, "target": t, "ligand_id": lid, "label": lab,
                                      "score": f"{table[key]:.6g}", "status": "ok",
                                      "higher_is_better": 1 if direction == "higher" else 0,
                                      "arm": "sensitivity" if method.endswith("singlechain") else "main"})
                    wide.setdefault(key, {})[f"score_{method}"] = table[key]
                else:
                    n_na += 1
                    long_rows.append({"method": method, "target": t, "ligand_id": lid, "label": lab,
                                      "score": "", "status": "NA",
                                      "higher_is_better": 1 if direction == "higher" else 0,
                                      "arm": "sensitivity" if method.endswith("singlechain") else "main"})
            coverage.append({"method": method, "target": t, "attempted": len(lids_t),
                             "ok": n_ok, "na": n_na})
    # long + wide outputs
    write_csv(os.path.join(B, "combined_method_axis_scores.csv"), long_rows,
              ["method", "target", "ligand_id", "label", "score", "status", "higher_is_better", "arm"])
    wide_rows = []
    for (t, lid), labs in labels.items():
        row = {"target": t, "ligand_id": lid, "label": labs}
        for m, _, _, _, _ in methods:
            row[f"score_{m}"] = wide.get((t, lid), {}).get(f"score_{m}", "")
        wide_rows.append(row)
    wide_fields = ["target", "ligand_id", "label"] + [f"score_{m}" for m, _, _, _, _ in methods]
    write_csv(os.path.join(B, "combined_method_axis_wide.csv"), wide_rows, wide_fields)
    write_csv(os.path.join(B, "combined_method_axis_coverage.csv"), coverage,
              ["method", "target", "attempted", "ok", "na"])
    # preview AUC (descriptive; D line owns final stats)
    preview = {}
    for method, table, direction, mtargets, universe in methods:
        pt = {}
        for t in mtargets:
            if t == "LepB_Ec":
                continue
            lids = [lid for lid in per_target_ligs[t]
                    if (t, lid) in table and (universe != "nested" or (t, lid) in boltz_ligs)]
            if len(lids) < 4:
                continue
            ys = [labels[(t, lid)] for lid in lids]
            xs = [table[(t, lid)] for lid in lids]
            a = auc(ys, xs)
            if a is not None:
                pt[t] = round(a, 4)
        preview[method] = {"per_target_auc": pt,
                           "pooled_mean_auc": round(float(np.mean(list(pt.values()))), 4) if pt else None,
                           "n_targets": len(pt)}
    # Lead pooling decisions (2026-11-20):
    #   main pooled for Vina/DiffDock/Boltz/DrugCLIP = 8 targets (exclude LepB)
    #   Uni-Dock pooled = 7 targets (exclude LepB + FabI; FabI coverage-limited)
    #   common 7-target sensitivity subset = exclude FabI for ALL methods
    MAIN8 = [t for t in TARGETS if t != "LepB_Ec"]
    NO_FABI = [t for t in TARGETS if t not in ("FabI_Ec", "LepB_Ec")]
    def pooled_for(method, target_list):
        pt = {}
        for t in target_list:
            lids = [lid for lid in per_target_ligs[t] if (t, lid) in method_table_ref[method]]
            if len(lids) < 4:
                continue
            ys = [labels[(t, lid)] for lid in lids]
            xs = [method_table_ref[method][(t, lid)] for lid in lids]
            a = auc(ys, xs)
            if a is not None:
                pt[t] = round(a, 4)
        return {"per_target_auc": pt,
                "pooled_mean_auc": round(float(np.mean(list(pt.values()))), 4) if pt else None,
                "n_targets": len(pt)}
    method_table_ref = {m: table for m, table, _, _, _ in methods}
    manifest_extra = {
        "pooling_rules_lead_2026_11_20": {
            "main_8targets": ["Vina", "DiffDock", "Boltz-2", "DrugCLIP"],
            "unidock_7targets_excluding_FabI": True,
            "common_7target_sensitivity_subset": NO_FABI,
            "note": "FabI_Ec excluded from Uni-Dock pooled (13/50 actives NA_positive_affinity, selective missingness).",
        },
        "pooled_preview_default": {m: pooled_for(m, MAIN8 if m != "unidock" else NO_FABI) for m in method_table_ref
                                    if m != "diffdock_singlechain"},
        "pooled_preview_common7": {m: pooled_for(m, NO_FABI) for m in method_table_ref if m != "diffdock_singlechain"},
    }
    manifest = {
        "addendum_12a": "ECFP4 subset reference withdrawn (leakage); not included here",
        "preview_only_warning": "preview_auc values are SUBSET-BALANCED descriptive previews for QC; "
                                "they are NOT the pre-registered final statistics (D line computes final "
                                "pooled/target-level bootstrap). Vina subset-balanced AUC differs strongly "
                                "from its full-panel value (decoy-ratio change), so do not compare with Table 1.",
        "unidock_positive_affinity_rule": "status != ok (incl NA_positive_affinity) excluded from scores/AUC",
        "methods": [m for m, _, _, _, _ in methods],
        "score_direction": {m: d for m, _, d, _, _ in methods},
        "n_rows_long": len(long_rows),
        "coverage": {m: {"attempted": sum(c["attempted"] for c in coverage if c["method"] == m),
                         "ok": sum(c["ok"] for c in coverage if c["method"] == m),
                         "na": sum(c["na"] for c in coverage if c["method"] == m)}
                     for m, _, _, _, _ in methods},
        "preview_auc": preview,
        **manifest_extra,
        "files": {},
    }
    for f in ["combined_method_axis_scores.csv", "combined_method_axis_wide.csv",
              "combined_method_axis_coverage.csv"]:
        p = os.path.join(B, f)
        manifest["files"][f] = {"path": os.path.relpath(p, ROOT), "sha256": sha256(p)}
    json.dump(manifest, open(os.path.join(B, "combined_method_axis_manifest.json"), "w"),
              indent=1, ensure_ascii=False)
    print(json.dumps(manifest["coverage"], indent=1))
    print(json.dumps({k: v["pooled_mean_auc"] for k, v in preview.items()}, indent=1))


if __name__ == "__main__":
    main()
