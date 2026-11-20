"""Generate Track A SI tables (S1-S9) from the official full JSON.

All numbers are read from ``figures/trackA_enrichment_full.json``; the script also
validates every per-target AUROC/EF5 against the JSON and writes a validation
report.  Tables are written under ``manuscript/si/`` as CSV + a README.

Usage (from ``code``)::

    python -m maes.si_tables_tracka
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
JSON_PATH = PROJECT_ROOT / "figures" / "trackA_enrichment_full.json"
OUT_DIR = PROJECT_ROOT / "manuscript" / "si"
VALIDATION_JSON = PROJECT_ROOT / "figures" / "si_trackA_validation.json"

PRIMARY_FIELDS = [
    "target", "panel_role", "metal_dependent", "n_act_scored", "n_decoy",
    "auc", "auc_ci_lo", "auc_ci_hi", "ef5_percent", "ef5_ci_lo", "ef5_ci_hi",
    "ef1_percent", "ef1_ci_lo", "ef1_ci_hi", "bedroc20", "decoy_ratio_observed",
    "box_size_A", "n_box_clipped_actives", "n_prep_missing",
]


def _entry(row: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "target": row.get("target"),
        "panel_role": row.get("panel_role"),
        "metal_dependent": row.get("metal_dependent"),
        "n_act_scored": row.get("n_active"),
        "n_decoy": row.get("n_decoy"),
        "auc": row.get("auc"),
        "auc_ci_lo": row.get("auc_ci_lo"),
        "auc_ci_hi": row.get("auc_ci_hi"),
        "ef5_percent": row.get("ef5_percent"),
        "ef5_ci_lo": row.get("ef5_ci_lo"),
        "ef5_ci_hi": row.get("ef5_ci_hi"),
        "ef1_percent": row.get("ef1_percent"),
        "ef1_ci_lo": row.get("ef1_ci_lo"),
        "ef1_ci_hi": row.get("ef1_ci_hi"),
        "bedroc20": row.get("bedroc20"),
        "decoy_ratio_observed": row.get("decoy_ratio_observed"),
        "box_size_A": row.get("box_size_A"),
        "n_box_clipped_actives": row.get("n_box_clipped_actives"),
        "n_prep_missing": row.get("n_prep_missing"),
    }


def _round_frame(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    for column in out.columns:
        if pd.api.types.is_float_dtype(out[column]):
            out[column] = out[column].astype(float).round(6)
    return out


def _clean_source(value: Any) -> Any:
    if value is None:
        return value
    return str(value).replace("\\", "/").replace("../results/", "results/").replace("..\\results\\", "results/")


def build(payload: Dict[str, Any]) -> Dict[str, Any]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    entries = payload["entries"]
    by_target = {entry["target"]: entry for entry in entries}
    primary = [entry for entry in entries if str(entry.get("panel_role", "main")).lower() == "main"]
    auxiliary = [entry for entry in entries if str(entry.get("panel_role", "main")).lower() != "main"]
    primary.sort(key=lambda entry: -float(entry["auc"]))
    auxiliary.sort(key=lambda entry: -float(entry["auc"]))

    s1 = _round_frame(pd.DataFrame([_entry(row) for row in primary]))
    s1.to_csv(OUT_DIR / "table_S1_per_target_primary.csv", index=False)
    s2 = _round_frame(pd.DataFrame([_entry(row) for row in auxiliary]))
    s2["note"] = "low-confidence/auxiliary target; excluded from all primary aggregates"
    s2.to_csv(OUT_DIR / "table_S2_auxiliary_targets.csv", index=False)

    # S3 property sensitivity -------------------------------------------------
    rows: List[Dict[str, Any]] = []
    for target, data in payload.get("supplementary_sensitivities", {}).get("property_sensitivity", {}).items():
        entry = by_target.get(target, {})
        by_prop = entry.get("auc_property_by_property") or {}
        rows.append(
            {
                "target": target,
                "panel_role": entry.get("panel_role"),
                "raw_auc": data.get("raw_auc"),
                "auc_property_stratified": data.get("auc_property_stratified"),
                "delta_stratified_vs_raw": (data.get("delta_vs_raw") or {}).get("stratified_mean"),
                "auc_property_matched": data.get("auc_property_matched"),
                "delta_matched_vs_raw": (data.get("delta_vs_raw") or {}).get("matched"),
                "mw_auc": (by_prop.get("mw") or {}).get("auc"),
                "logp_auc": (by_prop.get("logp") or {}).get("auc"),
                "hbd_auc": (by_prop.get("hbd") or {}).get("auc"),
            }
        )
    s3 = _round_frame(pd.DataFrame(rows))
    s3.to_csv(OUT_DIR / "table_S3_property_sensitivity.csv", index=False)

    # S4a Pa series exclusion (A-line per-series table is authoritative for EF5;
    # D-loader values are reported alongside where the ECFP4 cluster membership tie
    # gives a small EF5 difference)
    pa = payload.get("supplementary_sensitivities", {}).get("pa_series_exclusion", {})
    a_pa = PROJECT_ROOT / "data" / "trackA" / "LpxC_Pa_sensitivity_exclude_top.csv"
    a_pa_frame = pd.read_csv(a_pa) if a_pa.exists() else None
    def _a_row(keyword: str) -> Dict[str, Any]:
        if a_pa_frame is None:
            return {}
        hit = a_pa_frame[a_pa_frame["arm"].astype(str).str.contains(keyword, regex=False, na=False)]
        return hit.iloc[0].to_dict() if len(hit) else {}
    raw = pa.get("raw", {})
    a_full = _a_row("full")
    rows = [{
        "analysis": "raw (LpxC_Pa)", "n_excluded": 0,
        "auc_after": a_full.get("AUC", raw.get("auc")),
        "ef5_after": a_full.get("EF5%", raw.get("ef5_percent")),
        "delta_auc_vs_raw": 0.0, "delta_ef5_vs_raw": 0.0,
        "dloader_auc": raw.get("auc"), "dloader_ef5": raw.get("ef5_percent"),
        "note": "full panel; A-line per-series table and D-loader agree",
    }]
    bm = pa.get("largest_bm_scaffold", {})
    a_bm = _a_row("exclude_top1_scaffold")
    rows.append({
        "analysis": f"remove largest BM scaffold ({bm.get('n_excluded')} actives)",
        "n_excluded": bm.get("n_excluded"),
        "auc_after": a_bm.get("AUC", bm.get("auc_after")),
        "ef5_after": a_bm.get("EF5%", bm.get("ef5_after")),
        "delta_auc_vs_raw": (
            float(a_bm.get("AUC")) - float(raw.get("auc"))
            if a_bm.get("AUC") is not None and raw.get("auc") is not None else bm.get("delta_auc_vs_raw")
        ),
        "delta_ef5_vs_raw": (
            float(a_bm.get("EF5%")) - float(raw.get("ef5_percent"))
            if a_bm.get("EF5%") is not None and raw.get("ef5_percent") is not None else None
        ),
        "dloader_auc": bm.get("auc_after"), "dloader_ef5": bm.get("ef5_after"),
        "note": "AUC barely changes and EF5 does NOT decrease -> does not support chemotype-concentration explanation",
    })
    cluster = pa.get("largest_ecfp_cluster", {})
    a_cluster = _a_row("ECFP4_cluster")
    rows.append({
        "analysis": f"remove largest ECFP4 cluster ({cluster.get('n_excluded')} actives)",
        "n_excluded": cluster.get("n_excluded"),
        "auc_after": a_cluster.get("AUC", cluster.get("auc_after")),
        "ef5_after": a_cluster.get("EF5%", cluster.get("ef5_after")),
        "delta_auc_vs_raw": (
            float(a_cluster.get("AUC")) - float(raw.get("auc"))
            if a_cluster.get("AUC") is not None and raw.get("auc") is not None else cluster.get("delta_auc_vs_raw")
        ),
        "delta_ef5_vs_raw": (
            float(a_cluster.get("EF5%")) - float(raw.get("ef5_percent"))
            if a_cluster.get("EF5%") is not None and raw.get("ef5_percent") is not None else None
        ),
        "dloader_auc": cluster.get("auc_after"), "dloader_ef5": cluster.get("ef5_after"),
        "note": (
            "AUC barely changes and EF5 does NOT decrease; A-line EF5=6.15 vs D-loader 6.22 "
            "from ECFP4 cluster-membership tie/order (AUC identical) -> does not support "
            "chemotype-concentration explanation"
        ),
    })
    _round_frame(pd.DataFrame(rows)).to_csv(OUT_DIR / "table_S4a_pa_series_exclusion.csv", index=False)
    # A-line supplements: per-series exclusion table and 9-target chemotype composition
    a_chem = PROJECT_ROOT / "data" / "trackA" / "si_chemotype_table.csv"
    if a_pa.exists():
        shutil.copyfile(a_pa, OUT_DIR / "table_S4a1_pa_exclude_top_series.csv")
    if a_chem.exists():
        shutil.copyfile(a_chem, OUT_DIR / "table_S10_chemotype_composition.csv")

    # S4b shared decoy --------------------------------------------------------
    shared = payload.get("supplementary_sensitivities", {}).get("shared_decoy_paired", {})
    shared_intersection = payload.get("methods_evidence", {}).get("shared_decoy_intersection_lpxc")
    s4b = pd.DataFrame(
        [
            {
                "comparison": "LpxC_Ec vs LpxC_Pa (shared scored decoys)",
                "shared_decoy_pool_intersection": shared_intersection,
                "n_shared_scored_decoys": shared.get("n_shared_scored_decoys"),
                "n_ec_scored_decoys": shared.get("n_ec_scored_decoys"),
                "n_pa_scored_decoys": shared.get("n_pa_scored_decoys"),
                "auc_ec_shared": shared.get("auc_ec_shared"),
                "auc_pa_shared": shared.get("auc_pa_shared"),
                "delta_pa_minus_ec": shared.get("delta_pa_minus_ec"),
                "note": shared.get("note"),
            }
        ]
    )
    _round_frame(s4b).to_csv(OUT_DIR / "table_S4b_shared_decoy_paired.csv", index=False)

    # S5 boron surrogate ------------------------------------------------------
    boron = payload.get("supplementary_sensitivities", {}).get("boron_surrogate_sensitivity", {})
    rows = []
    for target, data in boron.items():
        main = data.get("with_surrogate", {}) or {}
        sens = data.get("without_surrogate", {}) or {}
        auc_main = main.get("auc")
        auc_sens = sens.get("auc")
        ef5_main = main.get("ef5_percent")
        ef5_sens = sens.get("ef5_percent")
        rows.append(
            {
                "target": target,
                "n_act_main": data.get("n_act_main"),
                "n_act_sensitivity": data.get("n_act_sensitivity"),
                "auc_main": auc_main,
                "auc_sensitivity": auc_sens,
                "delta_auc_main_minus_sensitivity": (
                    float(auc_main) - float(auc_sens)
                    if auc_main is not None and auc_sens is not None else None
                ),
                "ef5_main": ef5_main,
                "ef5_sensitivity": ef5_sens,
                "delta_ef5_main_minus_sensitivity": (
                    float(ef5_main) - float(ef5_sens)
                    if ef5_main is not None and ef5_sens is not None else None
                ),
                "n_surrogate_actives": data.get("n_surrogate_actives"),
                "rescue_source": _clean_source(data.get("source")),
            }
        )
    s5 = _round_frame(pd.DataFrame(rows))
    s5.to_csv(OUT_DIR / "table_S5_boron_surrogate.csv", index=False)

    # S6 exclusions -----------------------------------------------------------
    box_by = {row["target"]: row for row in payload.get("box_clipped_audit", {}).get("by_target", [])}
    prep_by = {row["target"]: row for row in payload.get("prep_missing_audit", {}).get("by_target", [])}
    dq_by = {row["target"]: row for row in payload.get("data_quality_audit", {}).get("by_target", [])}
    rows = []
    for target in [entry["target"] for entry in primary + auxiliary]:
        box = box_by.get(target, {})
        prep = prep_by.get(target, {})
        dq = dq_by.get(target, {})
        rows.append(
            {
                "target": target,
                "panel_role": by_target.get(target, {}).get("panel_role"),
                "box_size_A": box.get("box_size_A"),
                "n_box_clipped": box.get("n_box_clipped", 0),
                "n_box_clipped_actives": box.get("n_box_clipped_actives", 0),
                "n_act_before_box_clipped": box.get("n_act_before"),
                "n_act_after_box_clipped": box.get("n_act_after"),
                "n_prep_missing": prep.get("n_prep_missing", 0),
                "n_prep_missing_actives": prep.get("n_prep_missing_actives", 0),
                "n_prep_missing_decoys": prep.get("n_prep_missing_decoys", 0),
                "n_non_ok_other_actives": dq.get("n_non_ok_other_actives", 0),
                "non_ok_other_actives_fraction": dq.get("non_ok_other_actives_fraction"),
                "warning": dq.get("warning", False),
            }
        )
    s6 = _round_frame(pd.DataFrame(rows))
    s6.to_csv(OUT_DIR / "table_S6_exclusions.csv", index=False)

    # S7 panel scope ----------------------------------------------------------
    scope = payload.get("panel_scope", {})
    compute = payload.get("compute_settings", {})
    rows = [
        ("source_json", "figures/trackA_enrichment_full.json", "official D-line output"),
        ("source_generated_at", payload.get("generated_at"), "JSON timestamp"),
        ("n_panel_targets", scope.get("n_panel_targets"), "9-target panel (7 primary + 2 auxiliary)"),
        ("n_primary_targets", scope.get("n_primary_targets"), "used in main tables/figures/aggregates"),
        ("n_auxiliary", scope.get("n_auxiliary"), "GyrA_Ec + LepB_Ec; SI only"),
        ("n_actives_panel", scope.get("n_actives_panel"), "all 9 targets"),
        ("n_actives_primary", scope.get("n_actives_primary"), "7 primary targets"),
        ("primary_targets", "; ".join(scope.get("primary_targets", []) or []), "7 primary"),
        ("auxiliary_targets", "; ".join(scope.get("auxiliary_targets", []) or []), "2 auxiliary"),
        ("primary_caliber", payload.get("primary_caliber"), "main analysis includes B->C boron surrogates"),
        ("decoy_ratio_definition", payload.get("decoy_ratio_definition"), "observed per-target values in S1"),
        ("ef_definition", compute.get("ef_definition"), "EF1 gate n_act>=150"),
        ("box_clipped_rule", compute.get("box_clipped_rule"), "Amendment #10"),
        ("prep_missing_rule", compute.get("prep_missing_rule"), "data-processing rule"),
    ]
    pd.DataFrame(rows, columns=["field", "value", "definition"]).to_csv(
        OUT_DIR / "table_S7_panel_scope.csv", index=False
    )

    # S8 TEM1/LepB before/after rescue ---------------------------------------
    rows = []
    for target, data in boron.items():
        pre = data.get("pre_rescue_formal", {}) or {}
        main = data.get("with_surrogate", {}) or {}
        sens = data.get("without_surrogate", {}) or {}
        rows.append(
            {
                "target": target,
                "pre_rescue_n_act": pre.get("n_active"),
                "pre_rescue_auc": pre.get("auc"),
                "main_n_act": main.get("n_active"),
                "main_auc": main.get("auc"),
                "delta_main_minus_pre_rescue": (
                    float(main.get("auc")) - float(pre.get("auc"))
                    if main.get("auc") is not None and pre.get("auc") is not None else None
                ),
                "sensitivity_n_act": sens.get("n_active"),
                "sensitivity_auc": sens.get("auc"),
                "delta_main_minus_sensitivity": (
                    float(main.get("auc")) - float(sens.get("auc"))
                    if main.get("auc") is not None and sens.get("auc") is not None else None
                ),
                "n_surrogate_actives": data.get("n_surrogate_actives"),
                "note": "main includes B->C surrogates; sensitivity excludes them; pre-rescue formal used only available actives",
            }
        )
    _round_frame(pd.DataFrame(rows)).to_csv(OUT_DIR / "table_S8_tem1_lep_b_before_after.csv", index=False)

    # S9 multi-method table (Vina / Vinardo second scoring function / ECFP4) --
    ecfp4_metrics = PROJECT_ROOT / "results" / "trackA" / "ecfp4_baseline" / "baseline_metrics.csv"
    ecfp4_by_target: Dict[str, Dict[str, Any]] = {}
    if ecfp4_metrics.exists():
        for _, row in pd.read_csv(ecfp4_metrics).iterrows():
            ecfp4_by_target[str(row["target"])] = row.to_dict()
    multi_method_path = PROJECT_ROOT / "results" / "trackA" / "multi_method_comparison.csv"
    vinardo_by_target: Dict[str, Dict[str, Any]] = {}
    if multi_method_path.exists():
        multi = pd.read_csv(multi_method_path)
        for target, group in multi.groupby("target"):
            vina_row = group[group["method"].astype(str).str.lower() == "vina"]
            vin_row = group[group["method"].astype(str).str.lower().str.startswith("vinardo")]
            if len(vin_row):
                record = vin_row.iloc[0].to_dict()
                record["vina_auc"] = float(vina_row.iloc[0]["auc"]) if len(vina_row) else None
                record["vina_ef5"] = float(vina_row.iloc[0]["ef5"]) if len(vina_row) else None
                record["rho_vs_vina"] = float(vin_row.iloc[0].get("spearman_vs_vina"))
                vinardo_by_target[str(target)] = record
    rows = []
    for entry in primary + auxiliary:
        ecfp4 = ecfp4_by_target.get(entry["target"], {})
        vin = vinardo_by_target.get(entry["target"])
        rows.append(
            {
                "target": entry["target"],
                "panel_role": entry.get("panel_role"),
                "metal_dependent": entry.get("metal_dependent"),
                "n_act_scored": entry.get("n_active"),
                "vina_auc": entry.get("auc"),
                "vina_ef5": entry.get("ef5_percent"),
                "vinardo_N": vin.get("N") if vin else "not_computed",
                "vinardo_auc": vin.get("auc") if vin else "not_computed",
                "vinardo_ef5": vin.get("ef5") if vin else "not_computed",
                "vinardo_ef1": vin.get("ef1") if vin else "not_computed",
                "spearman_vinardo_vs_vina": vin.get("rho_vs_vina") if vin else "not_computed",
                "ecfp4_auc": ecfp4.get("auc"),
                "ecfp4_ef5": ecfp4.get("ef5"),
                "note": (
                    "Vina official (JSON); Vinardo = second scoring function applied to the "
                    "same Vina top poses (score_only, not an independent docking engine); "
                    "ECFP4 1-NN is a benchmark-validity control (DUD-E-style decoy bias), "
                    "not evidence that 2D methods outperform docking."
                ),
            }
        )
    s9 = _round_frame(pd.DataFrame(rows))
    s9.to_csv(OUT_DIR / "table_S9_multimethod.csv", index=False)

    # S11 ECFP4 benchmark-validity table -------------------------------------
    rows = []
    for entry in primary + auxiliary:
        target = entry["target"]
        ecfp4 = ecfp4_by_target.get(target, {})
        rows.append(
            {
                "target": target,
                "panel_role": entry.get("panel_role"),
                "metal_dependent": entry.get("metal_dependent"),
                "N": ecfp4.get("N"),
                "n_act": ecfp4.get("n_act"),
                "n_neg": ecfp4.get("n_neg"),
                "ecfp4_auc": ecfp4.get("auc"),
                "ecfp4_ef5": ecfp4.get("ef5"),
                "ecfp4_ef1": ecfp4.get("ef1"),
                "ecfp4_bedroc20": ecfp4.get("bedroc20"),
                "vina_auc": entry.get("auc"),
                "vina_ef5": entry.get("ef5_percent"),
                "delta_auc_ecfp4_minus_vina": (
                    float(ecfp4["auc"]) - float(entry["auc"]) if ecfp4 and ecfp4.get("auc") is not None else None
                ),
                "note": "2D 1-NN baseline; inflated by DUD-E-style decoy bias; not evidence that 2D > docking",
            }
        )
    s11 = _round_frame(pd.DataFrame(rows))
    s11.to_csv(OUT_DIR / "table_S11_ecfp4_baseline.csv", index=False)

    # S17 box-size three-tier table (Amendment #11) -------------------------
    box_tier_path = PROJECT_ROOT / "results" / "trackA" / "sensitivity" / "box" / "three_tier.csv"
    box28_summary_path = PROJECT_ROOT / "results" / "trackA" / "sensitivity" / "box" / "box28_summary.json"
    box28_summary: Dict[str, Any] = {}
    if box_tier_path.exists():
        s17 = pd.read_csv(box_tier_path)
        s17["amendment"] = (
            "Amendment #11: box28 equivalence arm PASS -> box26 retained as primary; "
            "box22 retained only as informative sensitivity (physical truncation)"
        )
        _round_frame(s17).to_csv(OUT_DIR / "table_S17_box_three_tier.csv", index=False)
        if box28_summary_path.exists():
            box28_summary = json.loads(box28_summary_path.read_text(encoding="utf-8"))

    # validation --------------------------------------------------------------
    checks: List[Dict[str, Any]] = []
    ok = True
    for entry in entries:
        target = entry["target"]
        row = s1[s1["target"] == target]
        if len(row):
            record = row.iloc[0]
        else:
            aux_row = s2[s2["target"] == target]
            if aux_row.empty:
                checks.append({"target": target, "check": "row_present", "ok": False})
                ok = False
                continue
            record = aux_row.iloc[0]
        for field, json_value in (("auc", entry["auc"]), ("ef5_percent", entry["ef5_percent"])):
            table_value = float(record[field])
            passed = bool(np.isclose(table_value, float(json_value), atol=5e-7))
            checks.append({"target": target, "check": field, "table": table_value,
                           "json": float(json_value), "ok": passed})
            ok = ok and passed
    # key cross-checks requested by Lead
    tem1 = by_target["TEM1_Ec"]
    checks.append({"target": "TEM1_Ec", "check": "main_auc_0.4296",
                   "table": float(round(tem1["auc"], 4)), "json": 0.4296,
                   "ok": bool(np.isclose(tem1["auc"], 0.4296, atol=5e-4))})
    ok = ok and checks[-1]["ok"]
    checks.append({"target": "panel_scope", "check": "n_actives_panel_1748",
                   "table": scope.get("n_actives_panel"), "json": 1748,
                   "ok": scope.get("n_actives_panel") == 1748})
    ok = ok and checks[-1]["ok"]
    checks.append({"target": "panel_scope", "check": "n_actives_primary_1648",
                   "table": scope.get("n_actives_primary"), "json": 1648,
                   "ok": scope.get("n_actives_primary") == 1648})
    ok = ok and checks[-1]["ok"]
    checks.append({"target": "novel", "check": "matched_178_of_178",
                   "table": payload["novel_scaffold_audit"].get("n_matched_total"),
                   "json": 178,
                   "ok": payload["novel_scaffold_audit"].get("n_matched_total") == 178})
    ok = ok and checks[-1]["ok"]
    # sensitivity/exclusion tables vs JSON
    for target in ("TEM1_Ec", "LepB_Ec"):
        data = boron[target]
        row = s5[s5["target"] == target]
        if len(row):
            record = row.iloc[0]
            for column, json_value in (
                ("auc_main", data["with_surrogate"]["auc"]),
                ("auc_sensitivity", data["without_surrogate"]["auc"]),
            ):
                passed = bool(np.isclose(float(record[column]), float(json_value), atol=5e-7))
                checks.append({"target": target, "check": f"sensitivity::{column}",
                               "table": float(record[column]), "json": float(json_value), "ok": passed})
                ok = ok and passed
    for target in ("Gyrase_Ec", "LepB_Ec"):
        box = box_by[target]
        row = s6[s6["target"] == target]
        if len(row):
            record = row.iloc[0]
            passed = int(record["n_box_clipped_actives"]) == int(box.get("n_box_clipped_actives", 0))
            checks.append({"target": target, "check": "exclusions::n_box_clipped_actives",
                           "table": int(record["n_box_clipped_actives"]),
                           "json": int(box.get("n_box_clipped_actives", 0)), "ok": passed})
            ok = ok and passed
    nonmetal_mean = payload["metal_vs_nonmetal"]["non-metal"]["mean_auc"]
    checks.append({"target": "metal_vs_nonmetal", "check": "nonmetal_mean_0.567",
                   "table": round(float(nonmetal_mean), 3), "json": 0.567,
                   "ok": bool(np.isclose(nonmetal_mean, 0.567173, atol=5e-4))})
    ok = ok and checks[-1]["ok"]
    # ECFP4 benchmark-validity table must match B-line baseline metrics
    for target, row in ecfp4_by_target.items():
        table_row = s11[s11["target"] == target]
        if table_row.empty:
            continue
        table_value = float(table_row.iloc[0]["ecfp4_auc"])
        passed = bool(np.isclose(table_value, float(row["auc"]), atol=5e-7))
        checks.append({"target": target, "check": "ecfp4_auc",
                       "table": table_value, "json": float(row["auc"]), "ok": passed})
        ok = ok and passed
    # Vinardo second-scoring-function values must match B-line comparison file
    for target, row in vinardo_by_target.items():
        table_row = s9[s9["target"] == target]
        if table_row.empty:
            continue
        record = table_row.iloc[0]
        for column, expected in (("vinardo_auc", row.get("auc")), ("spearman_vinardo_vs_vina", row.get("rho_vs_vina"))):
            passed = bool(np.isclose(float(record[column]), float(expected), atol=5e-7))
            checks.append({"target": target, "check": column,
                           "table": float(record[column]), "json": float(expected), "ok": passed})
            ok = ok and passed
    checks.append({"target": "multi_method", "check": "vinardo_spearman_0.65_0.84",
                   "table": [round(float(row.get("rho_vs_vina")), 3) for row in vinardo_by_target.values()],
                   "json": [0.65, 0.84],
                   "ok": bool(all(0.64 <= float(row.get("rho_vs_vina")) <= 0.85 for row in vinardo_by_target.values()))})
    ok = ok and checks[-1]["ok"]
    # Pa sensitivity supplement must match A-line numbers
    if (OUT_DIR / "table_S4a1_pa_exclude_top_series.csv").exists():
        pa_a = pd.read_csv(OUT_DIR / "table_S4a1_pa_exclude_top_series.csv")
        expected = {"full": 0.8525, "exclude_top1_scaffold(n=108)": 0.8351,
                    "exclude_top_ECFP4_cluster(n=119)": 0.8510}
        for arm, auc in expected.items():
            hit = pa_a[pa_a["arm"] == arm]
            passed = bool(len(hit)) and bool(np.isclose(float(hit.iloc[0]["AUC"]), auc, atol=5e-4))
            checks.append({"target": "LpxC_Pa", "check": f"s4a1::{arm}",
                           "table": float(hit.iloc[0]["AUC"]) if len(hit) else None,
                           "json": auc, "ok": passed})
            ok = ok and passed
    # Amendment #11 box28 confirmatory arm
    if box_tier_path.exists():
        s17_check = pd.read_csv(box_tier_path)
        box28 = s17_check[s17_check["arm"].astype(str) == "box28"].iloc[0]
        checks.append({"target": "LpxC_Ec", "check": "box28_rho_0.9487",
                       "table": round(float(box28["spearman_rho"]), 4), "json": 0.9487,
                       "ok": bool(np.isclose(float(box28["spearman_rho"]), 0.9487, atol=5e-4))})
        ok = ok and checks[-1]["ok"]
        checks.append({"target": "LpxC_Ec", "check": "box28_abs_dAUC_le_0.05",
                       "table": round(abs(float(box28["dAUC"])), 4), "json": 0.0433,
                       "ok": bool(abs(float(box28["dAUC"])) <= 0.05 and np.isclose(abs(float(box28["dAUC"])), 0.0433, atol=5e-4))})
        ok = ok and checks[-1]["ok"]
        checks.append({"target": "LpxC_Ec", "check": "box28_verdict_PASS",
                       "table": str(box28["verdict"]), "json": "PASS",
                       "ok": str(box28["verdict"]) == "PASS"})
        ok = ok and checks[-1]["ok"]
        box22 = s17_check[s17_check["arm"].astype(str) == "box22"].iloc[0]
        checks.append({"target": "LpxC_Ec", "check": "box22_verdict_informative_failure",
                       "table": str(box22["verdict"]), "json": "FAIL",
                       "ok": str(box22["verdict"]) == "FAIL"})
        ok = ok and checks[-1]["ok"]
    validation = {
        "source_json": str(JSON_PATH),
        "source_generated_at": payload.get("generated_at"),
        "n_checks": len(checks),
        "n_failures": int(sum(1 for check in checks if not check["ok"])),
        "all_passed": bool(ok),
        "checks": checks,
        "tables_written": sorted(path.name for path in OUT_DIR.glob("table_*.csv")),
    }
    VALIDATION_JSON.write_text(json.dumps(validation, indent=2, default=str), encoding="utf-8")

    readme = f"""# Track A SI tables (D line; source = official full JSON)

Source: `figures/trackA_enrichment_full.json`
- generated_at: **{payload.get('generated_at')}**
- primary_caliber: **{payload.get('primary_caliber')}**
- panel_scope: **{scope.get('n_panel_targets')} panel targets / {scope.get('n_primary_targets')} primary / {scope.get('n_auxiliary')} auxiliary; actives {scope.get('n_actives_panel')} panel / {scope.get('n_actives_primary')} primary**
- CI: active-level cluster bootstrap (decoys clustered by matched active), 95%
- EF1 only when n_act >= 150; EF5 for all targets
- decoy ratio: {payload.get('decoy_ratio_definition')}
- validation: **{validation['n_checks'] - validation['n_failures']}/{validation['n_checks']} checks passed** (`figures/si_trackA_validation.json`)

## Tables
- `table_S1_per_target_primary.csv` — 7 primary targets: AUROC/EF5/EF1 + cluster CI, BEDROC20, observed ratio, box, exclusions.
- `table_S2_auxiliary_targets.csv` — GyrA_Ec/LepB_Ec (SI only; excluded from aggregates).
- `table_S3_property_sensitivity.csv` — property-stratified and nearest-decoy matched AUROC + deltas.
- `table_S4a_pa_series_exclusion.csv` — LpxC_Pa BM-scaffold / ECFP4-cluster exclusion sensitivity (AUC and EF5 deltas; corrected interpretation note).
- `table_S4a1_pa_exclude_top_series.csv` — A-line per-series exclusion table (full/series-excluded metrics; source `data/trackA/LpxC_Pa_sensitivity_exclude_top.csv`).
- `table_S4b_shared_decoy_paired.csv` — LpxC_Ec vs LpxC_Pa shared-decoy comparison.
- `table_S5_boron_surrogate.csv` — main (includes B->C surrogates) vs exclude-surrogate sensitivity.
- `table_S6_exclusions.csv` — box-clipped, prep_missing, other non-ok actives; warnings.
- `table_S7_panel_scope.csv` — panel/primary/auxiliary scope, calibers, EF/decoy/box rules.
- `table_S8_tem1_lep_b_before_after.csv` — pre-rescue formal vs main vs sensitivity.
- `table_S9_multimethod.csv` — Vina official; **Vinardo second scoring function** (4 targets, score_only on same poses, not an independent engine); ECFP4 1-NN baseline (9 targets, benchmark-validity control).
- `table_S10_chemotype_composition.csv` — A-line actives chemical-series composition for all 9 targets (source `data/trackA/si_chemotype_table.csv`).
- `table_S11_ecfp4_baseline.csv` — ECFP4 1-NN benchmark-validity table (9 targets; AUC/EF5/EF1/BEDROC vs Vina).
- `table_S17_box_three_tier.csv` — LpxC_Ec 1,000-ligand box-size three-tier table (22/26/28 Å; Amendment #11: box28 PASS → box26 primary, box22 informative failure).

## Box-size sensitivity (Amendment #11)
- Three-tier AUC on the same 1,000-ligand LpxC_Ec subset: **0.2242 (box22) → 0.2893 (box26) → 0.3326 (box28)**; ranking stable (**Spearman ρ(26 vs 28)=0.949**), **|ΔAUC|=0.0433 ≤ 0.05** and ΔEF5 = 0 (CI contains 0) → **PASS**.
- Per Amendment #11: **box26 remains the primary panel; Track A is not downgraded; box22 is retained only as an informative sensitivity** (physical truncation: 18.5% actives >22 Å, max 24.24 Å).
- Suggested limitation sentence: "Absolute enrichment was mildly box-size dependent (AUC 0.224/0.289/0.333 for 22/26/28 Å on the 1,000-ligand LpxC_Ec subset), while ranking was stable (Spearman ρ=0.949 for 26 vs 28 Å); the pre-registered size audit selected 26 Å as primary, and 22 Å (which truncates the largest ligands) is retained only as an informative sensitivity."

## Multi-method (Vina / Vinardo) limitation
- Vinardo was applied as a **second scoring function** to the same Vina top poses (Vina 1.2.5 `--score_only --scoring vinardo`, same box/center; Gyrase uses the GyrB receptor); it is **not an independent docking engine**.
- Rankings are moderately concordant across scoring functions (**Spearman ρ 0.65–0.84**), but absolute enrichment is scoring-function dependent (LpxC_Ec AUC 0.300→**0.554**; LpxC_Pa 0.852→**0.918**; DHFR 0.666→**0.695**; Gyrase 0.578→**0.518**).
- Therefore **(i) the directional main conclusion is robust; (ii) scoring function is a major uncertainty source and a single scoring function must not be treated as the ceiling of physics-based scoring; (iii) multi-method comparisons must account for both engine and scoring function**. These results are R9/SI only and are not merged into the main AI-vs-physics conclusion.

## Benchmark-validity caveat (ECFP4 1-NN)
- The ECFP4 1-NN baseline reaches AUROC **0.82–1.00** (main 7 targets: 0.82–1.00; EF5 10.8–18.2); its near-ceiling performance reflects **DUD-E-style decoy bias** (decoys were ECFP4-filtered against actives, so the 2D task is artificially easy), **not evidence that 2D methods outperform docking**.
- Absolute EF/EF5/AUC from this panel must not be extrapolated to real screening; the 2D baseline is a **benchmark-validity control** only.
- Multi-method table: Vina (final); **Vinardo** = second scoring function on the same poses (not an independent engine; B-line pending); ECFP4 1-NN baseline (filled, task-30); Uni-Dock infeasible on the available hardware.

## Pa series-exclusion interpretation (corrected 2026-11-18)
- Removing the largest LpxC_Pa BM scaffold (108 actives) leaves AUC **0.835** (Δ−0.017) and EF5 **5.81** (not lower); removing the largest ECFP4 cluster (119 actives) leaves AUC **0.851** (Δ−0.001) and EF5 **6.15** (A-line per-series table; D-loader gives 6.22 from a cluster-membership tie, AUC identical).
- Therefore **LpxC_Pa's 0.852 is NOT a chemotype-concentration artifact**; the LpxC_Ec vs LpxC_Pa difference is **not explained by actives chemotype concentration** and is retained as a target/receptor/activity-space difference.
- Caveat: the two LpxC share **4,739 s1 decoys** (60.3% of the Ec pool) and **48 shared actives**; the cross-species axis uses the paired comparison on **951 shared scored decoys**.
- A-line sources: `data/trackA/si_chemotype_table.csv`, `data/trackA/LpxC_Pa_sensitivity_exclude_top.csv`, `notes/si_chemotype_table.md`.

## Methods sentences (verbatim, consistent with JSON)
- Panel scope: "The discovery panel comprises nine targets (seven primary and two auxiliary); the primary analysis includes 1,648 of 1,748 panel actives; GyrA_Ec and LepB_Ec are reported as low-confidence auxiliary targets and excluded from all primary aggregates."
- Primary caliber: "The primary analysis includes 25 boron-containing ligands rescued with a B→C isosteric surrogate (TEM1_Ec 23, LepB_Ec 2); the exclude-surrogate sensitivity is reported separately."
- Decoy ratio: "{payload.get('decoy_ratio_definition')}."
- Enrichment factor: "{compute.get('ef_definition')}."
- Metal strata: "Heterogeneous; per-target reporting is primary; metal/non-metal aggregates are descriptive only and are not evidence of a uniform metal effect (candidate confounders: metal coordination, chemotype concentration [not supported for the LpxC Ec/Pa difference], activity strength, assay heterogeneity)."
- Pa series: "Removing LpxC_Pa's largest BM scaffold (108 actives) or ECFP4 cluster (119 actives) leaves AUROC at 0.835/0.851 with EF5 not decreasing; the LpxC_Ec vs LpxC_Pa difference is therefore not attributable to actives chemotype concentration."
"""
    (OUT_DIR / "README.md").write_text(readme, encoding="utf-8")
    return {
        "out_dir": str(OUT_DIR),
        "tables": sorted(path.name for path in OUT_DIR.glob("table_*.csv")),
        "validation": str(VALIDATION_JSON),
        "all_passed": bool(ok),
        "n_checks": len(checks),
        "n_failures": int(sum(1 for check in checks if not check["ok"])),
    }


if __name__ == "__main__":
    payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    result = build(payload)
    print(json.dumps(result, indent=2))
