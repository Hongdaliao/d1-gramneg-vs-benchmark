"""LIT-PCBA anchor SI tables (S15 AC table, S16 per-set metrics, S16b subseed).

Anchor results are an **external consistency check (appendix)** after the
pre-registered AC2/AC3 failure; this script records both the pre-registered
9-evaluable-set primary judgement and the descriptive 15-set numbers without
cherry-picking.  All values are read from C/B-line anchor outputs.

Usage (from ``code``)::

    python -m maes.si_tables_anchor
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
ANCHOR = PROJECT_ROOT / "results" / "litpcba"
SI_DIR = PROJECT_ROOT / "manuscript" / "si"
VALIDATION_JSON = PROJECT_ROOT / "figures" / "si_anchor_validation.json"


def build() -> Dict[str, Any]:
    SI_DIR.mkdir(parents=True, exist_ok=True)
    ac = json.loads((ANCHOR / "anchor_AC.json").read_text(encoding="utf-8"))
    metrics = pd.read_csv(ANCHOR / "anchor_metrics.csv")
    subseed = pd.read_csv(ANCHOR / "anchor_metrics_by_subseed.csv")
    integrity = json.loads((ANCHOR / "anchor_integrity_check.json").read_text(encoding="utf-8"))

    n_sets = len(metrics)
    n_done = int((metrics["status"].astype(str) == "done").sum())
    ac2_eligible = metrics["n_act_val"].astype(float) >= 20
    n_eligible = int(ac2_eligible.sum())
    n_ef5_ge2_eligible = int(((metrics["ef5"] >= 2) & ac2_eligible).sum())
    ac2_primary = float(n_ef5_ge2_eligible / n_eligible) if n_eligible else float("nan")
    n_ef5_ge2_all = int((metrics["ef5"] >= 2).sum())
    ac2_descriptive = float(n_ef5_ge2_all / n_sets) if n_sets else float("nan")
    median_auc = float(metrics["auc"].median())
    median_ef5 = float(metrics["ef5"].median())
    max_cv = float(metrics["ef5_cv_pct"].max())
    max_cv_set = str(metrics.loc[metrics["ef5_cv_pct"].idxmax(), "set"])
    ef1_sets = metrics[metrics["n_act_ge150"].astype(bool)][["set", "ef1"]].to_dict(orient="records")

    ac_rows = [
        {
            "AC": "AC1",
            "criterion": "15/15 attempted; >=10/15 end-to-end completed",
            "result": f"{ac['AC1_attempted']}/15 attempted, {ac['AC1_completed']}/15 completed; {n_done}/15 rows status done",
            "status": "PASS" if ac["AC1_attempted"] == 15 and ac["AC1_completed"] == 15 and n_done == 15 else "FAIL",
            "notes": "no target cherry-picking; all 15 sets attempted",
        },
        {
            "AC": "AC2",
            "criterion": ">=30% of pre-registered evaluable sets (n_val_act>=20) have EF5%>=2",
            "result": f"{n_ef5_ge2_eligible}/{n_eligible} = {ac2_primary*100:.1f}%",
            "status": "FAIL",
            "notes": "primary pre-registered judgement; evaluable sets = " + ", ".join(sorted(metrics.loc[ac2_eligible, "set"].astype(str))),
        },
        {
            "AC": "AC3",
            "criterion": "EF5% coefficient of variation across 3 inactive subseeds <=10%",
            "result": f"max CV {max_cv:.2f}% ({max_cv_set})",
            "status": "FAIL",
            "notes": "CV computed over seed 42/43/44 inactive subsamples",
        },
        {
            "AC": "AC4",
            "criterion": "all 15 rows retained, failures with failure_reason",
            "result": f"15/15 rows; all status=done; integrity checks complete for {len(integrity.get('sets', {}))} sets",
            "status": "PASS",
            "notes": "no failed/omitted set",
        },
        {
            "AC": "AC5",
            "criterion": "AI arm >=5 targets if run; otherwise CPU-only",
            "result": "n/a — no AI arm run; CPU-only anchor",
            "status": "n/a (CPU-only)",
            "notes": "CPU-only anchor; not mixed into AI-vs-physics conclusions",
        },
        {
            "AC": "AC6",
            "criterion": "protocol freeze timestamp/hash precedes any anchor result",
            "result": "freeze 2026-11-16T21:54:28+08:00 precedes results (2026-11-18+)",
            "status": "PASS",
            "notes": "AC6 record in notes/review/litpcba_anchor_amendments.md",
        },
        {
            "AC": "AC7",
            "criterion": "one-command reproduction of a single set",
            "result": "demo command recorded; deterministic seeds",
            "status": "PASS",
            "notes": "ADRB2 demo re-run documented",
        },
    ]
    s15 = pd.DataFrame(ac_rows)
    s15.to_csv(SI_DIR / "table_S15_anchor_AC.csv", index=False)

    s16 = metrics.copy()
    s16["ac2_eligible"] = ac2_eligible.astype(bool)
    s16["descriptive_15set_only"] = ~s16["ac2_eligible"]
    keep = [
        "set", "n_act_val", "N_primary", "auc", "ef5", "ef1", "bedroc20",
        "ef5_mean_3sub", "ef5_cv_pct", "ef5_ge2", "ef1_ge2", "n_act_ge150",
        "ac2_eligible", "descriptive_only", "note", "status",
    ]
    s16 = s16[[column for column in keep if column in s16.columns]].round(6)
    s16.to_csv(SI_DIR / "table_S16_anchor_per_set.csv", index=False)
    subseed.round(6).to_csv(SI_DIR / "table_S16b_anchor_subseed.csv", index=False)

    validation = {
        "source_dir": str(ANCHOR),
        "n_checks": 12,
        "n_failures": 0,
        "all_passed": True,
        "checks": {},
    }
    checks: List[Dict[str, Any]] = []
    def check(name: str, table: Any, expected: Any) -> None:
        passed = bool(table == expected)
        checks.append({"check": name, "table": table, "expected": expected, "ok": passed})
    check("n_sets_15", n_sets, 15)
    check("all_attempted_completed", (ac["AC1_attempted"], ac["AC1_completed"]), (15, 15))
    check("AC2_eligible_9", n_eligible, 9)
    check("AC2_ef5_ge2_2", n_ef5_ge2_eligible, 2)
    check("AC2_proportion_0.222", round(ac2_primary, 3), 0.222)
    check("AC2_FAIL", ac["AC2_pass_ge_30pct"], False)
    check("AC3_max_cv_28.28", round(max_cv, 2), 28.28)
    check("AC3_FAIL", ac["AC3_pass_le_10pct"], False)
    check("descriptive_15_ef5_ge2_5", n_ef5_ge2_all, 5)
    check("descriptive_15_proportion_0.333", round(ac2_descriptive, 3), 0.333)
    check("median_auc_0.547", round(median_auc, 3), 0.547)
    check("median_ef5_1.19", round(median_ef5, 2), 1.19)
    validation["checks"] = {row["check"]: row["ok"] for row in checks}
    validation["n_checks"] = len(checks)
    validation["n_failures"] = int(sum(1 for row in checks if not row["ok"]))
    validation["all_passed"] = validation["n_failures"] == 0
    validation["ac_details"] = ac
    VALIDATION_JSON.write_text(json.dumps(validation, indent=2, default=str), encoding="utf-8")

    readme = f"""# LIT-PCBA external anchor — SI (appendix consistency check)

Source: `results/litpcba/anchor_AC.json`, `anchor_metrics.csv`, `anchor_metrics_by_subseed.csv`, `anchor_integrity_check.json`.
Validation: **{validation['n_checks'] - validation['n_failures']}/{validation['n_checks']} checks passed** (`figures/si_anchor_validation.json`).

## Positioning (read before using numbers)

- The anchor is an **external consistency check (appendix)**, **not external validation**. It does **not** satisfy the pre-registered anchor criteria (AC2 FAIL, AC3 FAIL), so it is excluded from the main statistical-power argument and from any validation claim.
- The anchor was **CPU-only** (no AI arm; AC5 n/a) and does not enter AI-vs-physics conclusions.

## Pre-registered AC judgement (primary)

| AC | result | status |
|---|---|---|
| AC1 | 15/15 attempted and completed | **PASS** |
| AC2 | **2/9 = 22.2% < 30%** (EF5%>=2, evaluable n_val>=20 sets) | **FAIL** |
| AC3 | **max CV 28.3% (MTORC1) > 10%** | **FAIL** |
| AC4 | 15/15 rows, no omitted/failed set | PASS |
| AC5 | n/a — CPU-only, no AI arm | n/a |
| AC6 | freeze before results | PASS |
| AC7 | one-command demo | PASS |

## Descriptive 15-set numbers (not used for pass/fail or cherry-picking)

- EF5%>=2: **5/15 = 33.3%** (GBA, KAT2A, OPRK1, PPARG, TP53) — **descriptive only; the primary AC2 judgement uses the 9 pre-registered evaluable sets (22.2%, FAIL)**.
- median AUC **0.547**, median EF5 **1.19**; EF1 only ALDH1 **1.21**, VDR **1.28**.
- AC2 eligible sets: ALDH1, ESR1_ant, FEN1, GBA, KAT2A, MAPK1, MTORC1, PKM2, VDR.

## Methods (anchor)

- Amendment A2 deviations from the original full anchor: inactives:active **20:1**, cap **1,000**/set (`min(20xn_val_act,1000,pool)`), **3 inactive subsample seeds**, Vina **exhaustiveness 4**, Vina seed 42; **15/15 sets attempted** (no target selection).
- AC2 primary judgement on the **9 sets with n_val_act>=20**; small sets descriptive only; EF1 only for ALDH1/VDR (n_val>=150); ADRB2 (4 actives) descriptive-only chain validation.
- The anchor is reported as an **appendix external consistency check**; it does not change the main target-specific conclusions.

## Tables

- `table_S15_anchor_AC.csv` — AC1–AC7 judgement.
- `table_S16_anchor_per_set.csv` — 15-set per-set metrics, AC2 eligibility/descriptive flags.
- `table_S16b_anchor_subseed.csv` — 45 rows (15 sets × 3 inactive subseeds) for AC3.
"""
    (SI_DIR / "README_anchor.md").write_text(readme, encoding="utf-8")
    return {
        "si_dir": str(SI_DIR),
        "tables": ["table_S15_anchor_AC.csv", "table_S16_anchor_per_set.csv", "table_S16b_anchor_subseed.csv"],
        "validation": str(VALIDATION_JSON),
        "all_passed": validation["all_passed"],
        "n_checks": validation["n_checks"],
        "n_failures": validation["n_failures"],
    }


if __name__ == "__main__":
    print(json.dumps(build(), indent=2))
