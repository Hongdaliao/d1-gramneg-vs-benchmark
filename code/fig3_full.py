"""Figure 3 (final) — Track A enrichment across the 9-target panel.

Rows = targets (7 primary ordered by AUC; 2 auxiliary shown greyed and excluded
from aggregates); columns = metal annotation, AUROC (with active-level cluster
bootstrap 95% CI), EF5% (with CI).  Outputs 300-dpi PNG + vector PDF plus a
machine-readable report that records the source JSON timestamp and the values
used (so the figure can be audited against the JSON).

Usage (from ``code``)::

    python -m maes.fig3_full
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Patch, Rectangle

PROJECT_ROOT = Path(__file__).resolve().parents[2]
FIGURES_DIR = PROJECT_ROOT / "figures"
JSON_PATH = FIGURES_DIR / "trackA_enrichment_full.json"

METAL_COLORS = {"Zn": "#E67E22", "Mg": "#8E44AD", "none": "#2E86C1", "unknown": "#7F8C8D"}


def _fmt(value: float, digits: int = 2) -> str:
    return "—" if value is None or not np.isfinite(float(value)) else f"{float(value):.{digits}f}"


def load_payload() -> Dict[str, Any]:
    if not JSON_PATH.exists():
        raise FileNotFoundError(f"official Track A JSON not found: {JSON_PATH}")
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


def ordered_entries(payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    entries = list(payload["entries"])
    primary = [entry for entry in entries if str(entry.get("panel_role", "main")).lower() == "main"]
    auxiliary = [entry for entry in entries if str(entry.get("panel_role", "main")).lower() != "main"]
    primary.sort(key=lambda entry: -float(entry["auc"]))
    auxiliary.sort(key=lambda entry: -float(entry["auc"]))
    for entry in primary:
        entry["_row_kind"] = "primary"
    for entry in auxiliary:
        entry["_row_kind"] = "auxiliary"
    return primary + auxiliary


def make_figure(json_path: Path = JSON_PATH) -> Dict[str, Any]:
    global JSON_PATH
    JSON_PATH = json_path
    payload = load_payload()
    entries = ordered_entries(payload)
    primary = [entry for entry in entries if entry["_row_kind"] == "primary"]
    auxiliary = [entry for entry in entries if entry["_row_kind"] == "auxiliary"]
    n_primary = len(primary)
    n_aux = len(auxiliary)

    labels = [
        f"{entry['target']}   {entry.get('metal_dependent', '?')}   n={entry['n_active']}"
        + ("   [aux; SI]" if entry["_row_kind"] == "auxiliary" else "")
        for entry in entries
    ]
    auc = np.array([float(entry["auc"]) for entry in entries])
    auc_lo = np.array([float(entry.get("auc_ci_lo", np.nan)) for entry in entries])
    auc_hi = np.array([float(entry.get("auc_ci_hi", np.nan)) for entry in entries])
    ef5 = np.array([float(entry["ef5_percent"]) for entry in entries])
    ef5_lo = np.array([float(entry.get("ef5_ci_lo", np.nan)) for entry in entries])
    ef5_hi = np.array([float(entry.get("ef5_ci_hi", np.nan)) for entry in entries])
    ef1 = [
        float(entry["ef1_percent"]) if entry.get("ef1_percent") is not None else np.nan
        for entry in entries
    ]
    bedroc = np.array([float(entry["bedroc20"]) for entry in entries])

    fig = plt.figure(figsize=(11.2, 5.9))
    grid = fig.add_gridspec(
        2, 4, height_ratios=[1.0, 0.075], width_ratios=[0.32, 1.05, 0.92, 0.035],
        left=0.235, right=0.965, top=0.82, bottom=0.26, wspace=0.20, hspace=0.42,
    )
    ax_metal = fig.add_subplot(grid[0, 0])
    ax_auc = fig.add_subplot(grid[0, 1], sharey=ax_metal)
    ax_ef5 = fig.add_subplot(grid[0, 2], sharey=ax_metal)
    ax_cb_auc = fig.add_subplot(grid[1, 1])
    ax_cb_ef5 = fig.add_subplot(grid[1, 2])

    rows = np.arange(len(entries))
    # metal strip ---------------------------------------------------------
    for row, entry in zip(rows, entries):
        metal = str(entry.get("metal_dependent", "unknown"))
        ax_metal.add_patch(
            Rectangle((0.0, row - 0.5), 1.0, 1.0,
                      facecolor=METAL_COLORS.get(metal, METAL_COLORS["unknown"]),
                      edgecolor="white", linewidth=1.2)
        )
        ax_metal.text(0.5, row, metal, ha="center", va="center", fontsize=8.5,
                      color="white", fontweight="bold")
    ax_metal.set_xlim(0, 1)
    ax_metal.set_ylim(len(entries) - 0.5, -0.5)
    ax_metal.set_xticks([])
    ax_metal.set_yticks(rows)
    ax_metal.set_yticklabels(labels, fontsize=8.5)
    for tick, entry in zip(ax_metal.get_yticklabels(), entries):
        if entry["_row_kind"] == "auxiliary":
            tick.set_color("0.45")
            tick.set_fontstyle("italic")
    ax_metal.set_title("metal", fontsize=9, pad=6)

    # AUC heatmap ---------------------------------------------------------
    norm_auc = TwoSlopeNorm(vmin=0.2, vcenter=0.5, vmax=0.9)
    im_auc = ax_auc.imshow(
        auc.reshape(-1, 1), cmap="RdBu_r", norm=norm_auc, aspect="auto"
    )
    # EF5 heatmap ---------------------------------------------------------
    norm_ef5 = Normalize(vmin=0.0, vmax=6.0)
    im_ef5 = ax_ef5.imshow(
        ef5.reshape(-1, 1), cmap="YlOrRd", norm=norm_ef5, aspect="auto"
    )

    for row, entry in enumerate(entries):
        is_aux = entry["_row_kind"] == "auxiliary"
        grey = dict(facecolor="0.92", alpha=0.82, edgecolor="none", zorder=3)
        if is_aux:
            for axis in (ax_auc, ax_ef5):
                axis.add_patch(Rectangle((-0.5, row - 0.5), 1.0, 1.0, **grey))
        auc_text = f"{auc[row]:.2f}"
        auc_ci_text = f"{auc_lo[row]:.2f}–{auc_hi[row]:.2f}"
        ax_auc.text(0.0, row - 0.14, auc_text, ha="center", va="center", fontsize=9.5,
                    fontweight="bold", color="0.25" if is_aux else "black", zorder=5)
        ax_auc.text(0.0, row + 0.20, auc_ci_text, ha="center", va="center", fontsize=7.2,
                    color="0.35", zorder=5)
        ax_ef5.text(0.0, row - 0.14, f"{ef5[row]:.2f}", ha="center", va="center", fontsize=9.5,
                    fontweight="bold", color="0.25" if is_aux else "black", zorder=5)
        ax_ef5.text(0.0, row + 0.20, f"{ef5_lo[row]:.2f}–{ef5_hi[row]:.2f}", ha="center",
                    va="center", fontsize=7.2, color="0.35", zorder=5)

    ax_auc.set_xticks([0])
    ax_auc.set_xticklabels(["AUROC\n(95% cluster CI)"], fontsize=9)
    ax_ef5.set_xticks([0])
    ax_ef5.set_xticklabels(["EF5%\n(95% cluster CI)"], fontsize=9)
    for axis in (ax_auc, ax_ef5):
        axis.tick_params(axis="y", left=False, labelleft=False)
        axis.set_ylim(len(entries) - 0.5, -0.5)
    ax_auc.set_title("AUROC (chance = 0.5)", fontsize=10, pad=8)
    ax_ef5.set_title("EF5%", fontsize=10, pad=8)

    # auxiliary separator; aux rows are labelled in the ytick labels
    if n_aux:
        boundary = n_primary - 0.5
        ax_auc.axhline(boundary, color="0.35", lw=1.1, ls="--")
        ax_ef5.axhline(boundary, color="0.35", lw=1.1, ls="--")

    cbar_auc = fig.colorbar(im_auc, cax=ax_cb_auc, orientation="horizontal")
    cbar_auc.set_label("AUROC (white = 0.5 chance)", fontsize=7.5)
    cbar_auc.ax.tick_params(labelsize=7)
    cbar_ef5 = fig.colorbar(im_ef5, cax=ax_cb_ef5, orientation="horizontal")
    cbar_ef5.set_label("EF5%", fontsize=7.5)
    cbar_ef5.ax.tick_params(labelsize=7)

    footer_note = (
        "Primary caliber includes 25 boron-surrogate ligands (TEM1 23, LepB 2; B→C isosteric). "
        "GyrA_Ec and LepB_Ec are auxiliary (grey/italic) and excluded from all primary aggregates. "
        "AUROC/EF5 CIs: active-level cluster bootstrap (decoys clustered by matched active). "
        "EF1 only for n_act≥150; box-clipped ligands excluded (Gyrase 17, LepB 21).\n"
        "Source: figures/trackA_enrichment_full.json "
        f"(generated_at {payload.get('generated_at', '?')}; "
        f"panel_scope: {payload.get('panel_scope', {}).get('n_panel_targets', '?')} panel / "
        f"{payload.get('panel_scope', {}).get('n_primary_targets', '?')} primary / "
        f"{payload.get('panel_scope', {}).get('n_auxiliary', '?')} auxiliary; "
        f"actives {payload.get('panel_scope', {}).get('n_actives_panel', '?')} panel / "
        f"{payload.get('panel_scope', {}).get('n_actives_primary', '?')} primary)."
    )
    fig.text(0.012, 0.018, footer_note, fontsize=6.5, color="0.35", va="bottom", wrap=True)

    fig.suptitle(
        "Figure 3 — Track A: Vina enrichment across the 9-target panel "
        f"(7 primary + {n_aux} auxiliary)",
        fontsize=12.5, y=0.955,
    )

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    png = FIGURES_DIR / "fig3_trackA_enrichment_full.png"
    pdf = FIGURES_DIR / "fig3_trackA_enrichment_full.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)

    values_payload = [
        {
            "target": entry["target"],
            "kind": entry["_row_kind"],
            "metal": entry.get("metal_dependent"),
            "n_act": entry.get("n_active"),
            "auc": entry.get("auc"),
            "auc_ci": [entry.get("auc_ci_lo"), entry.get("auc_ci_hi")],
            "ef5": entry.get("ef5_percent"),
            "ef5_ci": [entry.get("ef5_ci_lo"), entry.get("ef5_ci_hi")],
            "ef1": entry.get("ef1_percent"),
            "bedroc20": entry.get("bedroc20"),
        }
        for entry in entries
    ]
    report = {
        "figure_png": str(png),
        "figure_pdf": str(pdf),
        "source_json": str(JSON_PATH),
        "source_generated_at": payload.get("generated_at"),
        "primary_caliber": payload.get("primary_caliber"),
        "panel_scope": payload.get("panel_scope"),
        "pa_series_interpretation": payload.get("supplementary_sensitivities", {})
            .get("pa_series_exclusion", {}).get("interpretation"),
        "n_rows": len(entries),
        "n_primary": n_primary,
        "n_auxiliary": n_aux,
        "row_order": [entry["target"] for entry in entries],
        "values": values_payload,
        "notes": [
            "rows sorted by AUROC within primary/auxiliary; auxiliary rows greyed and excluded from aggregates",
            "metal column coloured Zn/Mg/none",
            "primary caliber includes B->C boron surrogates; sensitivity in SI",
            "EF1 only for n_act>=150",
            "LpxC_Pa series exclusion (BM scaffold/ECFP4 cluster) does NOT support a chemotype-concentration explanation for Pa's high AUC; see SI S4a/S10",
        ],
        "validation": {
            "values_match_json": all(
                np.isclose(item["auc"], entry["auc"])
                and np.isclose(item["ef5"], entry["ef5_percent"])
                for item, entry in zip(values_payload, entries)
            ),
            "n_values_checked": len(entries) * 2,
        },
    }
    report_path = FIGURES_DIR / "fig3_trackA_enrichment_full_report.json"
    report_path.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    report["report_path"] = str(report_path)
    return report


if __name__ == "__main__":
    out = make_figure()
    print(json.dumps({k: out[k] for k in
                      ("figure_png", "figure_pdf", "report_path", "source_generated_at", "n_rows")},
                     indent=2, default=str))
