"""Figure 3 draft — Track A enrichment: pilot heatmap + metal-control inset.

Panels:
(a) heatmap of available pilot entries (LpxC_Ec Zn/noZn, DHFR_Ec) × metrics
    (AUROC / EF5% / BEDROC20); pending full-panel rows are hatched;
(b) AUROC bars with ligand-bootstrap 95% CI, metal stratum colour-coded;
(c) inset: the committed metal-control ROC figure (pilot v2, three groups).

Everything is labelled **preliminary** (pilot data only; Track A full panel
pending). Uses the framework output ``figures/trackA_enrichment_pilot_v2.json``.

Usage (from the project ``code`` directory)::

    python -m maes.fig3
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
FIGURES_DIR = PROJECT_ROOT / "figures"

PENDING_TARGETS = [
    ("GyrA_Ec", "Mg"),
    ("GyrB_Ec", "Mg"),
    ("FabI_Ec", "none"),
    ("LepB_Ec", "none"),
    ("TEM1_Ec", "none"),
    ("LpxC_Pa", "Zn"),
]


def _load_framework(tag: str = "pilot_v2") -> Dict[str, Any]:
    path = FIGURES_DIR / f"trackA_enrichment_{tag}.json"
    if not path.exists():
        raise FileNotFoundError(f"run the framework first: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def make_figure(tag: str = "pilot_v2") -> Dict[str, Any]:
    payload = _load_framework(tag)
    entries = payload["entries"]
    labels = [
        f"{entry['target']} ({entry['condition']})" if entry["condition"] != "native" else entry["target"]
        for entry in entries
    ]
    metric_names = ["auc", "ef5_percent", "bedroc20"]
    metric_labels = ["AUROC", "EF5%", "BEDROC20"]
    matrix = np.array([[float(entry[name]) for name in metric_names] for entry in entries], dtype=float)

    fig = plt.figure(figsize=(16.2, 6.2))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.05, 1.25], wspace=0.32)

    # ---------------------------------------------------------------- (a) heatmap
    ax = fig.add_subplot(grid[0, 0])
    n_rows_available = len(entries)
    n_pending = len(PENDING_TARGETS)
    data = np.full((n_rows_available + n_pending, len(metric_names)), np.nan)
    data[:n_rows_available, :] = matrix
    im = ax.imshow(data, cmap="RdYlBu_r", vmin=0.0, vmax=1.0, aspect="auto")
    for i in range(n_rows_available):
        for j in range(len(metric_names)):
            value = data[i, j]
            ax.text(j, i, f"{value:.2f}" if np.isfinite(value) else "—", ha="center", va="center",
                    fontsize=8.5, color="black" if np.isfinite(value) else "0.4",
                    fontweight="bold" if j == 0 else "normal")
    for offset, (target, metal) in enumerate(PENDING_TARGETS):
        row = n_rows_available + offset
        ax.add_patch(plt.Rectangle((-0.5, row - 0.5), len(metric_names), 1.0,
                                   facecolor="white", edgecolor="0.85", hatch="///", zorder=2, alpha=0.65))
        ax.text(len(metric_names) - 0.55, row, "pending", ha="right", va="center", fontsize=7, color="0.45", zorder=3)
    row_labels = labels + [f"{target} [{metal}]" for target, metal in PENDING_TARGETS]
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_xticks(range(len(metric_names)))
    ax.set_xticklabels(metric_labels, fontsize=9)
    ax.set_title("(a) Pilot entries × metrics (full Track A panel pending)", fontsize=10)
    for i, entry in enumerate(entries):
        if entry["metal_stratum"] == "metal":
            ax.get_yticklabels()[i].set_color("#C0392B")
        else:
            ax.get_yticklabels()[i].set_color("#1F618D")
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label("metric value", fontsize=8)
    ax.set_xlabel("EF1% suppressed: n_active < 150 (R1)", fontsize=7.5)

    # ---------------------------------------------------------------- (b) AUROC bars
    ax2 = fig.add_subplot(grid[0, 1])
    x = np.arange(len(entries))
    aucs = np.array([entry["auc"] for entry in entries], dtype=float)
    lo = np.array([entry.get("auc_ci_lo", np.nan) for entry in entries], dtype=float)
    hi = np.array([entry.get("auc_ci_hi", np.nan) for entry in entries], dtype=float)
    colors = ["#C0392B" if entry["metal_stratum"] == "metal" else "#1F618D" for entry in entries]
    err = np.vstack([np.maximum(aucs - lo, 0), np.maximum(hi - aucs, 0)])
    ax2.bar(x, aucs, width=0.55, color=colors, alpha=0.85, yerr=err, capsize=4,
            error_kw={"ecolor": "0.3", "lw": 1.0})
    ax2.axhline(0.5, color="0.35", ls="--", lw=1)
    ax2.text(len(entries) - 0.45, 0.52, "chance", fontsize=7, color="0.35", ha="right")
    for xi, value in zip(x, aucs):
        ax2.text(xi, min(value + 0.2, 0.95), f"{value:.2f}", ha="center", fontsize=8)
    ax2.set_xticks(x)
    ax2.set_xticklabels([label.replace(" (", "\n(") for label in labels], fontsize=7.5)
    ax2.set_ylim(0.0, 1.0)
    ax2.set_ylabel("AUROC (95% bootstrap CI)", fontsize=9)
    ax2.set_title("(b) AUROC by metal stratum", fontsize=10)
    ax2.legend(handles=[Patch(color="#C0392B", label="metal"), Patch(color="#1F618D", label="non-metal")],
               fontsize=7.5, loc="upper left")

    # ---------------------------------------------------------------- (c) inset image
    image_path = FIGURES_DIR / "pilot_v2_zn_nozn_dhfr.png"
    ax3 = fig.add_subplot(grid[0, 2])
    if image_path.exists():
        img = plt.imread(image_path)
        ax3.imshow(img)
        ax3.set_title("(c) Known metal-coordination limitation (sanity control)", fontsize=10)
    else:
        ax3.text(0.5, 0.5, "metal-control figure pending", ha="center", va="center")
        ax3.set_title("(c) Known metal-coordination limitation (sanity control)", fontsize=10)
    ax3.axis("off")

    fig.suptitle(
        "Figure 3 (draft) — Track A enrichment, preliminary pilot evidence "
        f"(n = {payload['n_files']} tables; full Track A panel pending)",
        fontsize=11,
    )
    fig.text(
        0.99, 0.005,
        "preliminary — pilot data only; EF1% suppressed where n_active < 150 (R1); pending rows hatched",
        ha="right", fontsize=6.5, color="0.45",
    )
    fig.subplots_adjust(left=0.10, right=0.985, top=0.86, bottom=0.10, wspace=0.32)

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    png = FIGURES_DIR / "fig3_trackA_enrichment_draft.png"
    pdf = FIGURES_DIR / "fig3_trackA_enrichment_draft.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    report = {
        "figure_png": str(png),
        "figure_pdf": str(pdf),
        "tag": tag,
        "status": "draft_preliminary",
        "n_entries": len(entries),
        "n_pending_targets": len(PENDING_TARGETS),
        "framework_json": str(FIGURES_DIR / f"trackA_enrichment_{tag}.json"),
        "entries": entries,
        "notes": [
            "heatmap rows: available pilot entries + hatched pending Track A targets",
            "EF1% suppressed (n_active < 150); EF5%/BEDROC/AUROC reported",
            "metal stratum colour: LpxC=Zn (red), DHFR=none (blue)",
            "inset: KNOWN metal-coordination limitation (sanity control), pilot_v2 three-group ROC; not a new finding",
        ],
    }
    report_path = FIGURES_DIR / "fig3_trackA_enrichment_draft_report.json"
    report_path.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    report["report_path"] = str(report_path)
    return report


if __name__ == "__main__":
    out = make_figure()
    print(json.dumps({k: out[k] for k in ("figure_png", "report_path", "n_entries")}, indent=2, default=str))
