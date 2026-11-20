"""Figure 6 draft — MAES calibration on the first real Track B data.

Panels:
(a) reliability diagram (evaluation label RMSD < 1 Å; calibration label shown dashed)
(b) MAES score vs RMSD scatter with the mandatory decision lines
(c) ablation bars (AUROC left axis, Brier right axis)

Everything is labelled as **preliminary first-pass** (n = 36 analysis poses,
8 calibration positives, single method Vina, small-n, AcrB-dominated); the
figure must not be used to claim final performance.

Usage (from the project ``code`` directory)::

    python -m maes.fig6
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "figure.dpi": 120,
    }
)
import numpy as np
import pandas as pd

from .api import build_features_from_manifest, load_pose_manifest
from .calibration import MAESModel, build_calibration_labels, build_evaluation_labels
from .config import default_config
from .evaluation import auroc, brier_score, ece, reliability_curve, spearman_rho
from .ingest import DEFAULT_MANIFEST, DEFAULT_PREDICTIONS, FIGURES_DIR, normalise_predictions

ABLATION_ORDER = ("full", "no_M", "no_S", "no_strain", "uncalibrated")


def make_figure(
    predictions_path: Optional[str | Path] = None,
    manifest_path: Optional[str | Path] = None,
    tag: str = "draft",
) -> Dict[str, Any]:
    cfg = default_config()
    preds = normalise_predictions(pd.read_csv(predictions_path or DEFAULT_PREDICTIONS))
    manifest = load_pose_manifest(manifest_path or DEFAULT_MANIFEST)
    analysis = preds[preds["rmsd"].notna() & preds["posebusters_pass"].notna()].copy()
    features = build_features_from_manifest(manifest, analysis, subset="primary_main")

    model = MAESModel(cfg).fit(features)
    maes = model.predict(features)
    y_cal = build_calibration_labels(features, cfg.labels)
    y_eval = build_evaluation_labels(features, cfg.labels, which="primary")
    rmsd = pd.to_numeric(features["rmsd"], errors="coerce").to_numpy(dtype=float)

    metrics = {
        "n_analysis": int(len(features)),
        "n_calibration_positive": int(y_cal.sum()),
        "n_evaluation_positive": int(y_eval.sum()),
        "auroc_primary": float(auroc(y_eval, maes)),
        "brier_primary": float(brier_score(y_eval, maes)),
        "ece_primary": float(ece(y_eval, maes, n_bins=10)),
        "rho_maes_rmsd": float(spearman_rho(maes, rmsd)),
    }
    rel_eval = reliability_curve(y_eval, maes, n_bins=4)
    rel_cal = reliability_curve(y_cal, maes, n_bins=4)

    # ablation numbers: prefer the committed ingest report, fall back to nothing
    report_path = FIGURES_DIR / "ingest_report_holo.json"
    ablation = []
    if report_path.exists():
        payload = json.loads(report_path.read_text(encoding="utf-8"))
        ablation = [row for row in payload.get("ablation", []) if row.get("variant") in ABLATION_ORDER]
    if not ablation:  # pragma: no cover - defensive
        from .evaluation import ablation_study

        ablation = ablation_study(features, cfg).to_dict(orient="records")

    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.5), gridspec_kw={"width_ratios": [1.0, 1.15, 1.3]})

    # ------------------------------------------------------------------ (a) reliability
    ax = axes[0]
    ax.plot([0, 1], [0, 1], ls="--", color="0.55", lw=1, label="perfect calibration")
    valid_eval = rel_eval.dropna(subset=["mean_pred", "frac_pos"])
    valid_cal = rel_cal.dropna(subset=["mean_pred", "frac_pos"])
    sizes_eval = 18 + 40 * (rel_eval["count"].to_numpy(dtype=float) / max(1.0, rel_eval["count"].max()))
    ax.plot(
        valid_eval["mean_pred"], valid_eval["frac_pos"], "-", color="#0072B2", lw=1.6,
        marker="o", ms=5, label="evaluation label (RMSD < 1 Å)",
    )
    ax.scatter(
        valid_eval["mean_pred"], valid_eval["frac_pos"], s=sizes_eval[rel_eval["count"] > 0],
        color="#0072B2", alpha=0.25, zorder=1,
    )
    ax.plot(
        valid_cal["mean_pred"], valid_cal["frac_pos"], ":", color="#E69F00", lw=1.3,
        marker="s", ms=4, label="calibration label (RMSD < 2 Å ∧ PB)",
    )
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Mean predicted P(acceptable)")
    ax.set_ylabel("Observed fraction")
    ax.set_title("(a) Reliability (4 bins; evaluation label)")
    ax.legend(fontsize=7, loc="upper left")
    ax.text(
        0.98, 0.03,
        f"ECE (10-bin) = {metrics['ece_primary']:.3f}\nBrier = {metrics['brier_primary']:.3f}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.8", alpha=0.85),
    )

    # ------------------------------------------------------------------ (b) MAES vs RMSD
    ax = axes[1]
    acceptable = y_cal.astype(bool)
    ax.scatter(rmsd[~acceptable], maes[~acceptable], s=26, c="#8A8A8A", alpha=0.75,
               edgecolors="white", linewidths=0.4, label="not acceptable (calibration label)")
    ax.scatter(rmsd[acceptable], maes[acceptable], s=34, c="#009E73", alpha=0.9,
               edgecolors="white", linewidths=0.5, label="acceptable (calibration label)")
    ax.axvline(1.0, color="#D55E00", ls=":", lw=1.1)
    ax.axvline(2.0, color="#D55E00", ls="--", lw=1.1)
    ax.axhline(0.5, color="0.45", ls="--", lw=0.9)
    ax.text(1.02, 0.03, "RMSD = 1 Å", fontsize=6.5, color="#D55E00", rotation=90)
    ax.text(2.04, 0.03, "RMSD = 2 Å", fontsize=6.5, color="#D55E00", rotation=90)
    ax.text(8.9, 0.515, "MAES = 0.5", fontsize=6.5, color="0.35", ha="right")
    rho = metrics["rho_maes_rmsd"]
    ax.set_xlabel("Reference-pose RMSD (Å)")
    ax.set_ylabel("MAES P(acceptable)")
    ax.set_title(f"(b) MAES vs RMSD (Spearman ρ = {rho:.2f})", pad=10)
    ax.legend(fontsize=7, loc="lower right")
    ax.text(
        0.02, 0.97,
        f"n = {metrics['n_analysis']} ({metrics['n_calibration_positive']} acceptable)",
        transform=ax.transAxes, ha="left", va="top", fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.8", alpha=0.85),
    )

    # ------------------------------------------------------------------ (c) ablation bars
    ax = axes[2]
    rows = {str(row["variant"]): row for row in ablation}
    variants = [v for v in ABLATION_ORDER if v in rows]
    x = np.arange(len(variants))
    aurocs = [float(rows[v]["auroc"]) for v in variants]
    briers = [float(rows[v]["brier"]) for v in variants]
    bars = ax.bar(x, aurocs, width=0.55, color="#0072B2", alpha=0.85, label="AUROC (left)")
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("AUROC", color="#0072B2")
    ax.tick_params(axis="y", labelcolor="#0072B2")
    ax.set_xticks(x, variants, rotation=30, ha="right", fontsize=8)
    ax.axhline(0.65, color="#0072B2", ls=":", lw=1)
    ax.yaxis.grid(True, alpha=0.25)
    ax.set_axisbelow(True)

    ax2 = ax.twinx()
    ax2.plot(x, briers, "o--", color="#D55E00", lw=1.4, ms=5, label="Brier (right)")
    ax2.set_ylim(0.0, max(1.0, max(briers) * 1.1))
    ax2.set_ylabel("Brier score", color="#D55E00")
    ax2.tick_params(axis="y", labelcolor="#D55E00")
    for xi, val in zip(x, briers):
        ax2.annotate(f"{val:.3f}", (xi, val), textcoords="offset points", xytext=(10, 4),
                     ha="left", fontsize=6.5, color="#D55E00")
    ax.set_title("(c) Ablation (dotted: AUROC = 0.65 gate; uncalibrated Brier is the calibration argument)", fontsize=8.5)
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, fontsize=6.5, loc="center left")

    fig.suptitle(
        "Figure 6 (draft) — MAES calibration, preliminary first-pass "
        f"(n = {metrics['n_analysis']}, {metrics['n_calibration_positive']} calibration positives; "
        "single method Vina; small-n; strict physics PB)",
        fontsize=10.5,
    )
    fig.text(
        0.99, 0.005,
        "preliminary — do not quote in Abstract/Conclusions; revise under the W12 extension gate",
        ha="right", fontsize=6.5, color="0.45",
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.93))

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    png = FIGURES_DIR / f"fig6_maes_calibration_{tag}.png"
    pdf = FIGURES_DIR / f"fig6_maes_calibration_{tag}.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)

    report = {
        "figure_png": str(png),
        "figure_pdf": str(pdf),
        "tag": tag,
        "status": "preliminary_first_pass",
        "metrics": metrics,
        "reliability_evaluation": rel_eval.to_dict(orient="records"),
        "reliability_calibration": rel_cal.to_dict(orient="records"),
        "ablation": ablation,
        "caveats": [
            "n = 36 analysis poses, 8 calibration positives — small-n",
            "single method (Vina 1.2.5); no Boltz-2/DrugCLIP in this calibration",
            "AcrB-dominated pose set (34/40 inhibitor poses)",
            "reliability diagram uses 4 equal-width bins (small-n; avoids empty bins)",
            "Platt calibration only; isotonic is gated to n_pos >= 100",
            "IRT is descriptive only (36 obs x 15 clusters)",
            "15 clusters only — wild cluster bootstrap reported as robustness; CIs may be narrow",
            "|rho| = %.2f is below the preregistered 0.5-0.8 band — deviation reported, no post-hoc tuning" % abs(rho),
            "do not quote in Abstract/Conclusions; W12 gate: >=100 poses x >=3 methods + Track A before the formal performance claim",
        ],
    }
    report_path = FIGURES_DIR / f"fig6_maes_calibration_{tag}_report.json"
    report_path.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    report["report_path"] = str(report_path)
    return report


if __name__ == "__main__":
    out = make_figure()
    print(json.dumps({k: out[k] for k in ("figure_png", "report_path", "metrics")}, indent=2))
