"""Generate v3 figures Fig2/4/5/6 (300 dpi PNG+PDF)."""

from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "figures/v3"
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150})


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def fig2():
    tracka = json.loads((ROOT / "figures/trackA_enrichment_full.json").read_text(encoding="utf-8"))
    df = pd.DataFrame(tracka["entries"])
    df["role"] = np.where(df["auxiliary"].astype(bool), "auxiliary", "primary")
    order = df.sort_values(["role", "n_active"], ascending=[True, True])
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), gridspec_kw={"width_ratios": [1.6, 1.0]})
    colors = {"primary": "#225ea8", "auxiliary": "#969696"}
    y = np.arange(len(order))
    axes[0].barh(y, order["n_active"], color=[colors[r] for r in order["role"]], height=0.65)
    axes[0].set_yticks(y); axes[0].set_yticklabels(order["target"])
    axes[0].set_xlabel("Curated actives (n)"); axes[0].set_title("Track A panel composition")
    for yi, (_, row) in zip(y, order.iterrows()):
        axes[0].text(row["n_active"] + 6, yi, f"box {int(row['box_size_A'])} Å; ratio {row['decoy_ratio_observed']:.2f}",
                     va="center", fontsize=6.8, color="0.2")
    axes[0].set_xlim(0, max(order["n_active"]) * 1.55)
    axes[0].legend(handles=[plt.Rectangle((0, 0), 1, 1, color=colors["primary"], label="primary"),
                            plt.Rectangle((0, 0), 1, 1, color=colors["auxiliary"], label="auxiliary")],
                   frameon=False, fontsize=7, loc="lower right")
    axes[1].axis("off")
    lines = [
        "Track A (potency-driven)",
        "  9 targets: 7 primary + 2 auxiliary",
        "  1,748 curated actives (1,648 primary)",
        "  28,136 successfully docked ligands",
        "  1,710 scored actives in Table 1",
        "",
        "Track B (co-crystal poses)",
        "  40 poses from 5 of 9 membrane entries",
        "  AcrB 34 / PBP3 3 / LptB 1 / MexB 1 / MsbA 1",
        "  classes: inhibitor 4; antibiotic-substrate 24",
        "           covalent 3; fragment 8; unresolved 1",
        "  AcrB: 34 inhibitor-class / 55 non-polymer poses",
    ]
    axes[1].text(0.02, 0.98, "\n".join(lines), va="top", fontsize=7.5, family="monospace")
    axes[1].set_title("Data regimes and ligand classes")
    fig.tight_layout(); save(fig, "fig2_benchmark_design")


def fig4():
    temporal = pd.read_csv(ROOT / "results/trackA/axes/temporal.csv")
    cross = json.loads((ROOT / "results/trackA/axes/cross_species_bootstrap_corrected.json").read_text(encoding="utf-8"))
    structure = pd.read_csv(ROOT / "results/trackA/axes/structure_source_qcpassed.csv")
    fig, axes = plt.subplots(1, 4, figsize=(13.0, 3.7))
    y = np.arange(len(temporal))
    axes[0].errorbar(temporal["delta_auc"], y,
                     xerr=[temporal["delta_auc"] - temporal["delta_auc_ci_lo"], temporal["delta_auc_ci_hi"] - temporal["delta_auc"]],
                     fmt="o", color="#225ea8", ecolor="0.55", capsize=2, ms=3.5)
    axes[0].axvline(0, color="0.5", ls="--", lw=0.7)
    axes[0].set_yticks(y); axes[0].set_yticklabels(temporal["target"], fontsize=7)
    axes[0].set_xlabel("Delta AUC (novel − full)"); axes[0].set_title("Temporal (descriptive)")
    # cross-species with corrected CI
    p = cross["point"]; cilo, cihi = cross["reviewer_authoritative"]["decoy_level_delta_ci95"]
    axes[1].bar([0], [p["auc_ec"]], color="#9ecae1", width=0.5)
    axes[1].bar([1], [p["auc_pa"]], color="#3182bd", width=0.5)
    axes[1].set_xticks([0, 1], ["LpxC_Ec", "LpxC_Pa"]); axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("AUC (951 shared decoys)")
    axes[1].set_title(f"Paired comparison\nDelta = {p['delta_auc_pa_minus_ec']:.3f} (CI {cilo:.3f}–{cihi:.3f})", fontsize=7)
    # structure source QC-passed only
    x = np.arange(len(structure))
    axes[2].plot(x, structure["holo_rmsd"], color="#9ecae1", lw=0.8, marker="o", ms=2.5, label="holo")
    axes[2].plot(x, structure["afdb_rmsd"], color="#3182bd", lw=0.8, marker="o", ms=2.5, label="AFDB")
    axes[2].set_xlabel("QC-passed pair index"); axes[2].set_ylabel("Ligand RMSD (Å)")
    axes[2].set_title("Structure source: RMSD (n=8)")
    axes[2].legend(frameon=False, fontsize=6.5)
    axes[3].plot(x, structure["holo_plif"], color="#9ecae1", lw=0.8, marker="o", ms=2.5, label="holo")
    axes[3].plot(x, structure["afdb_plif"], color="#3182bd", lw=0.8, marker="o", ms=2.5, label="AFDB")
    axes[3].set_xlabel("QC-passed pair index"); axes[3].set_ylabel("PLIF recovery")
    axes[3].set_title("Structure source: PLIF (n=8)")
    axes[3].legend(frameon=False, fontsize=6.5)
    fig.tight_layout(); save(fig, "fig4_three_axis_transfer")


def fig5():
    lig = pd.read_csv(ROOT / "data/trackB/ligand_class_review.csv")
    counts = lig["suggested_class"].value_counts().reindex(["inhibitor", "antibiotic-substrate", "covalent", "fragment", "substrate"]).fillna(0)
    base = pd.read_csv(ROOT / "results/trackA/ecfp4_baseline/baseline_metrics_test_only.csv")
    vina = pd.DataFrame(json.loads((ROOT / "figures/trackA_enrichment_full.json").read_text(encoding="utf-8"))["entries"])[["target", "auc"]]
    merged = vina.merge(base[["target", "auc_test_only"]], on="target", how="left")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8))
    labels = {"inhibitor": "inhibitor", "antibiotic-substrate": "antibiotic\nsubstrate", "covalent": "covalent\nβ-lactam", "fragment": "fragment", "substrate": "unresolved"}
    axes[0].bar([labels[c] for c in counts.index], counts.values,
                color=["#225ea8", "#9ecae1", "#fdae6b", "#bdbdbd", "#fdd0a2"])
    axes[0].set_ylabel("Track B poses"); axes[0].set_title("Ligand-class composition (40 poses)")
    for i, v in enumerate(counts.values): axes[0].text(i, v + 0.4, str(int(v)), ha="center", fontsize=8)
    axes[1].scatter(merged["auc"], merged["auc_test_only"], s=28, color="#6a51a3")
    for _, row in merged.iterrows():
        if pd.notna(row["auc_test_only"]):
            axes[1].annotate(row["target"], (row["auc"], row["auc_test_only"]), fontsize=6, xytext=(2, 2), textcoords="offset points")
    axes[1].plot([0, 1], [0, 1], ls="--", color="0.6", lw=0.7)
    axes[1].set_xlim(0, 1); axes[1].set_ylim(0, 1)
    axes[1].set_xlabel("Vina AUC"); axes[1].set_ylabel("ECFP4 1-NN AUC (test-only)")
    axes[1].set_title("2D baseline vs Vina (test-only)")
    fig.tight_layout(); save(fig, "fig5_ligand_class_2d")


def fig6():
    report = json.loads((ROOT / "results/maes/w12/w12_report_formal.json").read_text(encoding="utf-8"))
    cc = json.loads((ROOT / "results/maes/w12/w12_complete_case.json").read_text(encoding="utf-8"))
    bins = report["primary"]["reliability_bins"]
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.7))
    axes[0].plot([0, 1], [0, 1], ls="--", color="0.6")
    if bins:
        axes[0].plot([b["mean_pred"] for b in bins], [b["observed"] for b in bins], "o-", color="#225ea8")
    axes[0].set_xlim(0, 1); axes[0].set_ylim(0, 1)
    axes[0].set_xlabel("Predicted MAES (primary V+M)"); axes[0].set_ylabel("Observed acceptable rate")
    axes[0].set_title("Calibration (full sample, n=146)")
    axes[1].bar([0, 1], [report["primary"]["metrics"]["auroc_primary"], 0.4956], color=["#225ea8", "#9ecae1"], width=0.55)
    axes[1].set_xticks([0, 1], ["full\nn=146", "exploratory\nV+M+S"]); axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("AUROC"); axes[1].set_title("AUROC")
    # complete-case bars: show AUROC/Brier
    labels = ["AUROC", "Brier"]
    full = [report["primary"]["metrics"]["auroc_primary"], report["primary"]["metrics"]["brier_primary"]]
    comp = [cc["reviewer_authoritative"]["auroc"], cc["reviewer_authoritative"]["brier"]]
    x = np.arange(2)
    axes[2].bar(x - 0.18, full, width=0.36, color="#9ecae1", label="full n=146")
    axes[2].bar(x + 0.18, comp, width=0.36, color="#3182bd", label="complete n=132")
    axes[2].set_xticks(x, labels); axes[2].set_title("Complete-case sensitivity")
    axes[2].legend(frameon=False, fontsize=6.5)
    fig.tight_layout(); save(fig, "fig6_w12_calibration_formal")


if __name__ == "__main__":
    fig2(); fig4(); fig5(); fig6()
    for name in ("fig2_benchmark_design", "fig4_three_axis_transfer", "fig5_ligand_class_2d", "fig6_w12_calibration_formal"):
        for suffix in (".png", ".pdf"):
            p = FIG / f"{name}{suffix}"; print(p, p.stat().st_size)
