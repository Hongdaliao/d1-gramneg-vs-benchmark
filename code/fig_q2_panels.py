"""Generate manuscript figures 2, 4 and 5 (300 dpi PNG + PDF)."""

from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "figures"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150})


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def load_tracka():
    return json.loads((ROOT / "figures/trackA_enrichment_full.json").read_text(encoding="utf-8"))


def fig2():
    entries = load_tracka()["entries"]
    df = pd.DataFrame(entries)
    df["role"] = np.where(df["auxiliary"].astype(bool), "auxiliary", "primary")
    order = df.sort_values(["role", "n_active"], ascending=[True, True])
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), gridspec_kw={"width_ratios": [1.6, 1.0]})
    colors = {"primary": "#225ea8", "auxiliary": "#969696"}
    y = np.arange(len(order))
    axes[0].barh(y, order["n_active"], color=[colors[r] for r in order["role"]], height=0.65)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(order["target"])
    axes[0].set_xlabel("Curated actives (n)")
    axes[0].set_title("Track A panel composition")
    for yi, (_, row) in zip(y, order.iterrows()):
        axes[0].text(row["n_active"] + 6, yi, f"box {int(row['box_size_A'])} Å; ratio {row['decoy_ratio_observed']:.2f}",
                     va="center", fontsize=6.8, color="0.2")
    axes[0].set_xlim(0, max(order["n_active"]) * 1.55)
    axes[0].legend(handles=[plt.Rectangle((0, 0), 1, 1, color=colors["primary"], label="primary"),
                            plt.Rectangle((0, 0), 1, 1, color=colors["auxiliary"], label="auxiliary")],
                   frameon=False, fontsize=7, loc="lower right")
    # right: regime/pose counts
    axes[1].axis("off")
    lines = [
        "Track A (potency-driven)",
        "  9 targets: 7 primary + 2 auxiliary",
        "  1,748 actives (1,648 primary)",
        "  28,136 successfully docked ligands",
        "  1,710 scored actives in Table 1",
        "",
        "Track B (structure-driven)",
        "  40 inhibitor cognate-redocking poses",
        "  from 5 of 9 curated membrane entries",
        "  AcrB 34 / PBP3 3 / LptB 1 / MexB 1 / MsbA 1",
        "  68 clean poses (7 targets) for sensitivity",
    ]
    axes[1].text(0.02, 0.98, "\n".join(lines), va="top", fontsize=8, family="monospace")
    axes[1].set_title("Data regimes")
    fig.tight_layout()
    save(fig, "fig2_benchmark_design")


def fig4():
    temporal = pd.read_csv(ROOT / "results/trackA/axes/temporal.csv")
    cross = pd.read_csv(ROOT / "results/trackA/axes/cross_species.csv").iloc[0]
    structure = pd.read_csv(ROOT / "results/trackA/axes/structure_source.csv")
    fig, axes = plt.subplots(1, 4, figsize=(12.5, 3.6))

    y = np.arange(len(temporal))
    axes[0].errorbar(temporal["delta_auc"], y,
                     xerr=[temporal["delta_auc"] - temporal["delta_auc_ci_lo"],
                           temporal["delta_auc_ci_hi"] - temporal["delta_auc"]],
                     fmt="o", color="#225ea8", ecolor="0.55", capsize=2, ms=3.5)
    axes[0].axvline(0, color="0.5", ls="--", lw=0.7)
    axes[0].set_yticks(y); axes[0].set_yticklabels(temporal["target"], fontsize=7)
    axes[0].set_xlabel("Delta AUC\n(novel − full)")
    axes[0].set_title("Temporal axis")

    axes[1].bar([0, 1], [cross["auc_ec"], cross["auc_pa"]], color=["#9ecae1", "#3182bd"], width=0.55)
    axes[1].errorbar([0, 1], [cross["auc_ec"], cross["auc_pa"]],
                     yerr=[[cross["auc_ec"] - cross["auc_ec_ci_lo"], cross["auc_pa"] - cross["auc_pa_ci_lo"]],
                           [cross["auc_ec_ci_hi"] - cross["auc_ec"], cross["auc_pa_ci_hi"] - cross["auc_pa"]]],
                     fmt="none", ecolor="0.3", capsize=3)
    axes[1].set_ylim(0, 1); axes[1].set_xticks([0, 1], ["LpxC_Ec", "LpxC_Pa"])
    axes[1].set_ylabel("AUC (951 shared decoys)")
    axes[1].set_title("Cross-species axis")

    for ax, holo, afdb, ylabel, title in (
        (axes[2], "holo_rmsd", "afdb_rmsd", "Ligand RMSD (Å)", "Structure source: RMSD"),
        (axes[3], "holo_plif", "afdb_plif", "PLIF recovery", "Structure source: PLIF"),
    ):
        x = np.arange(len(structure))
        ax.plot(x, structure[holo], color="#9ecae1", lw=0.7, marker="o", ms=2.5, label="holo")
        ax.plot(x, structure[afdb], color="#3182bd", lw=0.7, marker="o", ms=2.5, label="AFDB")
        ax.set_xlabel("Paired pose index")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=6.5)
    fig.tight_layout()
    save(fig, "fig4_three_axis_transfer")


def fig5():
    entries = pd.DataFrame(load_tracka()["entries"])
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.7))

    metal = entries[entries["metal_dependent"].isin(["Zn", "Mg"]) & (entries["panel_role"] == "main")]
    nonmetal = entries[(entries["metal_dependent"] == "none") & (entries["panel_role"] == "main")]
    axes[0].bar([0], [metal["auc"].mean()], color="#de2d26", width=0.45, label=f"metal n={len(metal)}")
    axes[0].bar([1], [nonmetal["auc"].mean()], color="#31a354", width=0.45, label=f"non-metal n={len(nonmetal)}")
    for i, (_, row) in enumerate(metal.iterrows()):
        axes[0].scatter(0 + np.linspace(-0.12, 0.12, len(metal))[i], row["auc"], color="k", s=8, zorder=3)
    for i, (_, row) in enumerate(nonmetal.iterrows()):
        axes[0].scatter(1 + np.linspace(-0.12, 0.12, len(nonmetal))[i], row["auc"], color="k", s=8, zorder=3)
    axes[0].set_xticks([0, 1], ["metal\nprimary", "non-metal\nprimary"])
    axes[0].set_ylim(0, 1); axes[0].set_ylabel("Mean AUC")
    axes[0].set_title("Metal stratification")
    axes[0].legend(frameon=False, fontsize=7)

    boron = pd.read_csv(ROOT / "manuscript/si/table_S5_boron_surrogate.csv")
    x = np.arange(len(boron))
    axes[1].bar(x - 0.18, boron["auc_main"], width=0.36, color="#3182bd", label="with B→C surrogate")
    axes[1].bar(x + 0.18, boron["auc_sensitivity"], width=0.36, color="#9ecae1", label="exclude surrogate")
    axes[1].set_xticks(x, boron["target"]); axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("AUC"); axes[1].set_title("Boron rescue sensitivity")
    axes[1].legend(frameon=False, fontsize=6.5)

    base = pd.read_csv(ROOT / "results/trackA/ecfp4_baseline/baseline_metrics.csv")
    merged = entries[["target", "auc"]].merge(base[["target", "auc"]], on="target", suffixes=("_vina", "_ecfp4"))
    axes[2].scatter(merged["auc_vina"], merged["auc_ecfp4"], s=24, color="#6a51a3")
    for _, row in merged.iterrows():
        axes[2].annotate(row["target"], (row["auc_vina"], row["auc_ecfp4"]), fontsize=6, xytext=(2, 2), textcoords="offset points")
    axes[2].plot([0, 1], [0, 1], ls="--", color="0.6", lw=0.7)
    axes[2].set_xlim(0, 1); axes[2].set_ylim(0, 1)
    axes[2].set_xlabel("Vina AUC"); axes[2].set_ylabel("ECFP4 1-NN AUC")
    axes[2].set_title("Benchmark validity")
    fig.tight_layout()
    save(fig, "fig5_failure_modes")


if __name__ == "__main__":
    fig2(); fig4(); fig5()
    for name in ("fig2_benchmark_design", "fig4_three_axis_transfer", "fig5_failure_modes"):
        for suffix in (".png", ".pdf"):
            p = FIG / f"{name}{suffix}"
            print(p, p.stat().st_size)
