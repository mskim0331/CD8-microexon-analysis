#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.stats import pearsonr, spearmanr


def parse_args():
    p = argparse.ArgumentParser(
        description="Create presentation-ready Whippet validation plots."
    )
    p.add_argument(
        "--concordance",
        default="analysis_results/06_whippet_candidate_concordance.tsv",
    )
    p.add_argument(
        "--summary",
        default="analysis_results/06_whippet_concordance_summary.tsv",
    )
    p.add_argument("--outdir", default="analysis_results")
    return p.parse_args()


def as_bool(series):
    return series.astype(str).str.lower().eq("true")


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.concordance, sep="\t")
    summary = pd.read_csv(args.summary, sep="\t")

    mapped = df.loc[as_bool(df["whippet_mapped"])].copy()
    x = pd.to_numeric(
        mapped["delta_PSI_activation_minus_ctrl"], errors="coerce"
    )
    y = pd.to_numeric(
        mapped["whippet_delta_activation_minus_ctrl"], errors="coerce"
    )
    ok = x.notna() & y.notna()
    mapped = mapped.loc[ok].copy()
    x = x.loc[ok].to_numpy(dtype=float)
    y = y.loc[ok].to_numpy(dtype=float)

    concordant = as_bool(mapped["direction_concordant"]).to_numpy()
    n = len(mapped)
    n_concordant = int(concordant.sum())
    concordance_pct = 100.0 * n_concordant / n if n else np.nan

    r = float(pearsonr(x, y).statistic)
    rho = float(spearmanr(x, y).statistic)

    lim = max(abs(x).max(), abs(y).max()) if n else 1.0
    lim = max(lim, 0.25) * 1.08

    # Scatter plot: MicroExonator vs Whippet DeltaPSI.
    fig, ax = plt.subplots(figsize=(7.8, 7.0))

    if (~concordant).any():
        ax.scatter(
            x[~concordant],
            y[~concordant],
            s=34,
            alpha=0.65,
            marker="x",
            label=f"Discordant direction (n={(~concordant).sum()})",
        )

    ax.scatter(
        x[concordant],
        y[concordant],
        s=34,
        alpha=0.65,
        label=f"Concordant direction (n={n_concordant})",
    )

    ax.plot(
        [-lim, lim],
        [-lim, lim],
        linestyle="--",
        linewidth=1.2,
        color="black",
        label="Perfect agreement (y = x)",
    )
    ax.axhline(0, linewidth=0.8, color="0.6")
    ax.axvline(0, linewidth=0.8, color="0.6")

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("MicroExonator ΔPSI\n(Activation − Control)")
    ax.set_ylabel("Whippet ΔPSI\n(Activation − Control)")
    ax.set_title("Independent Whippet Validation")

    stats_text = (
        f"Matched responses: {n}\n"
        f"Direction concordance: {n_concordant}/{n} ({concordance_pct:.1f}%)\n"
        f"Pearson r = {r:.3f}\n"
        f"Spearman ρ = {rho:.3f}"
    )
    ax.text(
        0.04,
        0.96,
        stats_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.9, edgecolor="0.75"),
    )
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(
        outdir / "06_whippet_validation_scatter.png",
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        outdir / "06_whippet_validation_scatter.svg",
        bbox_inches="tight",
    )
    plt.close(fig)

    # Optional presentation companion: concordance by time point.
    tp = summary.loc[summary["timepoint"] != "ALL"].copy()
    labels = tp["timepoint"].replace({
        "30min": "30 min",
        "3hr": "3 h",
        "12hr": "12 h",
        "24hr": "24 h",
        "48hr": "48 h",
        "72hr": "72 h",
        "7day": "7 d",
    }).tolist()
    fractions = pd.to_numeric(
        tp["direction_concordance_fraction_among_mapped"], errors="coerce"
    ).to_numpy(dtype=float) * 100.0
    mapped_n = pd.to_numeric(
        tp["mapped_to_whippet"], errors="coerce"
    ).to_numpy(dtype=int)
    concordant_n = pd.to_numeric(
        tp["direction_concordant"], errors="coerce"
    ).to_numpy(dtype=int)

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    bars = ax.bar(np.arange(len(labels)), fractions)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 110)
    ax.set_ylabel("Direction concordance among mapped responses (%)")
    ax.set_xlabel("Activation time point")
    ax.set_title("Whippet Direction Concordance by Time Point")

    for bar, pct, c, m in zip(bars, fractions, concordant_n, mapped_n):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 2,
            f"{pct:.1f}%\n({c}/{m})",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    ax.axhline(
        86.4,
        linestyle="--",
        linewidth=1.0,
        color="black",
        label="Overall = 86.4%",
    )
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(
        outdir / "06_whippet_direction_concordance_by_timepoint.png",
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        outdir / "06_whippet_direction_concordance_by_timepoint.svg",
        bbox_inches="tight",
    )
    plt.close(fig)

    report = [
        "WHIPPET VALIDATION FIGURES",
        "",
        f"Mapped responses plotted: {n}",
        f"Direction concordant: {n_concordant}/{n} ({concordance_pct:.1f}%)",
        f"Pearson r: {r:.6f}",
        f"Spearman rho: {rho:.6f}",
        "",
        "Scatter axes use the same sign convention:",
        "  DeltaPSI = Activation - Control",
    ]
    (outdir / "06_whippet_validation_plot_report.txt").write_text(
        "\n".join(report) + "\n",
        encoding="utf-8",
    )

    print("\n".join(report))


if __name__ == "__main__":
    main()
