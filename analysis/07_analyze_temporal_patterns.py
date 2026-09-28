#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch

from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform


TIMEPOINT_ORDER = ["30min", "3hr", "12hr", "24hr", "48hr", "72hr", "7day"]
TIMEPOINT_HOURS = np.array([0.5, 3.0, 12.0, 24.0, 48.0, 72.0, 168.0], dtype=float)


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Step 7: summarize and visualize temporal DeltaPSI patterns among "
            "Whippet-matched MicroExonator candidate microexons."
        )
    )
    p.add_argument(
        "--timecourse",
        default="analysis_results/02_ctrl_vs_timepoint_deltaPSI.tsv",
    )
    p.add_argument(
        "--whippet",
        default="analysis_results/06_whippet_candidate_concordance.tsv",
    )
    p.add_argument(
        "--gene-annotation",
        default="analysis_results/03_microexon_gene_annotation.tsv",
    )
    p.add_argument("--outdir", default="analysis_results")
    return p.parse_args()


def as_bool(series):
    return series.astype(str).str.lower().eq("true")


def normalize_timepoint(value):
    s = str(value).lower().replace("_", "").replace("-", "")
    aliases = {
        "30min": ["30min"],
        "3hr": ["3hr", "3h"],
        "12hr": ["12hr", "12h"],
        "24hr": ["24hr", "24h"],
        "48hr": ["48hr", "48h"],
        "72hr": ["72hr", "72h"],
        "7day": ["7day", "7d"],
    }
    for canonical, vals in aliases.items():
        if any(v in s for v in vals):
            return canonical
    return None


def interpolate_for_clustering(row):
    y = np.asarray(row, dtype=float)
    good = np.isfinite(y)
    if good.sum() == 0:
        return np.zeros_like(y)
    if good.sum() == 1:
        return np.repeat(y[good][0], len(y))
    return np.interp(
        TIMEPOINT_HOURS,
        TIMEPOINT_HOURS[good],
        y[good],
    )


def row_zscore(matrix):
    out = np.zeros_like(matrix, dtype=float)
    for i, row in enumerate(matrix):
        mean = float(np.mean(row))
        sd = float(np.std(row, ddof=0))
        if sd > 1e-12:
            out[i] = (row - mean) / sd
        else:
            out[i] = 0.0
    return out


def silhouette_mean(distance_matrix, labels):
    labels = np.asarray(labels)
    scores = []
    for i in range(len(labels)):
        own = labels == labels[i]
        own[i] = False
        if own.sum() == 0:
            scores.append(0.0)
            continue

        a = float(distance_matrix[i, own].mean())
        b_vals = []
        for lab in sorted(set(labels)):
            if lab == labels[i]:
                continue
            other = labels == lab
            if other.sum():
                b_vals.append(float(distance_matrix[i, other].mean()))

        if not b_vals:
            scores.append(0.0)
            continue

        b = min(b_vals)
        denom = max(a, b)
        scores.append(0.0 if denom == 0 else (b - a) / denom)
    return float(np.mean(scores))


def choose_cluster_count(z):
    if len(z) < 3:
        return np.ones(len(z), dtype=int), 1, np.nan, pd.DataFrame()

    condensed = pdist(z, metric="euclidean")
    tree = linkage(condensed, method="ward")
    dist_matrix = squareform(condensed)

    rows = []
    best = None
    for k in range(2, min(6, len(z) - 1) + 1):
        labels = fcluster(tree, t=k, criterion="maxclust")
        actual_k = len(set(labels))
        if actual_k < 2:
            continue
        score = silhouette_mean(dist_matrix, labels)
        rows.append({
            "requested_k": k,
            "actual_k": actual_k,
            "mean_silhouette": score,
        })
        candidate = (score, -k, labels)
        if best is None or candidate[:2] > best[:2]:
            best = candidate

    score_table = pd.DataFrame(rows)
    if best is None:
        labels = np.ones(len(z), dtype=int)
        return labels, 1, np.nan, score_table

    labels = best[2]
    chosen_k = len(set(labels))
    chosen_score = best[0]
    return labels, chosen_k, chosen_score, score_table


def first_strong_timepoint(sub):
    sub = sub.set_index("timepoint").reindex(TIMEPOINT_ORDER)
    strong = (
        as_bool(sub["valid_comparison"].fillna(False))
        & (pd.to_numeric(sub["abs_delta_PSI"], errors="coerce") >= 0.20)
        & as_bool(sub["complete_group_separation"].fillna(False))
    )
    vals = [tp for tp in TIMEPOINT_ORDER if bool(strong.loc[tp])]
    return vals[0] if vals else ""


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    tc = pd.read_csv(args.timecourse, sep="\t")
    wh = pd.read_csv(args.whippet, sep="\t")
    gene = pd.read_csv(args.gene_annotation, sep="\t", dtype=str)

    tc["timepoint"] = tc["timepoint"].map(normalize_timepoint)
    wh["timepoint"] = wh["timepoint"].map(normalize_timepoint)

    if tc["timepoint"].isna().any() or wh["timepoint"].isna().any():
        raise ValueError("Unrecognized timepoint found in Step 2 or Step 6 table.")

    wh_mapped = as_bool(wh["whippet_mapped"])
    candidate_mes = sorted(wh.loc[wh_mapped, "ME"].dropna().unique().tolist())
    if not candidate_mes:
        raise ValueError("No Whippet-matched candidate microexons found.")

    tc_sub = tc.loc[tc["ME"].isin(candidate_mes)].copy()
    tc_sub["delta_for_temporal_analysis"] = pd.to_numeric(
        tc_sub["delta_PSI_activation_minus_ctrl"], errors="coerce"
    )
    tc_sub.loc[
        ~as_bool(tc_sub["valid_comparison"]),
        "delta_for_temporal_analysis",
    ] = np.nan

    matrix = (
        tc_sub.pivot(
            index="ME",
            columns="timepoint",
            values="delta_for_temporal_analysis",
        )
        .reindex(index=candidate_mes, columns=TIMEPOINT_ORDER)
    )

    matrix.index.name = "ME"
    matrix.to_csv(
        outdir / "07_whippet_matched_microexon_deltaPSI_matrix.tsv",
        sep="\t",
        na_rep="NA",
    )

    clustering_input = np.vstack([
        interpolate_for_clustering(matrix.loc[me].to_numpy(dtype=float))
        for me in candidate_mes
    ])
    z = row_zscore(clustering_input)
    raw_labels, chosen_k, silhouette, score_table = choose_cluster_count(z)
    score_table.to_csv(
        outdir / "07_cluster_number_silhouette.tsv",
        sep="\t",
        index=False,
    )

    gene_cols = [
        c for c in ["ME", "gene_id", "gene_name", "gene_biotype"]
        if c in gene.columns
    ]
    gene_small = gene[gene_cols].drop_duplicates(subset=["ME"])

    records = []
    for me, raw_cluster in zip(candidate_mes, raw_labels):
        row = matrix.loc[me]
        valid = row.notna()
        vals = row[valid]

        if len(vals):
            peak_tp = vals.abs().idxmax()
            peak_delta = float(vals.loc[peak_tp])
            dominant_direction = (
                "increase" if peak_delta > 0
                else "decrease" if peak_delta < 0
                else "zero"
            )
        else:
            peak_tp = ""
            peak_delta = np.nan
            dominant_direction = "missing"

        tc_me = tc_sub.loc[tc_sub["ME"] == me].copy()
        strong_mask = (
            as_bool(tc_me["valid_comparison"])
            & (pd.to_numeric(tc_me["abs_delta_PSI"], errors="coerce") >= 0.20)
            & as_bool(tc_me["complete_group_separation"])
        )
        strong_tps = [
            tp for tp in TIMEPOINT_ORDER
            if tp in set(tc_me.loc[strong_mask, "timepoint"])
        ]

        wh_me = wh.loc[(wh["ME"] == me) & as_bool(wh["whippet_mapped"])].copy()
        n_wh = len(wh_me)
        n_conc = int(as_bool(wh_me["direction_concordant"]).sum())

        records.append({
            "ME": me,
            "raw_cluster": int(raw_cluster),
            "n_valid_timepoints": int(valid.sum()),
            "peak_timepoint": peak_tp,
            "peak_delta_PSI": peak_delta,
            "peak_abs_delta_PSI": abs(peak_delta) if np.isfinite(peak_delta) else np.nan,
            "dominant_direction": dominant_direction,
            "first_strong_separated_timepoint": strong_tps[0] if strong_tps else "",
            "n_strong_separated_timepoints": len(strong_tps),
            "strong_separated_timepoints": ";".join(strong_tps),
            "whippet_mapped_candidate_responses": n_wh,
            "whippet_direction_concordant_responses": n_conc,
            "whippet_direction_concordance_fraction": (
                n_conc / n_wh if n_wh else np.nan
            ),
        })

    meta = pd.DataFrame(records).merge(gene_small, on="ME", how="left")

    tp_index = {tp: i for i, tp in enumerate(TIMEPOINT_ORDER)}
    raw_cluster_order = (
        meta.assign(
            peak_order=meta["peak_timepoint"].map(tp_index).fillna(999),
            direction_order=meta["dominant_direction"].map(
                {"increase": 0, "decrease": 1, "zero": 2, "missing": 3}
            ).fillna(9),
        )
        .groupby("raw_cluster", as_index=False)
        .agg(
            median_peak_order=("peak_order", "median"),
            median_direction_order=("direction_order", "median"),
            median_peak_delta=("peak_delta_PSI", "median"),
        )
        .sort_values(
            ["median_peak_order", "median_direction_order", "median_peak_delta"]
        )
    )
    remap = {
        int(raw): i + 1
        for i, raw in enumerate(raw_cluster_order["raw_cluster"].tolist())
    }
    meta["temporal_cluster"] = meta["raw_cluster"].map(remap).astype(int)
    meta = meta.drop(columns=["raw_cluster"])

    sort_meta = meta.copy()
    sort_meta["peak_order"] = sort_meta["peak_timepoint"].map(tp_index).fillna(999)
    sort_meta["gene_sort"] = sort_meta["gene_name"].fillna("")
    sort_meta = sort_meta.sort_values(
        ["temporal_cluster", "peak_order", "gene_sort", "ME"]
    )
    ordered_mes = sort_meta["ME"].tolist()

    meta.to_csv(
        outdir / "07_temporal_clusters.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary_rows = []
    for cluster in sorted(meta["temporal_cluster"].unique()):
        mes = meta.loc[meta["temporal_cluster"] == cluster, "ME"].tolist()
        sub = matrix.loc[mes]
        row = {
            "temporal_cluster": int(cluster),
            "n_microexons": len(mes),
            "n_peak_increase": int(
                (meta.loc[meta["temporal_cluster"] == cluster, "dominant_direction"]
                 == "increase").sum()
            ),
            "n_peak_decrease": int(
                (meta.loc[meta["temporal_cluster"] == cluster, "dominant_direction"]
                 == "decrease").sum()
            ),
            "median_whippet_direction_concordance_fraction": float(
                meta.loc[
                    meta["temporal_cluster"] == cluster,
                    "whippet_direction_concordance_fraction",
                ].median()
            ),
        }
        for tp in TIMEPOINT_ORDER:
            row[f"median_delta_{tp}"] = float(sub[tp].median(skipna=True))
        summary_rows.append(row)

    cluster_summary = pd.DataFrame(summary_rows)
    cluster_summary.to_csv(
        outdir / "07_temporal_cluster_summary.tsv",
        sep="\t",
        index=False,
    )

    onset_rows = []
    for tp in TIMEPOINT_ORDER:
        sub = meta.loc[meta["first_strong_separated_timepoint"] == tp]
        onset_rows.append({
            "first_strong_timepoint": tp,
            "n_microexons": len(sub),
            "n_peak_increase": int((sub["dominant_direction"] == "increase").sum()),
            "n_peak_decrease": int((sub["dominant_direction"] == "decrease").sum()),
        })
    onset_summary = pd.DataFrame(onset_rows)
    onset_summary.to_csv(
        outdir / "07_first_strong_response_summary.tsv",
        sep="\t",
        index=False,
    )

    peak_rows = []
    for tp in TIMEPOINT_ORDER:
        sub = meta.loc[meta["peak_timepoint"] == tp]
        peak_rows.append({
            "peak_timepoint": tp,
            "n_microexons": len(sub),
            "n_peak_increase": int((sub["dominant_direction"] == "increase").sum()),
            "n_peak_decrease": int((sub["dominant_direction"] == "decrease").sum()),
        })
    peak_summary = pd.DataFrame(peak_rows)
    peak_summary.to_csv(
        outdir / "07_peak_time_summary.tsv",
        sep="\t",
        index=False,
    )

    heat = matrix.loc[ordered_mes]
    label_df = sort_meta.set_index("ME")
    row_labels = []
    for me in ordered_mes:
        g = label_df.loc[me, "gene_name"] if "gene_name" in label_df.columns else ""
        gid = label_df.loc[me, "gene_id"] if "gene_id" in label_df.columns else ""
        if pd.isna(g) or str(g).strip() == "":
            g = gid if not pd.isna(gid) else "unannotated"
        row_labels.append(f"{g} | {me}")

    arr = heat.to_numpy(dtype=float)
    finite = np.abs(arr[np.isfinite(arr)])
    vmax = float(np.quantile(finite, 0.98)) if len(finite) else 1.0
    vmax = max(vmax, 0.20)

    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("#d9d9d9")
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    fig_h = max(15.0, 0.24 * len(ordered_mes))
    fig, ax = plt.subplots(figsize=(11, fig_h))
    im = ax.imshow(
        np.ma.masked_invalid(arr),
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        norm=norm,
    )
    ax.set_xticks(np.arange(len(TIMEPOINT_ORDER)))
    ax.set_xticklabels(["30 min", "3 h", "12 h", "24 h", "48 h", "72 h", "7 d"])
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=6)
    ax.set_xlabel("Time after CD3/CD28 activation")
    ax.set_ylabel("Whippet-matched candidate microexons")
    ax.set_title(
        f"Temporal microexon ΔPSI patterns (n={len(ordered_mes)})\n"
        "Raw ΔPSI shown; rows grouped by exploratory temporal clustering"
    )

    cluster_ids = sort_meta["temporal_cluster"].to_numpy()
    for i in range(1, len(cluster_ids)):
        if cluster_ids[i] != cluster_ids[i - 1]:
            ax.axhline(i - 0.5, color="black", linewidth=1.0)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("ΔPSI (Activation − Control)")
    ax.legend(
        handles=[Patch(facecolor="#d9d9d9", edgecolor="none", label="Insufficient coverage / NA")],
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(outdir / "07_temporal_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "07_temporal_heatmap.svg", bbox_inches="tight")
    plt.close(fig)

    onset_meta = meta.copy()
    onset_meta["onset_order"] = onset_meta[
        "first_strong_separated_timepoint"
    ].map(tp_index).fillna(999)
    onset_meta["peak_order"] = onset_meta["peak_timepoint"].map(tp_index).fillna(999)
    onset_meta["direction_order"] = onset_meta["dominant_direction"].map(
        {"increase": 0, "decrease": 1, "zero": 2, "missing": 3}
    ).fillna(9)
    onset_meta["gene_sort"] = onset_meta["gene_name"].fillna("")
    onset_meta = onset_meta.sort_values(
        ["onset_order", "direction_order", "peak_order", "gene_sort", "ME"]
    )
    onset_mes = onset_meta["ME"].tolist()
    onset_heat = matrix.loc[onset_mes]
    onset_labels = []
    onset_label_df = onset_meta.set_index("ME")
    for me in onset_mes:
        g = onset_label_df.loc[me, "gene_name"]
        gid = onset_label_df.loc[me, "gene_id"]
        if pd.isna(g) or str(g).strip() == "":
            g = gid if not pd.isna(gid) else "unannotated"
        onset_labels.append(f"{g} | {me}")

    onset_arr = onset_heat.to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(11, fig_h))
    im = ax.imshow(
        np.ma.masked_invalid(onset_arr),
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        norm=norm,
    )
    ax.set_xticks(np.arange(len(TIMEPOINT_ORDER)))
    ax.set_xticklabels(["30 min", "3 h", "12 h", "24 h", "48 h", "72 h", "7 d"])
    ax.set_yticks(np.arange(len(onset_labels)))
    ax.set_yticklabels(onset_labels, fontsize=6)
    ax.set_xlabel("Time after CD3/CD28 activation")
    ax.set_ylabel("Whippet-matched candidate microexons")
    ax.set_title(
        f"Temporal microexon ΔPSI patterns (n={len(onset_mes)})\n"
        "Rows ordered by first strong reproducible response"
    )

    onset_groups = onset_meta["first_strong_separated_timepoint"].to_numpy()
    for i in range(1, len(onset_groups)):
        if onset_groups[i] != onset_groups[i - 1]:
            ax.axhline(i - 0.5, color="black", linewidth=1.0)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("ΔPSI (Activation − Control)")
    ax.legend(
        handles=[Patch(facecolor="#d9d9d9", edgecolor="none", label="Insufficient coverage / NA")],
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(
        outdir / "07_temporal_heatmap_by_onset.png",
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        outdir / "07_temporal_heatmap_by_onset.svg",
        bbox_inches="tight",
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(TIMEPOINT_ORDER))
    for cluster in sorted(meta["temporal_cluster"].unique()):
        mes = meta.loc[meta["temporal_cluster"] == cluster, "ME"].tolist()
        sub = matrix.loc[mes].to_numpy(dtype=float)
        med = np.nanmedian(sub, axis=0)
        q1 = np.nanquantile(sub, 0.25, axis=0)
        q3 = np.nanquantile(sub, 0.75, axis=0)
        ax.plot(
            x,
            med,
            marker="o",
            linewidth=2.2,
            label=f"Cluster {cluster} (n={len(mes)})",
        )
        ax.fill_between(x, q1, q3, alpha=0.15)

    ax.axhline(0, linewidth=1.0, color="black")
    ax.set_xticks(x)
    ax.set_xticklabels(["30 min", "3 h", "12 h", "24 h", "48 h", "72 h", "7 d"])
    ax.set_xlabel("Time after CD3/CD28 activation")
    ax.set_ylabel("Median ΔPSI (Activation − Control)")
    ax.set_title("Exploratory temporal clusters of Whippet-matched microexons")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(
        outdir / "07_temporal_cluster_trajectories.png",
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        outdir / "07_temporal_cluster_trajectories.svg",
        bbox_inches="tight",
    )
    plt.close(fig)

    report = [
        "TEMPORAL PATTERN ANALYSIS - STEP 7",
        "",
        f"Unique Whippet-matched candidate microexons: {len(candidate_mes)}",
        f"Temporal matrix dimensions: {matrix.shape[0]} microexons x {matrix.shape[1]} timepoints",
        f"Observed valid DeltaPSI values: {int(matrix.notna().sum().sum())}",
        f"Missing/insufficient-coverage DeltaPSI values: {int(matrix.isna().sum().sum())}",
        "",
        "Temporal clustering:",
        "  Missing values are not replaced in the reported matrix or heatmap.",
        "  For clustering only, missing timepoints are linearly interpolated",
        "  across actual sampling times, with edge gaps filled by the nearest",
        "  observed value. Each trajectory is then row-standardized so clustering",
        "  emphasizes temporal shape rather than absolute effect magnitude.",
        "  Ward hierarchical clustering is evaluated for k=2..6, and the k with",
        "  the highest mean silhouette score is selected.",
        f"Chosen number of temporal clusters: {chosen_k}",
        f"Mean silhouette score: {silhouette}",
        "",
        "For timing-oriented interpretation, a second heatmap orders rows by",
        "the first timepoint meeting the existing Step 5 criteria:",
        "valid comparison, |DeltaPSI| >= 0.20, and complete group separation.",
        "This uses no new response threshold and makes early versus late onset",
        "directly visible without forcing extra temporal clusters.",
        "",
        "The heatmap displays raw MicroExonator DeltaPSI values for all matched",
        "candidate microexons. Gray cells indicate comparisons that did not meet",
        "the Step 2 minimum-valid-replicate requirement.",
        "",
        "Whippet information is used to define the matched candidate set and",
        "to summarize direction concordance; the temporal trajectories themselves",
        "use the full MicroExonator time course so non-Step-5 timepoints remain visible.",
    ]
    (outdir / "07_temporal_analysis_report.txt").write_text(
        "\n".join(report) + "\n",
        encoding="utf-8",
    )

    print("\n".join(report))
    print("\nCluster summary:")
    print(cluster_summary.to_string(index=False))


if __name__ == "__main__":
    main()
