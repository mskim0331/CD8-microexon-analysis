#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


TIMEPOINT_ORDER = ["30min", "3hr", "12hr", "24hr", "48hr", "72hr", "7day"]
DISPLAY_TIMEPOINTS = ["30 min", "3 h", "12 h", "24 h", "48 h", "72 h", "7 d"]
PHASE_MAP = {
    "30min": "early",
    "3hr": "early",
    "12hr": "intermediate",
    "24hr": "intermediate",
    "48hr": "late",
    "72hr": "late",
    "7day": "late",
}
PHASE_ORDER = ["early", "intermediate", "late"]


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Finalize Step 7 temporal response descriptors, robustness checks, "
            "and representative PSI trajectories."
        )
    )
    p.add_argument(
        "--timecourse",
        default="analysis_results/02_ctrl_vs_timepoint_deltaPSI.tsv",
    )
    p.add_argument(
        "--step5",
        default="analysis_results/02_strong_responses_complete_separation.tsv",
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


def strong_status(row):
    valid = str(row["valid_comparison"]).lower() == "true"
    sep = str(row["complete_group_separation"]).lower() == "true"
    delta = pd.to_numeric(row["delta_PSI_activation_minus_ctrl"], errors="coerce")
    if not valid or not sep or pd.isna(delta) or abs(delta) < 0.20:
        return 0
    return 1 if delta > 0 else -1 if delta < 0 else 0


def longest_consecutive(indices):
    if not indices:
        return 0
    best = cur = 1
    for a, b in zip(indices[:-1], indices[1:]):
        if b == a + 1:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def phase_label(tp):
    return PHASE_MAP.get(tp, "unknown")


def pattern_family(onset_phase, direction, n_strong, strong_at_7day, mixed):
    if mixed:
        return "direction-switching"
    if n_strong == 1:
        return f"{onset_phase} single-timepoint {direction}"
    if strong_at_7day:
        return f"{onset_phase} multi-timepoint, retained at 7 d {direction}"
    return f"{onset_phase} multi-timepoint, not retained at 7 d {direction}"


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    tc = pd.read_csv(args.timecourse, sep="\t")
    step5 = pd.read_csv(args.step5, sep="\t")
    wh = pd.read_csv(args.whippet, sep="\t")
    gene = pd.read_csv(args.gene_annotation, sep="\t", dtype=str)

    for df in [tc, step5, wh]:
        df["timepoint"] = df["timepoint"].map(normalize_timepoint)
        if df["timepoint"].isna().any():
            raise ValueError("Unrecognized timepoint in an input table.")

    mapped_wh = wh.loc[as_bool(wh["whippet_mapped"])].copy()
    matched_mes = sorted(mapped_wh["ME"].dropna().unique())
    if len(matched_mes) != 75:
        raise ValueError(
            f"Expected 75 Whippet-matched candidate microexons, found {len(matched_mes)}."
        )

    gene_cols = [
        c for c in ["ME", "gene_id", "gene_name", "gene_biotype"]
        if c in gene.columns
    ]
    gene_small = gene[gene_cols].drop_duplicates("ME")

    rows = []
    for me in matched_mes:
        sub = (
            tc.loc[tc["ME"] == me]
            .copy()
            .set_index("timepoint")
            .reindex(TIMEPOINT_ORDER)
        )
        statuses = [strong_status(sub.loc[tp]) for tp in TIMEPOINT_ORDER]
        strong_idx = [i for i, s in enumerate(statuses) if s != 0]
        strong_tps = [TIMEPOINT_ORDER[i] for i in strong_idx]
        signs = sorted(set(statuses[i] for i in strong_idx))
        mixed = len(signs) > 1

        if strong_idx:
            onset_tp = TIMEPOINT_ORDER[strong_idx[0]]
            onset_phase = phase_label(onset_tp)
            direction = (
                "increase" if signs == [1]
                else "decrease" if signs == [-1]
                else "mixed"
            )
        else:
            onset_tp = ""
            onset_phase = "none"
            direction = "none"

        valid = as_bool(sub["valid_comparison"].fillna(False))
        delta = pd.to_numeric(
            sub["delta_PSI_activation_minus_ctrl"], errors="coerce"
        )
        valid_delta = delta.loc[valid & delta.notna()]
        if len(valid_delta):
            peak_tp = valid_delta.abs().idxmax()
            peak_delta = float(valid_delta.loc[peak_tp])
        else:
            peak_tp = ""
            peak_delta = np.nan

        wh_me = mapped_wh.loc[mapped_wh["ME"] == me]
        n_wh = len(wh_me)
        n_conc = int(as_bool(wh_me["direction_concordant"]).sum())

        strong_at_7d = statuses[-1] != 0
        n_strong = len(strong_idx)

        rows.append({
            "ME": me,
            "onset_timepoint": onset_tp,
            "onset_phase": onset_phase,
            "strong_direction": direction,
            "n_strong_timepoints": n_strong,
            "strong_timepoints": ";".join(strong_tps),
            "longest_consecutive_strong_run": longest_consecutive(strong_idx),
            "strong_at_7day": strong_at_7d,
            "mixed_strong_directions": mixed,
            "pattern_family": pattern_family(
                onset_phase, direction, n_strong, strong_at_7d, mixed
            ),
            "n_valid_timepoints": int(valid.sum()),
            "peak_timepoint": peak_tp,
            "peak_delta_PSI": peak_delta,
            "peak_abs_delta_PSI": abs(peak_delta) if np.isfinite(peak_delta) else np.nan,
            "whippet_mapped_candidate_responses": n_wh,
            "whippet_direction_concordant_responses": n_conc,
            "whippet_direction_concordance_fraction": (
                n_conc / n_wh if n_wh else np.nan
            ),
        })

    patterns = pd.DataFrame(rows).merge(gene_small, on="ME", how="left")
    patterns.to_csv(
        outdir / "07_temporal_pattern_labels.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary = (
        patterns.groupby(
            ["onset_phase", "strong_direction", "strong_at_7day"],
            dropna=False,
            as_index=False,
        )
        .agg(
            n_microexons=("ME", "nunique"),
            median_n_strong_timepoints=("n_strong_timepoints", "median"),
            median_peak_abs_delta_PSI=("peak_abs_delta_PSI", "median"),
        )
    )
    summary["phase_order"] = summary["onset_phase"].map(
        {p: i for i, p in enumerate(PHASE_ORDER)}
    ).fillna(99)
    summary["direction_order"] = summary["strong_direction"].map(
        {"increase": 0, "decrease": 1, "mixed": 2, "none": 3}
    ).fillna(9)
    summary = summary.sort_values(
        ["phase_order", "direction_order", "strong_at_7day"],
        ascending=[True, True, False],
    ).drop(columns=["phase_order", "direction_order"])
    summary.to_csv(
        outdir / "07_temporal_pattern_summary.tsv",
        sep="\t",
        index=False,
    )

    # Sensitivity analysis: retain the same Step 5 effect and separation criteria
    # but require all 3 biological replicates to be non-missing in both groups.
    step5_all3 = step5.loc[as_bool(step5["all_3_replicates_both_groups"])].copy()
    mapped_keys = set(
        zip(mapped_wh["ME"].astype(str), mapped_wh["timepoint"].astype(str))
    )
    step5_all3["whippet_mapped"] = [
        (str(me), str(tp)) in mapped_keys
        for me, tp in zip(step5_all3["ME"], step5_all3["timepoint"])
    ]

    robustness_rows = []
    for tp in TIMEPOINT_ORDER:
        base_tp = step5.loc[step5["timepoint"] == tp]
        strict_tp = step5_all3.loc[step5_all3["timepoint"] == tp]
        robustness_rows.append({
            "timepoint": tp,
            "original_step5_responses_min2of3": len(base_tp),
            "strict_all3_responses": len(strict_tp),
            "strict_fraction_of_original": (
                len(strict_tp) / len(base_tp) if len(base_tp) else np.nan
            ),
            "strict_whippet_mapped_responses": int(strict_tp["whippet_mapped"].sum()),
        })
    robustness = pd.DataFrame(robustness_rows)
    robustness.to_csv(
        outdir / "07_replicate_coverage_sensitivity.tsv",
        sep="\t",
        index=False,
    )

    robustness_overall = pd.DataFrame([{
        "original_step5_responses": len(step5),
        "original_step5_unique_microexons": int(step5["ME"].nunique()),
        "strict_all3_responses": len(step5_all3),
        "strict_all3_unique_microexons": int(step5_all3["ME"].nunique()),
        "strict_all3_whippet_mapped_responses": int(step5_all3["whippet_mapped"].sum()),
        "strict_all3_whippet_mapped_unique_microexons": int(
            step5_all3.loc[step5_all3["whippet_mapped"], "ME"].nunique()
        ),
    }])
    robustness_overall.to_csv(
        outdir / "07_replicate_coverage_sensitivity_overall.tsv",
        sep="\t",
        index=False,
    )

    # Data-driven representative trajectories. Selection is deliberately
    # independent of literature: complete time-course coverage, perfect observed
    # Whippet direction concordance, >=2 strong timepoints, then maximum response
    # breadth and effect magnitude within each onset-phase/direction stratum.
    reps = []
    eligible = patterns.loc[
        (patterns["n_valid_timepoints"] == 7)
        & (patterns["whippet_direction_concordance_fraction"] == 1.0)
        & (patterns["n_strong_timepoints"] >= 2)
        & patterns["gene_name"].notna()
        & patterns["gene_name"].astype(str).ne("")
        & patterns["onset_phase"].isin(PHASE_ORDER)
        & patterns["strong_direction"].isin(["increase", "decrease"])
    ].copy()

    for phase in PHASE_ORDER:
        for direction in ["increase", "decrease"]:
            sub = eligible.loc[
                (eligible["onset_phase"] == phase)
                & (eligible["strong_direction"] == direction)
            ].copy()
            if sub.empty:
                continue
            sub = sub.sort_values(
                ["n_strong_timepoints", "peak_abs_delta_PSI", "ME"],
                ascending=[False, False, True],
            )
            chosen = sub.iloc[0].to_dict()
            chosen["representative_stratum"] = f"{phase}_{direction}"
            reps.append(chosen)

    reps = pd.DataFrame(reps)
    reps.to_csv(
        outdir / "07_data_driven_representative_candidates.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # Figure 1: onset-phase/direction counts.
    count_table = (
        patterns.loc[
            patterns["onset_phase"].isin(PHASE_ORDER)
            & patterns["strong_direction"].isin(["increase", "decrease"])
        ]
        .groupby(["onset_phase", "strong_direction"])
        .size()
        .unstack(fill_value=0)
        .reindex(PHASE_ORDER)
    )

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(PHASE_ORDER))
    inc = count_table.get("increase", pd.Series(0, index=PHASE_ORDER)).to_numpy()
    dec = count_table.get("decrease", pd.Series(0, index=PHASE_ORDER)).to_numpy()
    ax.bar(x, inc, label="Increase")
    ax.bar(x, dec, bottom=inc, label="Decrease")
    ax.set_xticks(x)
    ax.set_xticklabels(["Early\n30 min–3 h", "Intermediate\n12–24 h", "Late\n48 h–7 d"])
    ax.set_ylabel("Number of microexons")
    ax.set_title("Timing and direction of first strong reproducible response")
    ax.legend(frameon=False)
    for i, (a, b) in enumerate(zip(inc, dec)):
        ax.text(i, a + b + 0.5, str(int(a + b)), ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "07_onset_phase_direction_counts.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "07_onset_phase_direction_counts.svg", bbox_inches="tight")
    plt.close(fig)

    # Figure 2: six representative PSI trajectories (when all strata are available).
    if len(reps):
        n = len(reps)
        ncols = 2
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(
            nrows, ncols, figsize=(12, 3.8 * nrows), squeeze=False
        )
        x = np.arange(8)
        xlabels = ["Ctrl"] + DISPLAY_TIMEPOINTS

        for ax, (_, rep) in zip(axes.flat, reps.iterrows()):
            me = rep["ME"]
            sub = (
                tc.loc[tc["ME"] == me]
                .copy()
                .set_index("timepoint")
                .reindex(TIMEPOINT_ORDER)
            )

            ctrl_mean = pd.to_numeric(sub["ctrl_mean_PSI"], errors="coerce").dropna()
            ctrl_sd = pd.to_numeric(sub["ctrl_sd_PSI"], errors="coerce").dropna()
            ctrl_y = float(ctrl_mean.iloc[0]) if len(ctrl_mean) else np.nan
            ctrl_err = float(ctrl_sd.iloc[0]) if len(ctrl_sd) else np.nan

            act_y = pd.to_numeric(sub["activation_mean_PSI"], errors="coerce").to_numpy()
            act_err = pd.to_numeric(sub["activation_sd_PSI"], errors="coerce").to_numpy()
            y = np.concatenate([[ctrl_y], act_y])
            yerr = np.concatenate([[ctrl_err], act_err])

            ax.errorbar(x, y, yerr=yerr, marker="o", capsize=3, linewidth=1.8)
            ax.set_xticks(x)
            ax.set_xticklabels(xlabels, rotation=45, ha="right")
            ax.set_ylim(-0.05, 1.05)
            ax.set_ylabel("PSI")
            ax.set_title(
                f"{rep['gene_name']} | {rep['representative_stratum']}\n{me}",
                fontsize=10,
            )
            ax.grid(axis="y", alpha=0.2)

        for ax in axes.flat[n:]:
            ax.axis("off")

        fig.suptitle(
            "Data-driven representative microexon PSI trajectories",
            y=1.01,
            fontsize=14,
        )
        fig.tight_layout()
        fig.savefig(
            outdir / "07_representative_candidate_PSI_trajectories.png",
            dpi=300,
            bbox_inches="tight",
        )
        fig.savefig(
            outdir / "07_representative_candidate_PSI_trajectories.svg",
            bbox_inches="tight",
        )
        plt.close(fig)

    report = [
        "STEP 7 FINAL TEMPORAL SUMMARY",
        "",
        f"Whippet-matched candidate microexons analyzed: {len(patterns)}",
        "",
        "Response definition:",
        "  A strong reproducible response uses the existing Step 5 definition:",
        "  valid comparison, |DeltaPSI| >= 0.20, and complete group separation.",
        "  No new DeltaPSI threshold is introduced in Step 7.",
        "",
        "Timing descriptors:",
        "  early = first strong response at 30 min or 3 h",
        "  intermediate = first strong response at 12 h or 24 h",
        "  late = first strong response at 48 h, 72 h, or 7 d",
        "  'retained at 7 d' means the 7 d comparison itself meets Step 5 criteria;",
        "  it does not imply that every intervening timepoint was significant/strong.",
        "",
        "Robustness check:",
        "  The original Step 5 allowed >=2/3 valid biological replicates per group.",
        "  A sensitivity analysis re-counts the same responses requiring all 3/3",
        "  replicates in both control and activation groups.",
        "",
        "Representative trajectories:",
        "  Selection is data-driven, not literature-driven: complete 7-timepoint",
        "  coverage, 100% observed Whippet direction concordance, >=2 strong",
        "  timepoints, then response breadth and peak effect size within each",
        "  onset-phase/direction stratum.",
    ]
    (outdir / "07_final_temporal_summary_report.txt").write_text(
        "\n".join(report) + "\n",
        encoding="utf-8",
    )

    print("\n".join(report))
    print("\nPattern summary:")
    print(summary.to_string(index=False))
    print("\nRobustness overall:")
    print(robustness_overall.to_string(index=False))
    print("\nRepresentative candidates:")
    if len(reps):
        print(
            reps[
                [
                    "representative_stratum",
                    "gene_name",
                    "ME",
                    "n_strong_timepoints",
                    "peak_abs_delta_PSI",
                ]
            ].to_string(index=False)
        )


if __name__ == "__main__":
    main()
