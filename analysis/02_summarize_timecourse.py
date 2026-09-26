#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


CONDITION_ORDER = [
    "Ctrl",
    "30min",
    "3hr",
    "12hr",
    "24hr",
    "48hr",
    "72hr",
    "7day",
]

FOCUS_MICROEXONS = {
    "8_-_27290155_27290178": "TRIM35",
    "3_-_108050577_108050602": "CD47",
    "19_+_18122106_18122127": "MAST3",
    "5_+_57221375_57221396": "GPBP1",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--matrix",
        default="analysis_results/01_high_quality_PSI_matrix.tsv",
    )
    p.add_argument(
        "--annotation",
        default="analysis_results/01_high_quality_annotation.tsv",
    )
    p.add_argument("--outdir", default="analysis_results")
    p.add_argument("--min-valid-replicates", type=int, default=2)
    p.add_argument("--exploratory-delta", type=float, default=0.10)
    p.add_argument("--strong-delta", type=float, default=0.20)
    return p.parse_args()


def condition_from_sample(sample):
    if sample.startswith("CD8_ctrl_"):
        return "Ctrl"
    for condition in CONDITION_ORDER[1:]:
        if sample.startswith(f"CD8_{condition}_activation_"):
            return condition
    raise ValueError(f"Unrecognized sample name: {sample}")


def summarize_values(values):
    values = pd.Series(values, dtype=float).dropna()
    n = len(values)

    return {
        "n_valid": n,
        "mean_PSI": values.mean() if n else np.nan,
        "median_PSI": values.median() if n else np.nan,
        "sd_PSI": values.std(ddof=1) if n >= 2 else np.nan,
        "min_PSI": values.min() if n else np.nan,
        "max_PSI": values.max() if n else np.nan,
        "range_PSI": (values.max() - values.min()) if n else np.nan,
    }


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    annotation = pd.read_csv(args.annotation, sep="\t", dtype={"ME": str})

    if matrix.index.has_duplicates:
        raise ValueError("PSI matrix has duplicated microexon coordinates.")

    sample_conditions = {
        sample: condition_from_sample(sample)
        for sample in matrix.columns
    }

    condition_samples = {
        condition: [
            sample for sample in matrix.columns
            if sample_conditions[sample] == condition
        ]
        for condition in CONDITION_ORDER
    }

    for condition in CONDITION_ORDER:
        if len(condition_samples[condition]) != 3:
            raise ValueError(
                f"Expected 3 samples for {condition}, found "
                f"{len(condition_samples[condition])}: "
                f"{condition_samples[condition]}"
            )

    condition_records = []
    for me, row in matrix.iterrows():
        for condition in CONDITION_ORDER:
            stats = summarize_values(row[condition_samples[condition]])
            condition_records.append({
                "ME": me,
                "condition": condition,
                **stats,
            })

    condition_summary = pd.DataFrame(condition_records)

    # One control summary per microexon, reused for all seven comparisons.
    ctrl = (
        condition_summary.loc[condition_summary["condition"] == "Ctrl"]
        .set_index("ME")
    )

    delta_records = []
    for condition in CONDITION_ORDER[1:]:
        act = (
            condition_summary.loc[
                condition_summary["condition"] == condition
            ]
            .set_index("ME")
        )

        for me in matrix.index:
            c = ctrl.loc[me]
            a = act.loc[me]

            valid = (
                c["n_valid"] >= args.min_valid_replicates
                and a["n_valid"] >= args.min_valid_replicates
            )

            delta = np.nan
            abs_delta = np.nan
            complete_separation = False
            effect_class = "not_testable"

            if valid:
                delta = a["mean_PSI"] - c["mean_PSI"]
                abs_delta = abs(delta)

                if delta > 0:
                    complete_separation = a["min_PSI"] > c["max_PSI"]
                elif delta < 0:
                    complete_separation = a["max_PSI"] < c["min_PSI"]

                if abs_delta >= args.strong_delta:
                    effect_class = "strong"
                elif abs_delta >= args.exploratory_delta:
                    effect_class = "exploratory"
                else:
                    effect_class = "small"

            delta_records.append({
                "ME": me,
                "timepoint": condition,
                "ctrl_n_valid": int(c["n_valid"]),
                "activation_n_valid": int(a["n_valid"]),
                "ctrl_mean_PSI": c["mean_PSI"],
                "activation_mean_PSI": a["mean_PSI"],
                "delta_PSI_activation_minus_ctrl": delta,
                "abs_delta_PSI": abs_delta,
                "ctrl_sd_PSI": c["sd_PSI"],
                "activation_sd_PSI": a["sd_PSI"],
                "ctrl_min_PSI": c["min_PSI"],
                "ctrl_max_PSI": c["max_PSI"],
                "activation_min_PSI": a["min_PSI"],
                "activation_max_PSI": a["max_PSI"],
                "valid_comparison": valid,
                "all_3_replicates_both_groups": (
                    c["n_valid"] == 3 and a["n_valid"] == 3
                ),
                "complete_group_separation": complete_separation,
                "effect_class": effect_class,
            })

    delta = pd.DataFrame(delta_records)

    exploratory = delta.loc[
        delta["valid_comparison"]
        & (delta["abs_delta_PSI"] >= args.exploratory_delta)
    ].copy()

    strong = delta.loc[
        delta["valid_comparison"]
        & (delta["abs_delta_PSI"] >= args.strong_delta)
    ].copy()

    strong_separated = strong.loc[
        strong["complete_group_separation"]
    ].copy()

    # Summarize how dynamic each microexon is over the complete time course.
    dynamic_records = []
    for me, sub in delta.groupby("ME", sort=False):
        valid_sub = sub.loc[sub["valid_comparison"]].copy()

        if valid_sub.empty:
            dynamic_records.append({
                "ME": me,
                "n_testable_timepoints": 0,
                "n_exploratory_timepoints": 0,
                "n_strong_timepoints": 0,
                "n_strong_separated_timepoints": 0,
                "max_abs_delta_PSI": np.nan,
                "timepoint_of_max_abs_delta": "",
                "delta_at_max": np.nan,
                "direction_at_max": "",
            })
            continue

        idx = valid_sub["abs_delta_PSI"].idxmax()
        top = valid_sub.loc[idx]
        d = top["delta_PSI_activation_minus_ctrl"]

        dynamic_records.append({
            "ME": me,
            "n_testable_timepoints": len(valid_sub),
            "n_exploratory_timepoints": int(
                (valid_sub["abs_delta_PSI"] >= args.exploratory_delta).sum()
            ),
            "n_strong_timepoints": int(
                (valid_sub["abs_delta_PSI"] >= args.strong_delta).sum()
            ),
            "n_strong_separated_timepoints": int(
                (
                    (valid_sub["abs_delta_PSI"] >= args.strong_delta)
                    & valid_sub["complete_group_separation"]
                ).sum()
            ),
            "max_abs_delta_PSI": top["abs_delta_PSI"],
            "timepoint_of_max_abs_delta": top["timepoint"],
            "delta_at_max": d,
            "direction_at_max": (
                "increased" if d > 0
                else "decreased" if d < 0
                else "unchanged"
            ),
        })

    dynamic = pd.DataFrame(dynamic_records)
    dynamic = annotation.merge(dynamic, on="ME", how="right")

    # Counts by time point make the effect-size filters auditable.
    count_records = []
    for condition in CONDITION_ORDER[1:]:
        sub = delta.loc[delta["timepoint"] == condition]
        valid = sub.loc[sub["valid_comparison"]]

        count_records.append({
            "timepoint": condition,
            "total_microexons": len(sub),
            "testable_microexons": len(valid),
            "not_testable_microexons": int((~sub["valid_comparison"]).sum()),
            "abs_delta_ge_0.10": int(
                (valid["abs_delta_PSI"] >= args.exploratory_delta).sum()
            ),
            "abs_delta_ge_0.20": int(
                (valid["abs_delta_PSI"] >= args.strong_delta).sum()
            ),
            "abs_delta_ge_0.20_and_complete_separation": int(
                (
                    (valid["abs_delta_PSI"] >= args.strong_delta)
                    & valid["complete_group_separation"]
                ).sum()
            ),
            "increased_ge_0.20": int(
                (
                    valid["delta_PSI_activation_minus_ctrl"]
                    >= args.strong_delta
                ).sum()
            ),
            "decreased_le_minus_0.20": int(
                (
                    valid["delta_PSI_activation_minus_ctrl"]
                    <= -args.strong_delta
                ).sum()
            ),
        })

    filter_counts = pd.DataFrame(count_records)

    # Recheck the four candidates identified in the earlier analysis only
    # after the full high-confidence set has been processed.
    focus_rows = delta.loc[
        delta["ME"].isin(FOCUS_MICROEXONS)
    ].copy()
    focus_rows.insert(
        1,
        "gene_label_from_previous_analysis",
        focus_rows["ME"].map(FOCUS_MICROEXONS),
    )

    condition_summary.to_csv(
        outdir / "02_condition_PSI_summary.tsv",
        sep="\t",
        index=False,
    )
    delta.to_csv(
        outdir / "02_ctrl_vs_timepoint_deltaPSI.tsv",
        sep="\t",
        index=False,
    )
    exploratory.to_csv(
        outdir / "02_exploratory_responses_absDelta_ge_0.10.tsv",
        sep="\t",
        index=False,
    )
    strong.to_csv(
        outdir / "02_strong_responses_absDelta_ge_0.20.tsv",
        sep="\t",
        index=False,
    )
    strong_separated.to_csv(
        outdir / "02_strong_responses_complete_separation.tsv",
        sep="\t",
        index=False,
    )
    dynamic.to_csv(
        outdir / "02_microexon_timecourse_summary.tsv",
        sep="\t",
        index=False,
    )
    filter_counts.to_csv(
        outdir / "02_filter_counts_by_timepoint.tsv",
        sep="\t",
        index=False,
    )
    focus_rows.to_csv(
        outdir / "02_previous_candidate_recheck.tsv",
        sep="\t",
        index=False,
    )

    report = []
    report.append("MICROEXON TIME-COURSE ANALYSIS - STEP 2")
    report.append("")
    report.append(
        f"Minimum valid replicates per group: {args.min_valid_replicates}/3"
    )
    report.append(
        f"Exploratory effect-size threshold: |delta PSI| >= "
        f"{args.exploratory_delta:.2f}"
    )
    report.append(
        f"Strong effect-size threshold: |delta PSI| >= "
        f"{args.strong_delta:.2f}"
    )
    report.append("")
    report.append(
        "Delta PSI convention: activation mean PSI - control mean PSI"
    )
    report.append("  positive = increased inclusion after activation")
    report.append("  negative = decreased inclusion after activation")
    report.append("")
    report.append("Counts by time point:")

    for _, row in filter_counts.iterrows():
        report.append(
            f"  {row['timepoint']}: "
            f"testable={int(row['testable_microexons'])}, "
            f"|dPSI|>=0.10={int(row['abs_delta_ge_0.10'])}, "
            f"|dPSI|>=0.20={int(row['abs_delta_ge_0.20'])}, "
            f"strong+separated="
            f"{int(row['abs_delta_ge_0.20_and_complete_separation'])}"
        )

    report.append("")
    report.append("Previous candidate recheck:")
    for me, gene in FOCUS_MICROEXONS.items():
        sub = focus_rows.loc[focus_rows["ME"] == me]
        if sub.empty:
            report.append(f"  {gene} ({me}): not present")
            continue

        best = sub.loc[sub["abs_delta_PSI"].idxmax()]
        report.append(
            f"  {gene} ({me}): max |dPSI|="
            f"{best['abs_delta_PSI']:.4f} at {best['timepoint']}, "
            f"dPSI={best['delta_PSI_activation_minus_ctrl']:.4f}"
        )

    (outdir / "02_timecourse_report.txt").write_text(
        "\n".join(report) + "\n",
        encoding="utf-8",
    )

    print("\n".join(report))


if __name__ == "__main__":
    main()
