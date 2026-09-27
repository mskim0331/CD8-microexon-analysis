#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


TIMEPOINT_ORDER = ["30min", "3hr", "12hr", "24hr", "48hr", "72hr", "7day"]
TIMEPOINT_ALIASES = {
    "30min": ["30min"],
    "3hr": ["3hr", "3h"],
    "12hr": ["12hr", "12h"],
    "24hr": ["24hr", "24h"],
    "48hr": ["48hr", "48h"],
    "72hr": ["72hr", "72h"],
    "7day": ["7day", "7d"],
}


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Validate MicroExonator candidate responses with paired-end Whippet "
            "differential-splicing results."
        )
    )
    p.add_argument(
        "--candidates",
        default="analysis_results/02_strong_responses_complete_separation.tsv",
    )
    p.add_argument(
        "--gene-annotation",
        default="analysis_results/03_microexon_gene_annotation.tsv",
    )
    p.add_argument(
        "--whippet-exons",
        required=True,
        help="Whippet index exon table: whippet.jls.exons.tab.gz",
    )
    p.add_argument(
        "--delta-root",
        required=True,
        help="Directory containing the seven paired-end Ctrl-vs-activation .diff.gz files.",
    )
    p.add_argument(
        "--outdir",
        default="analysis_results",
    )
    return p.parse_args()


def normalize_timepoint(value):
    s = str(value).lower().replace("_", "").replace("-", "")
    for canonical, aliases in TIMEPOINT_ALIASES.items():
        for alias in aliases:
            if alias.replace("_", "").replace("-", "") in s:
                return canonical
    return None


def discover_delta_files(root):
    root = Path(root)
    files = sorted(root.rglob("*.diff.gz"))
    if not files:
        raise FileNotFoundError(f"No .diff.gz files found under: {root}")

    by_timepoint = {tp: [] for tp in TIMEPOINT_ORDER}
    for path in files:
        name = path.name.lower()
        # Step 6 must use control-referenced comparisons only.
        if "ctrl" not in name and "control" not in name:
            continue
        tp = normalize_timepoint(name)
        if tp:
            by_timepoint[tp].append(path)

    problems = []
    selected = {}
    for tp in TIMEPOINT_ORDER:
        matches = by_timepoint[tp]
        if len(matches) == 1:
            selected[tp] = matches[0]
        elif len(matches) == 0:
            problems.append(f"{tp}: no Ctrl-vs-{tp} .diff.gz found")
        else:
            problems.append(
                f"{tp}: multiple candidate files found: "
                + ", ".join(str(x) for x in matches)
            )

    if problems:
        raise RuntimeError(
            "Could not uniquely identify the seven paired-end Whippet comparisons:\n  "
            + "\n  ".join(problems)
        )

    return selected


def load_whippet_node_exons(path):
    node_exons = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Potential_Exon", "Whippet_Nodes", "Gene", "Is_Annotated"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Missing columns in Whippet exon table: {sorted(missing)}"
            )

        for row in reader:
            for node in row["Whippet_Nodes"].split(","):
                node_exons[(row["Gene"], node)] = (
                    row["Potential_Exon"],
                    row["Is_Annotated"],
                )
    return node_exons


def whippet_row_to_microexon_id(row, node_exons):
    # Reproduces the coordinate conversion used by the pinned
    # MicroExonator whippet_delta_to_ME.py script.
    chrom, pos = row["Coord"].split(":")
    estart, eend = pos.split("-")
    exon_id = "_".join([chrom, row["Strand"], str(int(estart) - 1), eend])

    key = (row["Gene"], row["Node"])
    if key not in node_exons:
        return exon_id

    potential_exon, _ = node_exons[key]

    if row["Type"] == "AD":
        _, _, nstart, nend = exon_id.split("_")
        echrom, eloci, estrand = potential_exon.split(":")
        pstart, pend = eloci.split("-")

        if estrand == "+" and pend == nend:
            exon_id = "_".join(
                [echrom, estrand, str(int(pstart) - 1), pend]
            )
        if estrand == "-" and str(int(pstart) - 1) == nstart:
            exon_id = "_".join(
                [echrom, estrand, str(int(pstart) - 1), pend]
            )

    elif row["Type"] == "AA":
        _, _, nstart, nend = exon_id.split("_")
        echrom, eloci, estrand = potential_exon.split(":")
        pstart, pend = eloci.split("-")

        if estrand == "-" and pend == nend:
            exon_id = "_".join(
                [echrom, estrand, str(int(pstart) - 1), pend]
            )
        if estrand == "+" and str(int(pstart) - 1) == nstart:
            exon_id = "_".join(
                [echrom, estrand, str(int(pstart) - 1), pend]
            )

    return exon_id


def extract_whippet_candidates(diff_path, timepoint, candidate_mes, node_exons):
    records = []
    with gzip.open(diff_path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Gene", "Node", "Coord", "Strand", "Type",
            "Psi_A", "Psi_B", "DeltaPsi", "Probability"
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Missing columns in {diff_path}: {sorted(missing)}"
            )

        for row in reader:
            me = whippet_row_to_microexon_id(row, node_exons)
            if me not in candidate_mes:
                continue

            raw_delta = pd.to_numeric(row["DeltaPsi"], errors="coerce")
            corrected_delta = -raw_delta if pd.notna(raw_delta) else np.nan

            records.append({
                "ME": me,
                "timepoint": timepoint,
                "whippet_Gene": row["Gene"],
                "whippet_Node": row["Node"],
                "whippet_Coord": row["Coord"],
                "whippet_Strand": row["Strand"],
                "whippet_Type": row["Type"],
                "whippet_Psi_A": pd.to_numeric(row["Psi_A"], errors="coerce"),
                "whippet_Psi_B": pd.to_numeric(row["Psi_B"], errors="coerce"),
                "whippet_delta_raw_ctrl_minus_activation": raw_delta,
                "whippet_delta_activation_minus_ctrl": corrected_delta,
                "whippet_Probability": pd.to_numeric(
                    row["Probability"], errors="coerce"
                ),
                "whippet_diff_file": str(diff_path),
            })
    return records


def sign_label(x):
    if pd.isna(x):
        return "missing"
    if x > 0:
        return "increase"
    if x < 0:
        return "decrease"
    return "zero"


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(args.candidates, sep="\t")
    gene = pd.read_csv(args.gene_annotation, sep="\t", dtype=str)

    required_candidates = {
        "ME", "timepoint", "delta_PSI_activation_minus_ctrl",
        "abs_delta_PSI", "valid_comparison", "complete_group_separation"
    }
    missing = required_candidates - set(candidates.columns)
    if missing:
        raise ValueError(
            f"Candidate table is missing columns: {sorted(missing)}"
        )

    candidates["timepoint"] = candidates["timepoint"].map(normalize_timepoint)
    if candidates["timepoint"].isna().any():
        bad = candidates.loc[candidates["timepoint"].isna()]
        raise ValueError(
            "Unrecognized timepoint(s) in candidate table: "
            + ", ".join(sorted(set(bad["timepoint"].astype(str))))
        )

    # The Step 5 file is already filtered, but enforce the intended criteria.
    step5 = candidates.loc[
        candidates["valid_comparison"].astype(str).str.lower().eq("true")
        & (pd.to_numeric(candidates["abs_delta_PSI"], errors="coerce") >= 0.20)
        & candidates["complete_group_separation"].astype(str).str.lower().eq("true")
    ].copy()

    candidate_mes = set(step5["ME"])
    delta_files = discover_delta_files(args.delta_root)
    node_exons = load_whippet_node_exons(args.whippet_exons)

    whippet_records = []
    for tp in TIMEPOINT_ORDER:
        whippet_records.extend(
            extract_whippet_candidates(
                delta_files[tp],
                tp,
                candidate_mes,
                node_exons,
            )
        )

    whippet = pd.DataFrame(whippet_records)

    # A microexon/timepoint should map to at most one Whippet event.
    if not whippet.empty:
        duplicate_mask = whippet.duplicated(
            subset=["ME", "timepoint"], keep=False
        )
        if duplicate_mask.any():
            dup = whippet.loc[
                duplicate_mask,
                ["ME", "timepoint", "whippet_Gene", "whippet_Node", "whippet_Type"]
            ]
            dup.to_csv(
                outdir / "06_whippet_duplicate_mappings.tsv",
                sep="\t",
                index=False,
            )
            raise ValueError(
                "Multiple Whippet events mapped to the same candidate "
                "microexon/timepoint. See 06_whippet_duplicate_mappings.tsv."
            )

    gene_cols = [
        c for c in ["ME", "gene_id", "gene_name", "gene_biotype"]
        if c in gene.columns
    ]
    step5 = step5.merge(gene[gene_cols], on="ME", how="left")

    result = step5.merge(
        whippet,
        on=["ME", "timepoint"],
        how="left",
        validate="one_to_one",
    )

    result["whippet_mapped"] = result[
        "whippet_delta_activation_minus_ctrl"
    ].notna()

    me_delta = pd.to_numeric(
        result["delta_PSI_activation_minus_ctrl"], errors="coerce"
    )
    wh_delta = pd.to_numeric(
        result["whippet_delta_activation_minus_ctrl"], errors="coerce"
    )

    result["microexonator_direction"] = me_delta.map(sign_label)
    result["whippet_direction"] = wh_delta.map(sign_label)
    result["direction_concordant"] = (
        result["whippet_mapped"]
        & result["microexonator_direction"].isin(["increase", "decrease"])
        & (result["microexonator_direction"] == result["whippet_direction"])
    )
    result["delta_difference_whippet_minus_microexonator"] = wh_delta - me_delta
    result["abs_delta_difference"] = (
        result["delta_difference_whippet_minus_microexonator"].abs()
    )

    output = outdir / "06_whippet_candidate_concordance.tsv"
    result.to_csv(output, sep="\t", index=False)

    summary_rows = []
    for tp in TIMEPOINT_ORDER:
        sub = result.loc[result["timepoint"] == tp].copy()
        mapped = sub.loc[sub["whippet_mapped"]].copy()

        pearson = np.nan
        spearman = np.nan
        if len(mapped) >= 3:
            x = pd.to_numeric(
                mapped["delta_PSI_activation_minus_ctrl"], errors="coerce"
            )
            y = pd.to_numeric(
                mapped["whippet_delta_activation_minus_ctrl"], errors="coerce"
            )
            valid = x.notna() & y.notna()
            if int(valid.sum()) >= 3:
                pearson = x[valid].corr(y[valid], method="pearson")
                spearman = x[valid].rank().corr(
                    y[valid].rank(), method="pearson"
                )

        summary_rows.append({
            "timepoint": tp,
            "step5_candidate_responses": len(sub),
            "mapped_to_whippet": int(sub["whippet_mapped"].sum()),
            "not_mapped_to_whippet": int((~sub["whippet_mapped"]).sum()),
            "direction_concordant": int(
                sub["direction_concordant"].sum()
            ),
            "direction_concordance_fraction_among_mapped": (
                float(mapped["direction_concordant"].mean())
                if len(mapped) else np.nan
            ),
            "mean_abs_delta_difference": (
                float(mapped["abs_delta_difference"].mean())
                if len(mapped) else np.nan
            ),
            "median_abs_delta_difference": (
                float(mapped["abs_delta_difference"].median())
                if len(mapped) else np.nan
            ),
            "pearson_delta_correlation": pearson,
            "spearman_delta_correlation": spearman,
        })

    summary = pd.DataFrame(summary_rows)

    mapped_all = result.loc[result["whippet_mapped"]].copy()
    overall = {
        "timepoint": "ALL",
        "step5_candidate_responses": len(result),
        "mapped_to_whippet": int(result["whippet_mapped"].sum()),
        "not_mapped_to_whippet": int((~result["whippet_mapped"]).sum()),
        "direction_concordant": int(result["direction_concordant"].sum()),
        "direction_concordance_fraction_among_mapped": (
            float(mapped_all["direction_concordant"].mean())
            if len(mapped_all) else np.nan
        ),
        "mean_abs_delta_difference": (
            float(mapped_all["abs_delta_difference"].mean())
            if len(mapped_all) else np.nan
        ),
        "median_abs_delta_difference": (
            float(mapped_all["abs_delta_difference"].median())
            if len(mapped_all) else np.nan
        ),
        "pearson_delta_correlation": np.nan,
        "spearman_delta_correlation": np.nan,
    }

    if len(mapped_all) >= 3:
        x = pd.to_numeric(
            mapped_all["delta_PSI_activation_minus_ctrl"], errors="coerce"
        )
        y = pd.to_numeric(
            mapped_all["whippet_delta_activation_minus_ctrl"], errors="coerce"
        )
        valid = x.notna() & y.notna()
        if int(valid.sum()) >= 3:
            overall["pearson_delta_correlation"] = x[valid].corr(
                y[valid], method="pearson"
            )
            overall["spearman_delta_correlation"] = x[valid].rank().corr(
                y[valid].rank(), method="pearson"
            )

    summary = pd.concat(
        [summary, pd.DataFrame([overall])],
        ignore_index=True,
    )
    summary.to_csv(
        outdir / "06_whippet_concordance_summary.tsv",
        sep="\t",
        index=False,
    )

    unique_candidate_mes = int(step5["ME"].nunique())
    unique_mapped_mes = int(
        result.loc[result["whippet_mapped"], "ME"].nunique()
    )

    report = [
        "WHIPPET VALIDATION - STEP 6",
        "",
        f"Step 5 candidate responses: {len(step5)}",
        f"Unique Step 5 microexons: {unique_candidate_mes}",
        f"Candidate responses mapped to Whippet: {int(result['whippet_mapped'].sum())}",
        f"Unique candidate microexons mapped to Whippet: {unique_mapped_mes}",
        f"Direction-concordant mapped responses: {int(result['direction_concordant'].sum())}",
        "",
        "Delta PSI convention after harmonization:",
        "  MicroExonator = activation - control",
        "  Whippet corrected = -(Whippet raw DeltaPsi)",
        "",
        "Effect-size agreement is reported as the absolute DeltaPSI difference",
        "and as Pearson/Spearman correlations. No arbitrary agreement cutoff",
        "is imposed in this step.",
        "",
        "Whippet Probability is retained as Whippet output and is not treated",
        "as a p-value or FDR.",
        "",
        "Whippet comparison files:",
    ]

    for tp in TIMEPOINT_ORDER:
        report.append(f"  {tp}: {delta_files[tp]}")

    (outdir / "06_whippet_validation_report.txt").write_text(
        "\n".join(report) + "\n",
        encoding="utf-8",
    )
    print("\n".join(report))


if __name__ == "__main__":
    main()
