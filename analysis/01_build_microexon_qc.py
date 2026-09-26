#!/usr/bin/env python3

import argparse
from itertools import combinations
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


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--high-quality", default="out.high_quality.txt")
    p.add_argument("--psi", default="out_filtered_ME.PSI.txt")
    p.add_argument("--outdir", default="analysis_results")
    return p.parse_args()


def condition_from_sample(sample):
    if sample.startswith("CD8_ctrl_"):
        return "Ctrl"
    for condition in CONDITION_ORDER[1:]:
        if sample.startswith(f"CD8_{condition}_activation_"):
            return condition
    raise ValueError(f"Unrecognized sample name: {sample}")


def replicate_from_sample(sample):
    try:
        return int(sample.rsplit("_", 1)[1])
    except Exception as exc:
        raise ValueError(f"Cannot parse replicate number from: {sample}") from exc


def ordered_samples(samples):
    return sorted(
        samples,
        key=lambda s: (
            CONDITION_ORDER.index(condition_from_sample(s)),
            replicate_from_sample(s),
        ),
    )


def write_key_value_table(path, records):
    pd.DataFrame(records, columns=["metric", "value"]).to_csv(
        path, sep="\t", index=False
    )


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    hq = pd.read_csv(
        args.high_quality,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    psi = pd.read_csv(
        args.psi,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    required_hq = {
        "ME", "Transcript", "ME_length", "U2_score",
        "Mean_conservation", "ME_P_value", "ME_type"
    }
    required_psi = {
        "File", "ME_coords", "PSI", "CI_Lo", "CI_Hi",
        "ME_coverages", "SJ_coverages"
    }

    missing_hq = required_hq - set(hq.columns)
    missing_psi = required_psi - set(psi.columns)
    if missing_hq:
        raise ValueError(f"Missing columns in high-quality table: {sorted(missing_hq)}")
    if missing_psi:
        raise ValueError(f"Missing columns in PSI table: {sorted(missing_psi)}")

    # Convert only the fields needed for numerical summaries.
    hq["ME_length_num"] = pd.to_numeric(hq["ME_length"], errors="coerce")
    hq["U2_score_num"] = pd.to_numeric(hq["U2_score"], errors="coerce")
    hq["Mean_conservation_num"] = pd.to_numeric(
        hq["Mean_conservation"], errors="coerce"
    )
    hq["ME_P_value_num"] = pd.to_numeric(hq["ME_P_value"], errors="coerce")

    psi["PSI_num"] = pd.to_numeric(psi["PSI"], errors="coerce")
    psi["condition"] = psi["File"].map(condition_from_sample)
    psi["replicate"] = psi["File"].map(replicate_from_sample)

    hq_ids = pd.Index(hq["ME"].drop_duplicates())
    hq_id_set = set(hq_ids)

    # Basic structural checks.
    n_hq_rows = len(hq)
    n_hq_unique = hq["ME"].nunique()
    duplicate_hq_rows = int(hq.duplicated(subset=["ME"], keep=False).sum())

    n_psi_rows = len(psi)
    n_psi_unique_me = psi["ME_coords"].nunique()
    samples = ordered_samples(psi["File"].unique().tolist())
    n_samples = len(samples)

    duplicate_psi_pairs = int(
        psi.duplicated(subset=["File", "ME_coords"], keep=False).sum()
    )

    rows_per_sample = psi.groupby("File").size()
    unique_me_per_sample = psi.groupby("File")["ME_coords"].nunique()

    hq_psi = psi.loc[psi["ME_coords"].isin(hq_id_set)].copy()
    hq_seen = hq_psi["ME_coords"].nunique()
    hq_missing_entirely = sorted(hq_id_set - set(hq_psi["ME_coords"]))

    duplicate_hq_sample_pairs = int(
        hq_psi.duplicated(subset=["File", "ME_coords"], keep=False).sum()
    )

    # For the PSI matrix, every sample/microexon pair must be unique.
    if duplicate_hq_sample_pairs:
        raise ValueError(
            f"Found {duplicate_hq_sample_pairs} duplicated high-quality "
            "sample/microexon rows; pivoting would be ambiguous."
        )

    matrix = hq_psi.pivot(
        index="ME_coords",
        columns="File",
        values="PSI_num",
    )
    matrix = matrix.reindex(index=hq_ids, columns=samples)
    matrix.index.name = "ME"

    # Sample-level QC.
    sample_records = []
    for sample in samples:
        all_s = psi.loc[psi["File"] == sample]
        hq_s = hq_psi.loc[hq_psi["File"] == sample]

        sample_records.append({
            "sample": sample,
            "condition": condition_from_sample(sample),
            "replicate": replicate_from_sample(sample),
            "all_candidate_rows": len(all_s),
            "all_candidate_unique_ME": all_s["ME_coords"].nunique(),
            "all_candidate_valid_PSI": int(all_s["PSI_num"].notna().sum()),
            "all_candidate_NA_PSI": int(all_s["PSI_num"].isna().sum()),
            "all_candidate_NA_fraction": float(all_s["PSI_num"].isna().mean()),
            "high_quality_rows": len(hq_s),
            "high_quality_unique_ME": hq_s["ME_coords"].nunique(),
            "high_quality_valid_PSI": int(hq_s["PSI_num"].notna().sum()),
            "high_quality_NA_PSI": int(hq_s["PSI_num"].isna().sum()),
            "high_quality_NA_fraction": float(hq_s["PSI_num"].isna().mean()),
        })

    sample_qc = pd.DataFrame(sample_records)

    # Pairwise sample correlations based only on high-confidence microexons
    # with valid PSI in both samples.
    corr_records = []
    for s1, s2 in combinations(samples, 2):
        x = matrix[s1]
        y = matrix[s2]
        valid = x.notna() & y.notna()
        n_complete = int(valid.sum())

        pearson = np.nan
        spearman = np.nan
        if n_complete >= 3:
            pearson = x[valid].corr(y[valid], method="pearson")
            spearman = x[valid].corr(y[valid], method="spearman")

        corr_records.append({
            "sample_1": s1,
            "sample_2": s2,
            "condition_1": condition_from_sample(s1),
            "condition_2": condition_from_sample(s2),
            "same_condition": condition_from_sample(s1) == condition_from_sample(s2),
            "n_complete_microexons": n_complete,
            "pearson_r": pearson,
            "spearman_rho": spearman,
        })

    correlations = pd.DataFrame(corr_records)

    # High-confidence annotation summaries.
    type_counts = (
        hq.groupby("ME_type", dropna=False)
        .size()
        .rename("n_rows")
        .reset_index()
        .sort_values("ME_type")
    )

    length_counts = (
        hq.groupby("ME_length_num", dropna=False)
        .size()
        .rename("n_rows")
        .reset_index()
        .rename(columns={"ME_length_num": "ME_length"})
        .sort_values("ME_length")
    )

    duplicate_me_table = (
        hq.groupby("ME", as_index=False)
        .agg(
            n_rows=("ME", "size"),
            n_transcripts=("Transcript", "nunique"),
            transcripts=("Transcript", lambda x: ";".join(sorted(set(x)))),
        )
    )
    duplicate_me_table = duplicate_me_table.loc[
        duplicate_me_table["n_rows"] > 1
    ].sort_values(["n_rows", "ME"], ascending=[False, True])

    # A compact annotation table with one row per unique coordinate.
    unique_annotation = (
        hq.groupby("ME", as_index=False)
        .agg(
            transcripts=("Transcript", lambda x: ";".join(sorted(set(x)))),
            n_transcripts=("Transcript", "nunique"),
            ME_length=("ME_length_num", "first"),
            U2_score=("U2_score_num", "first"),
            Mean_conservation=("Mean_conservation_num", "first"),
            ME_P_value=("ME_P_value_num", "first"),
            ME_type=("ME_type", lambda x: ";".join(sorted(set(x)))),
        )
        .set_index("ME")
        .reindex(hq_ids)
        .reset_index()
    )

    # Save analysis tables.
    sample_qc.to_csv(outdir / "01_sample_qc.tsv", sep="\t", index=False)
    correlations.to_csv(
        outdir / "01_sample_correlations.tsv", sep="\t", index=False
    )
    matrix.to_csv(outdir / "01_high_quality_PSI_matrix.tsv", sep="\t")
    unique_annotation.to_csv(
        outdir / "01_high_quality_annotation.tsv", sep="\t", index=False
    )
    type_counts.to_csv(
        outdir / "01_high_quality_type_counts.tsv", sep="\t", index=False
    )
    length_counts.to_csv(
        outdir / "01_high_quality_length_distribution.tsv",
        sep="\t",
        index=False,
    )
    duplicate_me_table.to_csv(
        outdir / "01_duplicate_high_quality_coordinates.tsv",
        sep="\t",
        index=False,
    )

    summary_records = [
        ("high_quality_table_rows", n_hq_rows),
        ("high_quality_unique_ME", n_hq_unique),
        ("high_quality_duplicate_rows_flagged", duplicate_hq_rows),
        ("PSI_table_rows", n_psi_rows),
        ("PSI_unique_candidate_ME", n_psi_unique_me),
        ("PSI_sample_count", n_samples),
        ("PSI_duplicate_sample_ME_rows_flagged", duplicate_psi_pairs),
        ("high_quality_ME_found_in_PSI", hq_seen),
        ("high_quality_ME_missing_entirely_from_PSI", len(hq_missing_entirely)),
        ("expected_high_quality_grid_rows", n_hq_unique * n_samples),
        ("observed_high_quality_grid_rows", len(hq_psi)),
        ("matrix_rows", matrix.shape[0]),
        ("matrix_columns", matrix.shape[1]),
        ("matrix_valid_PSI_values", int(matrix.notna().sum().sum())),
        ("matrix_NA_PSI_values", int(matrix.isna().sum().sum())),
        ("matrix_NA_fraction", float(matrix.isna().mean().mean())),
        ("minimum_ME_length", hq["ME_length_num"].min()),
        ("maximum_ME_length", hq["ME_length_num"].max()),
    ]
    write_key_value_table(outdir / "01_structure_summary.tsv", summary_records)

    # Human-readable report.
    same_cond = correlations.loc[correlations["same_condition"]].copy()
    report = []
    report.append("MICROEXON DATA QC - STEP 1")
    report.append("")
    report.append(f"High-quality table rows: {n_hq_rows}")
    report.append(f"Unique high-quality microexon coordinates: {n_hq_unique}")
    report.append(f"PSI table rows: {n_psi_rows}")
    report.append(f"Unique quantified candidate microexons: {n_psi_unique_me}")
    report.append(f"Biological samples: {n_samples}")
    report.append(f"High-quality microexons present in PSI table: {hq_seen}")
    report.append(
        f"High-quality microexons absent from PSI table: {len(hq_missing_entirely)}"
    )
    report.append(
        f"High-quality PSI matrix: {matrix.shape[0]} microexons x "
        f"{matrix.shape[1]} samples"
    )
    report.append(
        f"Overall high-quality PSI missing fraction: "
        f"{matrix.isna().mean().mean():.4f}"
    )
    report.append("")
    report.append("Rows per sample (all candidates):")
    for s in samples:
        report.append(
            f"  {s}: rows={rows_per_sample[s]}, "
            f"unique_ME={unique_me_per_sample[s]}"
        )

    if not same_cond.empty:
        report.append("")
        report.append("Within-condition replicate correlations:")
        for _, row in same_cond.iterrows():
            report.append(
                f"  {row['sample_1']} vs {row['sample_2']}: "
                f"n={int(row['n_complete_microexons'])}, "
                f"Pearson={row['pearson_r']:.4f}, "
                f"Spearman={row['spearman_rho']:.4f}"
            )

    if hq_missing_entirely:
        report.append("")
        report.append("High-quality microexons absent from PSI table:")
        report.extend(f"  {x}" for x in hq_missing_entirely)

    (outdir / "01_qc_report.txt").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )

    print("\n".join(report))


if __name__ == "__main__":
    main()
