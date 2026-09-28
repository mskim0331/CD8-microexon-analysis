#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import requests

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


API_BASE = "https://biit.cs.ut.ee/gprofiler/api"
SOURCES = ["GO:BP", "REAC"]


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Functional enrichment for Step 7 using g:Profiler with a custom "
            "background consisting of genes represented among the 656 "
            "high-confidence MicroExonator microexons."
        )
    )
    p.add_argument(
        "--patterns",
        default="analysis_results/07_temporal_pattern_labels.tsv",
    )
    p.add_argument(
        "--gene-annotation",
        default="analysis_results/03_microexon_gene_annotation.tsv",
    )
    p.add_argument("--outdir", default="analysis_results")
    return p.parse_args()


def clean_ids(series):
    vals = []
    for x in series.dropna().astype(str):
        x = x.strip()
        if x and x.upper() != "NA":
            vals.append(x)
    return sorted(set(vals))


def query_gprofiler(name, genes, background):
    payload = {
        "organism": "hsapiens",
        "query": genes,
        "sources": SOURCES,
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "domain_scope": "custom",
        "background": background,
        "no_evidences": False,
        "ordered": False,
        "all_results": False,
    }
    response = requests.post(
        f"{API_BASE}/gost/profile/",
        json=payload,
        headers={"User-Agent": "CD8-microexon-analysis/Step7"},
        timeout=120,
    )
    response.raise_for_status()
    data = response.json()
    return payload, data


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    patterns = pd.read_csv(args.patterns, sep="\t", dtype=str)
    annot = pd.read_csv(args.gene_annotation, sep="\t", dtype=str)

    background = clean_ids(annot["gene_id"])
    if len(background) < 100:
        raise ValueError(
            f"Custom background unexpectedly small: {len(background)} unique genes."
        )

    patterns["peak_delta_num"] = pd.to_numeric(
        patterns["peak_delta_PSI"], errors="coerce"
    )

    query_sets = {
        "all_matched": clean_ids(patterns["gene_id"]),
        "increase": clean_ids(
            patterns.loc[patterns["strong_direction"] == "increase", "gene_id"]
        ),
        "decrease": clean_ids(
            patterns.loc[patterns["strong_direction"] == "decrease", "gene_id"]
        ),
        "early_onset": clean_ids(
            patterns.loc[patterns["onset_phase"] == "early", "gene_id"]
        ),
        "intermediate_onset": clean_ids(
            patterns.loc[patterns["onset_phase"] == "intermediate", "gene_id"]
        ),
        "late_onset": clean_ids(
            patterns.loc[patterns["onset_phase"] == "late", "gene_id"]
        ),
    }

    membership_rows = []
    for group, genes in query_sets.items():
        for gene_id in genes:
            membership_rows.append({"query_set": group, "gene_id": gene_id})
    pd.DataFrame(membership_rows).to_csv(
        outdir / "07_gprofiler_query_sets.tsv",
        sep="\t",
        index=False,
    )

    version_response = requests.get(
        f"{API_BASE}/util/data_versions/",
        params={"organism": "hsapiens"},
        headers={"User-Agent": "CD8-microexon-analysis/Step7"},
        timeout=60,
    )
    version_response.raise_for_status()
    versions = version_response.json()
    (outdir / "07_gprofiler_data_versions.json").write_text(
        json.dumps(versions, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    raw_results = {}
    result_rows = []
    summary_rows = []

    for group, genes in query_sets.items():
        if len(genes) < 3:
            raw_results[group] = {
                "skipped": True,
                "reason": "fewer than 3 unique genes",
                "genes": genes,
            }
            summary_rows.append({
                "query_set": group,
                "n_unique_query_genes": len(genes),
                "n_background_genes": len(background),
                "n_significant_terms": 0,
                "smallest_adjusted_p_value": np.nan,
            })
            continue

        payload, data = query_gprofiler(group, genes, background)
        raw_results[group] = {
            "payload": payload,
            "response": data,
        }

        results = data.get("result", [])
        query_ensgs = (
            data.get("meta", {})
            .get("genes_metadata", {})
            .get("query", {})
            .get("query_1", {})
            .get("ensgs", [])
        )
        gene_symbol_map = (
            annot[["gene_id", "gene_name"]]
            .dropna(subset=["gene_id"])
            .drop_duplicates("gene_id")
            .set_index("gene_id")["gene_name"]
            .to_dict()
        )
        pvals = []
        for item in results:
            pval = item.get("p_value")
            pvals.append(pval)
            intersections = item.get("intersections") or []
            intersection_ids = [
                ensg for ensg, evidence in zip(query_ensgs, intersections)
                if evidence
            ]
            intersection_symbols = [
                str(gene_symbol_map.get(ensg, "") or "")
                for ensg in intersection_ids
            ]
            result_rows.append({
                "query_set": group,
                "source": item.get("source"),
                "term_id": item.get("native"),
                "term_name": item.get("name"),
                "adjusted_p_value": pval,
                "significant": item.get("significant"),
                "query_size": item.get("query_size"),
                "intersection_size": item.get("intersection_size"),
                "term_size": item.get("term_size"),
                "effective_domain_size": item.get("effective_domain_size"),
                "precision": item.get("precision"),
                "recall": item.get("recall"),
                "intersection_gene_ids": ";".join(intersection_ids),
                "intersection_gene_symbols": ";".join(
                    x for x in intersection_symbols if x and x.upper() != "NA"
                ),
            })

        summary_rows.append({
            "query_set": group,
            "n_unique_query_genes": len(genes),
            "n_background_genes": len(background),
            "n_significant_terms": len(results),
            "smallest_adjusted_p_value": min(pvals) if pvals else np.nan,
        })

    (outdir / "07_gprofiler_raw_results.json").write_text(
        json.dumps(raw_results, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    result_df = pd.DataFrame(result_rows)
    if result_df.empty:
        result_df = pd.DataFrame(columns=[
            "query_set", "source", "term_id", "term_name",
            "adjusted_p_value", "significant", "query_size",
            "intersection_size", "term_size", "effective_domain_size",
            "precision", "recall", "intersection_gene_ids",
            "intersection_gene_symbols",
        ])
    result_df.to_csv(
        outdir / "07_gprofiler_enrichment.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(
        outdir / "07_gprofiler_enrichment_summary.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # Compact figure: up to five strongest terms per query set.
    plot_df = result_df.copy()
    plot_df["adjusted_p_value"] = pd.to_numeric(
        plot_df["adjusted_p_value"], errors="coerce"
    )
    plot_df = plot_df.loc[
        plot_df["adjusted_p_value"].notna()
        & (plot_df["adjusted_p_value"] > 0)
    ].copy()
    if not plot_df.empty:
        selected = (
            plot_df.sort_values(["query_set", "adjusted_p_value"])
            .groupby("query_set", as_index=False, group_keys=False)
            .head(5)
            .copy()
        )
        selected["minus_log10_fdr"] = -np.log10(selected["adjusted_p_value"])
        selected["label"] = (
            selected["query_set"]
            + " | "
            + selected["source"]
            + " | "
            + selected["term_name"].astype(str)
        )
        selected = selected.sort_values("minus_log10_fdr", ascending=True)

        fig_h = max(6.0, 0.34 * len(selected))
        fig, ax = plt.subplots(figsize=(12, fig_h))
        y = np.arange(len(selected))
        ax.barh(y, selected["minus_log10_fdr"].to_numpy())
        ax.set_yticks(y)
        ax.set_yticklabels(selected["label"].tolist(), fontsize=7)
        ax.set_xlabel("-log10(FDR-adjusted p-value)")
        ax.set_title(
            "g:Profiler functional enrichment of Whippet-matched microexon genes\n"
            "Custom background: genes represented among 656 high-confidence microexons"
        )
        fig.tight_layout()
        fig.savefig(
            outdir / "07_gprofiler_top_terms.png",
            dpi=300,
            bbox_inches="tight",
        )
        fig.savefig(
            outdir / "07_gprofiler_top_terms.svg",
            bbox_inches="tight",
        )
        plt.close(fig)
    else:
        fig, ax = plt.subplots(figsize=(9, 3))
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "No significant GO:BP or Reactome enrichment\nat FDR < 0.05 with the custom background.",
            ha="center",
            va="center",
            fontsize=12,
        )
        fig.tight_layout()
        fig.savefig(
            outdir / "07_gprofiler_top_terms.png",
            dpi=300,
            bbox_inches="tight",
        )
        fig.savefig(
            outdir / "07_gprofiler_top_terms.svg",
            bbox_inches="tight",
        )
        plt.close(fig)

    report = [
        "STEP 7 FUNCTIONAL ENRICHMENT",
        "",
        f"Unique genes in high-confidence custom background: {len(background)}",
        "Sources queried: GO:BP, Reactome",
        "Multiple-testing correction: FDR",
        "Significance threshold: adjusted p < 0.05",
        "Domain scope: custom",
        "",
        "Rationale for custom background:",
        "  The candidate genes arose from a restricted set of genes represented",
        "  among the 656 high-confidence microexons, rather than from an unbiased",
        "  genome-wide gene screen. Using those represented genes as background",
        "  reduces enrichment bias from the upstream detection/selection process.",
        "",
    ]
    for _, row in summary_df.iterrows():
        report.append(
            f"{row['query_set']}: {int(row['n_unique_query_genes'])} genes, "
            f"{int(row['n_significant_terms'])} significant terms"
        )
    (outdir / "07_functional_enrichment_report.txt").write_text(
        "\n".join(report) + "\n",
        encoding="utf-8",
    )

    print("\n".join(report))


if __name__ == "__main__":
    main()
