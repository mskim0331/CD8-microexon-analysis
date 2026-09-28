#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Map MicroExonator transcript IDs to Ensembl gene IDs and gene names."
    )
    p.add_argument(
        "--annotation",
        default="analysis_results/01_high_quality_annotation.tsv",
        help="High-confidence microexon annotation table from Step 1.",
    )
    p.add_argument(
        "--gtf",
        required=True,
        help="Ensembl GTF file (.gtf or .gtf.gz).",
    )
    p.add_argument(
        "--outdir",
        default="analysis_results",
    )
    return p.parse_args()


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def parse_gtf_attributes(text):
    attrs = {}
    for field in text.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        parts = field.split(" ", 1)
        if len(parts) != 2:
            continue
        key, value = parts
        attrs[key] = value.strip().strip('"')
    return attrs


def split_transcripts(value):
    return [x for x in str(value).split(";") if x]


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    me = pd.read_csv(args.annotation, sep="\t", dtype=str, keep_default_na=False)

    required = {"ME", "transcripts"}
    missing = required - set(me.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    wanted_transcripts = {
        transcript
        for value in me["transcripts"]
        for transcript in split_transcripts(value)
    }

    transcript_map = {}

    with open_text(args.gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            feature = fields[2]
            if feature != "transcript":
                continue

            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id", "")

            if transcript_id not in wanted_transcripts:
                continue

            transcript_map[transcript_id] = {
                "transcript_id": transcript_id,
                "gene_id": attrs.get("gene_id", ""),
                "gene_name": attrs.get("gene_name", ""),
                "gene_biotype": attrs.get(
                    "gene_biotype",
                    attrs.get("gene_type", ""),
                ),
                "transcript_biotype": attrs.get(
                    "transcript_biotype",
                    attrs.get("transcript_type", ""),
                ),
                "chromosome": fields[0],
                "strand": fields[6],
            }

    missing_transcripts = sorted(wanted_transcripts - set(transcript_map))

    map_rows = []
    for transcript in sorted(wanted_transcripts):
        record = transcript_map.get(
            transcript,
            {
                "transcript_id": transcript,
                "gene_id": "",
                "gene_name": "",
                "gene_biotype": "",
                "transcript_biotype": "",
                "chromosome": "",
                "strand": "",
            },
        )
        map_rows.append(record)

    transcript_table = pd.DataFrame(map_rows)

    def annotate_row(row):
        transcripts = split_transcripts(row["transcripts"])
        records = [
            transcript_map[t]
            for t in transcripts
            if t in transcript_map
        ]

        def unique_join(key):
            values = sorted({r[key] for r in records if r.get(key, "")})
            return ";".join(values)

        return pd.Series({
            "gene_id": unique_join("gene_id"),
            "gene_name": unique_join("gene_name"),
            "gene_biotype": unique_join("gene_biotype"),
            "n_mapped_transcripts": len(records),
            "n_genes": len({
                r["gene_id"] for r in records if r.get("gene_id", "")
            }),
        })

    mapped = me.apply(annotate_row, axis=1)

    insert_at = list(me.columns).index("transcripts") + 1
    result = me.copy()
    for column in reversed([
        "gene_id",
        "gene_name",
        "gene_biotype",
        "n_mapped_transcripts",
        "n_genes",
    ]):
        result.insert(insert_at, column, mapped[column])

    output_path = outdir / "03_microexon_gene_annotation.tsv"
    transcript_path = outdir / "03_transcript_gene_map.tsv"
    report_path = outdir / "03_gene_annotation_report.txt"

    result.to_csv(output_path, sep="\t", index=False)
    transcript_table.to_csv(transcript_path, sep="\t", index=False)

    n_me = len(result)
    n_me_with_gene = int(result["gene_id"].ne("").sum())
    n_me_with_name = int(result["gene_name"].ne("").sum())
    n_multigene = int((pd.to_numeric(result["n_genes"]) > 1).sum())

    report = [
        "MICROEXON GENE ANNOTATION - STEP 3",
        "",
        f"Input microexons: {n_me}",
        f"Unique transcript IDs requested: {len(wanted_transcripts)}",
        f"Transcript IDs mapped in GTF: {len(transcript_map)}",
        f"Transcript IDs missing from GTF: {len(missing_transcripts)}",
        f"Microexons with gene_id: {n_me_with_gene}",
        f"Microexons with gene_name: {n_me_with_name}",
        f"Microexons mapping to more than one gene: {n_multigene}",
        "",
        "Mapping logic:",
        "  microexon -> MicroExonator transcript ID -> Ensembl GTF transcript record",
        "  -> gene_id and gene_name",
    ]

    if missing_transcripts:
        report.extend(["", "Missing transcript IDs:"])
        report.extend(f"  {t}" for t in missing_transcripts)

    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))


if __name__ == "__main__":
    main()
