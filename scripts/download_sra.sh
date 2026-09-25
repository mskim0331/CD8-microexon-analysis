#!/usr/bin/env bash
set -euo pipefail

SAMPLE_SHEET="${1:-config/local_samples.example.tsv}"
OUTPUT_DIR="${2:-data/raw_fastq}"
TMP_DIR="${TMPDIR:-./tmp_fasterq}"

command -v prefetch >/dev/null 2>&1 || {
    echo "ERROR: prefetch not found in PATH" >&2
    exit 1
}
command -v fasterq-dump >/dev/null 2>&1 || {
    echo "ERROR: fasterq-dump not found in PATH" >&2
    exit 1
}
command -v gzip >/dev/null 2>&1 || {
    echo "ERROR: gzip not found in PATH" >&2
    exit 1
}

mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

mapfile -t RUNS < <(
    awk -F '\t' 'NR > 1 {print $1}' "$SAMPLE_SHEET" |
    xargs -n1 basename |
    sed -E 's/_[12]\.fastq\.gz$//' |
    sort -u
)

echo "SRA runs to download: ${#RUNS[@]}"

for run in "${RUNS[@]}"; do
    echo "Downloading $run"

    prefetch "$run"

    fasterq-dump "$run"         --split-files         --threads 8         --temp "$TMP_DIR"         --outdir "$OUTPUT_DIR"

    gzip -f "$OUTPUT_DIR/${run}_1.fastq"
    gzip -f "$OUTPUT_DIR/${run}_2.fastq"
done

echo "SRA download complete."
