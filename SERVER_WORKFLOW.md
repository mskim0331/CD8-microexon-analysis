# Server workflow

This note records the upstream processing that produced the MicroExonator and Whippet files used by the downstream analysis.

## Input data

- GEO: `GSE247647`
- BioProject: `PRJNA1039860`
- SRA: `SRR26818905`–`SRR26818928`
- 24 biological samples
- paired-end RNA-seq
- Ctrl and seven activation time points, n=3 per condition

The original FASTQ-to-sample mapping is in `config/local_samples.tsv`. R1/R2 biological pairs are in `config/Paired_samples.tsv`.

Raw FASTQ files are not stored in this repository.

## MicroExonator

The completed server run used:

- MicroExonator commit `3f8d4aa9c8ace8d1fcc99e8e5cb14f782eccb5c6`
- Snakemake 3.12.0
- Python 3.6.7 in the MicroExonator environment
- Julia 1.6.5
- GRCh38 primary assembly
- Ensembl release 114 annotation
- hg38 phyloP100way conservation track

The exact completed-run configuration is archived in `config/config.yaml`. Its important analysis parameters are:

```yaml
ME_len: 30
min_number_files_detected: 3
min_reads_PSI: 5
```

The source tree used for the run is represented by `source_snapshot/`, and the required local compatibility changes are recorded in `patches/microexonator_required_compatibility.patch`.

Two compact outputs are retained because they are the direct inputs to the downstream analysis:

```text
microexonator_outputs/out.high_quality.txt
microexonator_outputs/out_filtered_ME.PSI.txt
```

The completed server run contained 656 high-confidence microexons. The full `out.high_quality.gtf` and large intermediate files are not included.

The original shell history was not part of the archived bundle, so this repository does not claim an exact byte-for-byte reconstruction of the command line used to launch Snakemake. The pinned source, run configuration, software versions, compatibility patch, sample registry, and resulting analysis inputs are preserved.

## Whippet paired-end validation

Whippet was run independently from the MicroExonator PSI analysis.

Recorded source state:

- Whippet commit `f9cf2e4270af96361670b8b74b13476b3fe18529`
- Whippet v1.6.2 in completed logs
- local change recorded in `patches/whippet_local_changes.patch`

The local source tree was not pristine; this is intentional and documented rather than hidden.

### 1. Pair-aware N filtering

`scripts/filter_whippet_paired_noN.py` reads both mates together. If either mate contains an `N`, the complete pair is discarded. Bases are not substituted.

### 2. Paired quantification

`scripts/run_whippet_paired.py` runs both mates in the same Whippet quantification command and produces one PSI file per biological sample.

The server completed 24 paired-end quantifications.

### 3. Differential splicing

`scripts/run_paired_delta.sh` performs seven comparisons:

```text
Ctrl vs 30 min
Ctrl vs 3 h
Ctrl vs 12 h
Ctrl vs 24 h
Ctrl vs 48 h
Ctrl vs 72 h
Ctrl vs 7 d
```

Each side contains three biological samples.

The seven `.diff.gz` files and the exon/node table required for Step 6 are retained under `whippet_validation_inputs/`. Their SHA-256 hashes are recorded in `provenance/whippet_step6_SHA256SUMS.txt`.

Whippet reports `DeltaPsi = Psi_A - Psi_B`. Here A is Ctrl and B is activation, so downstream Step 6 negates Whippet's raw value to use the project-wide convention:

```text
ΔPSI = activation - control
```

## Environment records

See `provenance/` for:

- system and software versions
- reference-file manifest
- MicroExonator and Whippet commit information
- Whippet source status and modifications
- server-side output counts
- checksums and file manifests

These files are retained so the server-side processing and the downstream analysis can be audited separately.
