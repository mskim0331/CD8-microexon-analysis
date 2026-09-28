# CD8 T-cell microexon analysis

Reproducible analysis of microexon regulation during human CD8⁺ T-cell activation using public bulk RNA-seq data from GEO **GSE247647** (BioProject **PRJNA1039860**).

The analysis starts from 24 biological samples (Ctrl and seven activation time points, n=3 per condition), uses **MicroExonator** for microexon discovery and PSI quantification, and uses an independent **Whippet** analysis for cross-validation.

## Dataset

| Condition | Time point | Biological replicates |
|---|---:|---:|
| Control | 0 h | 3 |
| Activated | 30 min | 3 |
| Activated | 3 h | 3 |
| Activated | 12 h | 3 |
| Activated | 24 h | 3 |
| Activated | 48 h | 3 |
| Activated | 72 h | 3 |
| Activated | 7 d | 3 |

SRA accessions are **SRR26818905–SRR26818928**. The exact accession-to-sample mapping used on the server is preserved in `config/local_samples.tsv`.

## Analysis overview

```text
Public paired-end RNA-seq
        |
        v
MicroExonator
  de novo discovery + PSI
        |
        v
656 high-confidence microexons
        |
        v
QC + gene annotation + Ctrl-referenced ΔPSI
        |
        v
83 candidate microexons
279 microexon-timepoint responses
        |
        v
Independent paired-end Whippet analysis
        |
        v
75/83 candidate microexons matched
216/250 matched responses direction-concordant
        |
        v
Temporal-pattern and functional-enrichment analysis
```

Repository-wide ΔPSI is defined as:

```text
ΔPSI = PSI(activation) - PSI(control)
```

## Main results

- **656** high-confidence microexons were retained from the MicroExonator output.
- **83** unique candidate microexons met the reproducible response filter across one or more activation time points.
- **75/83** candidate microexons could be matched to Whippet events.
- **250/279** candidate microexon-timepoint responses were matched to Whippet.
- **216/250 (86.4%)** matched responses showed the same direction of change.
- MicroExonator and Whippet ΔPSI estimates were correlated (**Pearson r = 0.762; Spearman ρ = 0.729**).
- Among the 75 Whippet-matched candidates, first strong responses were classified as **early (27)**, **intermediate (29)**, or **late (19)**.
- Functional enrichment was not significant for the full matched set, early group, or late group. The **intermediate-onset** group showed significant GO Biological Process enrichment, including **regulation of actin filament bundle assembly** (FDR = 0.0247).

Functional enrichment is interpreted at the **host-gene level**. It is not direct evidence that an individual microexon controls the enriched pathway.

## Repository structure

```text
analysis/                    downstream analysis scripts
analysis_results/            result tables, reports, and figures
audit_previous_analysis/     archived material from the earlier exploratory workflow
config/                      sample registries and archived MicroExonator configuration
microexonator_outputs/       MicroExonator tables used by downstream analysis
patches/                     local changes required in the analysis environment
provenance/                  versions, commits, manifests, checksums, and server records
scripts/                     paired-end Whippet preprocessing/quantification scripts
source_snapshot/             archived MicroExonator source files used in the run
whippet_validation_inputs/   exon table and seven Ctrl-vs-activation Whippet diff files
.github/workflows/           reproducible downstream GitHub Actions workflows
```

See `analysis/README.md` for the downstream analysis commands and output definitions.

## Server-side MicroExonator analysis

The MicroExonator run was performed on Ubuntu 22.04.3 and is preserved here through the exact configuration, source snapshot, software commit, compatibility patch, and output tables used downstream.

Key settings:

- MicroExonator commit: `3f8d4aa9c8ace8d1fcc99e8e5cb14f782eccb5c6`
- GRCh38 primary assembly
- Ensembl release 114 GTF
- hg38 phyloP100way conservation track
- `ME_len: 30`
- `min_number_files_detected: 3`
- `min_reads_PSI: 5`
- Snakemake 3.12.0
- Python 3.6.7 in the MicroExonator environment
- Julia 1.6.5

The archived `config/config.yaml` contains the original server paths. These paths are provenance, not portable defaults. Reference files themselves are not stored in the repository.

The required local MicroExonator changes are recorded in `patches/microexonator_required_compatibility.patch`.

The two compact MicroExonator outputs used by the downstream analysis are stored in `microexonator_outputs/`.

## Independent Whippet analysis

Whippet was run separately on paired-end reads, with one PSI file per biological sample. Read pairs were filtered conservatively: if either mate contained an `N`, the entire pair was removed; sequence bases were not altered.

The seven differential comparisons are all **Ctrl vs one activation time point**, with three biological replicates per group.

- Whippet commit: `f9cf2e4270af96361670b8b74b13476b3fe18529`
- Whippet v1.6.2 in completed logs
- Local Whippet modification: `patches/whippet_local_changes.patch`

Whippet reports raw `DeltaPsi = Psi_A - Psi_B`. Because group A was Ctrl and group B was activation, Step 6 reverses the sign before comparison so that both methods use `activation - control`.

The archived Whippet validation inputs and their integrity checks are described in `whippet_validation_inputs/README.md`.

## Candidate filtering

A Ctrl-vs-activation comparison was considered testable when at least 2/3 biological replicates had non-missing PSI in both groups.

The primary candidate filter required:

1. a valid comparison;
2. `|ΔPSI| >= 0.20`; and
3. complete separation of the observed Ctrl and activation PSI ranges.

This is an exploratory effect-size/reproducibility filter, not a p-value threshold. A sensitivity analysis requiring all 3/3 replicates in both groups is reported in the Step 7 outputs.

## Data and large files

The repository does **not** contain raw FASTQ files, BAM/SAM files, the reference FASTA, the full Ensembl GTF, the full MicroExonator high-quality GTF, Whippet quantification files, or the full Whippet index. These files are either publicly retrievable or too large for this repository.

Their identities, versions, server paths, and relevant manifests are recorded under `config/` and `provenance/`.

## Previous exploratory analysis

`audit_previous_analysis/` is retained only to document the earlier workflow and the reasons it was not used as the final method. In particular, the earlier Whippet setup treated R1 and R2 outputs as separate inputs rather than one paired-end biological sample. Nothing in that directory is used to generate the final Step 1–7 results.
