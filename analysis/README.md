# Analysis

The scripts in this directory can be run locally or through GitHub Actions. The GitHub Actions run used Python 3.11.16 with the package versions listed in `requirements.txt`.

## Local environment

From the repository root:

```bash
python3.11 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -r analysis/requirements.txt
```

On macOS, if `python3.11` is not available, install Python 3.11 first or use a Conda environment with Python 3.11.

## Step 1: structural QC and high-confidence PSI matrix

```bash
python analysis/01_build_microexon_qc.py \
  --high-quality out.high_quality.txt \
  --psi out_filtered_ME.PSI.txt \
  --outdir analysis_results
```

This step:

1. verifies the expected columns in both MicroExonator outputs;
2. counts rows, unique microexon coordinates, and biological samples;
3. checks for duplicated sample/microexon records;
4. restricts the PSI table to coordinates present in `out.high_quality.txt`;
5. builds a 656 x 24 PSI matrix;
6. measures missing PSI values for every sample; and
7. calculates pairwise Pearson and Spearman correlations among biological samples.

Main outputs:

```text
analysis_results/01_structure_summary.tsv
analysis_results/01_sample_qc.tsv
analysis_results/01_sample_correlations.tsv
analysis_results/01_high_quality_PSI_matrix.tsv
analysis_results/01_high_quality_annotation.tsv
analysis_results/01_qc_report.txt
```

No activation-responsive microexon is selected in Step 1.

## Step 2: control-referenced time-course summary

```bash
python analysis/02_summarize_timecourse.py \
  --matrix analysis_results/01_high_quality_PSI_matrix.tsv \
  --annotation analysis_results/01_high_quality_annotation.tsv \
  --outdir analysis_results \
  --min-valid-replicates 2 \
  --exploratory-delta 0.10 \
  --strong-delta 0.20
```

For each microexon and condition, the script calculates the number of valid replicates, mean PSI, median PSI, standard deviation, minimum, maximum, and range.

For each activation time point:

```text
Delta PSI = mean PSI(activation) - mean PSI(control)
```

A comparison is considered testable only when at least two of the three replicates have non-missing PSI in both control and the activation condition.

The two effect-size labels are descriptive filters:

- exploratory: `|Delta PSI| >= 0.10`
- strong: `|Delta PSI| >= 0.20`

They are not p-value or FDR thresholds.

`complete_group_separation` is TRUE when the observed PSI ranges of the two groups do not overlap in the direction of the mean effect. For an increase, every valid activation PSI is greater than every valid control PSI; for a decrease, every valid activation PSI is lower than every valid control PSI.

Main outputs:

```text
analysis_results/02_condition_PSI_summary.tsv
analysis_results/02_ctrl_vs_timepoint_deltaPSI.tsv
analysis_results/02_filter_counts_by_timepoint.tsv
analysis_results/02_exploratory_responses_absDelta_ge_0.10.tsv
analysis_results/02_strong_responses_absDelta_ge_0.20.tsv
analysis_results/02_strong_responses_complete_separation.tsv
analysis_results/02_microexon_timecourse_summary.tsv
analysis_results/02_timecourse_report.txt
```

## Step 3: gene annotation

```bash
python analysis/03_annotate_microexons_to_genes.py \
  --high-quality out.high_quality.txt \
  --gtf Homo_sapiens.GRCh38.114.gtf \
  --outdir analysis_results
```

This step maps each high-quality MicroExonator coordinate through its associated Ensembl transcript ID to Ensembl gene ID, gene name, and gene biotype. The repository stores the resulting annotation table and transcript-to-gene map in `analysis_results/`.

Main outputs:

```text
analysis_results/03_microexon_gene_annotation.tsv
analysis_results/03_transcript_gene_map.tsv
analysis_results/03_gene_annotation_report.txt
```

## Step 6: independent validation with paired-end Whippet

```bash
python analysis/06_validate_with_whippet.py \
  --candidates analysis_results/02_strong_responses_complete_separation.tsv \
  --gene-annotation analysis_results/03_microexon_gene_annotation.tsv \
  --whippet-exons whippet_validation_inputs/whippet.jls.exons.tab.gz \
  --delta-root whippet_validation_inputs/Delta \
  --outdir analysis_results
```

The validation uses the seven paired-end Ctrl-vs-activation Whippet differential-splicing files. Candidate Whippet events are converted to MicroExonator coordinates using the exon table generated with the Whippet index. Whippet raw `DeltaPsi` is defined as `Psi_A - Psi_B`; because A is Ctrl and B is activation here, the sign is reversed so that both methods use:

```text
Delta PSI = activation - control
```

If the same genomic microexon maps to multiple Whippet gene models, the ambiguity is resolved only when exactly one Whippet `Gene` matches the Ensembl `gene_id` assigned to that microexon in Step 3. The diagnostic table records all such cases and the selected mapping.

Effect-size agreement is reported using the absolute difference between the two Delta PSI estimates and Pearson/Spearman correlations. No arbitrary effect-agreement threshold is imposed. Whippet `Probability` is retained as a Whippet output and is not treated as a p-value or FDR.

Main outputs:

```text
analysis_results/06_whippet_candidate_concordance.tsv
analysis_results/06_whippet_concordance_summary.tsv
analysis_results/06_whippet_duplicate_mappings.tsv
analysis_results/06_whippet_validation_report.txt
```

The GitHub Actions workflow `.github/workflows/step6-whippet-validation.yml` verifies the archived Step 6 inputs against their recorded SHA-256 checksums before running the validation.
