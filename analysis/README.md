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

The Whippet cross-validation step will be added separately after microexon-specific Whippet differential events are extracted.
