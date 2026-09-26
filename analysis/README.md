# Analysis

The analysis scripts in this directory are written so that the same calculations can be run locally or in GitHub Actions.

## Step 1: structural QC and PSI matrix

```bash
python -m pip install pandas numpy

python analysis/01_build_microexon_qc.py \
  --high-quality out.high_quality.txt \
  --psi out_filtered_ME.PSI.txt \
  --outdir analysis_results
```

This step does not select activation-responsive candidates. It first checks the structure of the MicroExonator outputs, restricts the PSI table to microexons present in `out.high_quality.txt`, builds the 24-sample PSI matrix, quantifies missing PSI values, and calculates pairwise sample correlations.

Candidate filtering is performed only after these QC results are reviewed.
