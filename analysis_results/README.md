# Analysis results

This directory contains the tracked outputs from the downstream analysis.

## Core tables

| File | Description |
|---|---|
| `01_high_quality_PSI_matrix.tsv` | 656 high-confidence microexons × 24 biological samples |
| `01_sample_qc.tsv` | sample-level PSI missingness and QC |
| `01_sample_correlations.tsv` | pairwise Pearson/Spearman sample correlations |
| `02_ctrl_vs_timepoint_deltaPSI.tsv` | Ctrl-referenced ΔPSI for all high-confidence microexons |
| `02_strong_responses_complete_separation.tsv` | final response-level candidate table |
| `03_microexon_gene_annotation.tsv` | microexon-to-gene annotation |
| `06_whippet_candidate_concordance.tsv` | MicroExonator/Whippet candidate-level comparison |
| `06_whippet_concordance_summary.tsv` | Whippet validation summary |
| `07_temporal_pattern_labels.tsv` | temporal onset/direction labels for 75 matched candidates |
| `07_gprofiler_enrichment.tsv` | GO:BP/Reactome enrichment output |

## Main figures

- `06_whippet_validation_scatter.png` — MicroExonator vs Whippet ΔPSI
- `06_whippet_direction_concordance_by_timepoint.png` — direction concordance by activation time point
- `07_temporal_heatmap_by_onset.png` — 75 × 7 ΔPSI heatmap ordered by first strong response
- `07_onset_phase_direction_counts.png` — early/intermediate/late onset counts by direction
- `07_representative_candidate_PSI_trajectories.png` — representative PSI trajectories
- `07_gprofiler_top_terms.png` — significant functional-enrichment terms

SVG versions are retained alongside the PNG files for editing and publication-quality export.

## Key counts

```text
High-confidence microexons                    656
Unique candidate microexons                    83
Candidate microexon-timepoint responses       279
Whippet-matched candidate microexons           75
Whippet-matched responses                     250
Direction-concordant matched responses        216 (86.4%)
```

The detailed interpretation and limitations for Step 7 are in `07_biological_relevance_report.md`.
