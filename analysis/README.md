# Downstream analysis

These scripts start from the compact MicroExonator result tables and the archived Whippet differential-splicing files. They do not rerun raw-read alignment or the full MicroExonator workflow.

## Environment

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r analysis/requirements.txt
```

The recorded GitHub Actions environment used Python 3.11.16. Package versions are pinned in `analysis/requirements.txt`.

## Conceptual workflow and scripts

| Analysis stage | Script |
|---|---|
| PSI matrix and sample QC | `01_build_microexon_qc.py` |
| Ctrl-referenced time-course summary and candidate filtering | `02_summarize_timecourse.py` |
| Ensembl gene annotation | `03_annotate_microexons_to_genes.py` |
| Independent Whippet validation | `06_validate_with_whippet.py` |
| Whippet validation figures | `06_plot_whippet_validation.py` |
| Temporal matrix, heatmaps, exploratory clustering | `07_analyze_temporal_patterns.py` |
| Temporal labels, sensitivity analysis, representative trajectories | `07_finalize_temporal_biology.py` |
| GO:BP / Reactome enrichment | `07_functional_enrichment.py` |

The file prefixes reflect the order in which the scripts were developed; they are not meant to define a separate biological step numbering scheme.

## 1. Build the high-confidence PSI matrix and run QC

```bash
python analysis/01_build_microexon_qc.py \
  --high-quality microexonator_outputs/out.high_quality.txt \
  --psi microexonator_outputs/out_filtered_ME.PSI.txt \
  --outdir analysis_results
```

This restricts the long PSI table to the 656 high-confidence MicroExonator coordinates, builds the 656 × 24 PSI matrix, summarizes missingness, and calculates pairwise Pearson/Spearman sample correlations.

## 2. Summarize the activation time course and select candidate responses

```bash
python analysis/02_summarize_timecourse.py \
  --matrix analysis_results/01_high_quality_PSI_matrix.tsv \
  --annotation analysis_results/01_high_quality_annotation.tsv \
  --outdir analysis_results \
  --min-valid-replicates 2 \
  --exploratory-delta 0.10 \
  --strong-delta 0.20
```

For every activation time point, `ΔPSI = mean PSI(activation) - mean PSI(control)`.

The final candidate-response table is `02_strong_responses_complete_separation.tsv`: at least 2/3 valid replicates in both groups, `|ΔPSI| >= 0.20`, and non-overlapping observed PSI ranges in the direction of the effect.

## 3. Map microexons to Ensembl genes

```bash
python analysis/03_annotate_microexons_to_genes.py \
  --annotation analysis_results/01_high_quality_annotation.tsv \
  --gtf Homo_sapiens.GRCh38.114.gtf.gz \
  --outdir analysis_results
```

The script maps MicroExonator-associated transcript IDs to Ensembl gene IDs, gene names, and gene biotypes using Ensembl release 114.

## 4. Validate candidate responses with Whippet

```bash
python analysis/06_validate_with_whippet.py \
  --candidates analysis_results/02_strong_responses_complete_separation.tsv \
  --gene-annotation analysis_results/03_microexon_gene_annotation.tsv \
  --whippet-exons whippet_validation_inputs/whippet.jls.exons.tab.gz \
  --delta-root whippet_validation_inputs/Delta \
  --outdir analysis_results
```

Whippet raw `DeltaPsi` is negated so that both methods use `activation - control`. If one genomic microexon maps to overlapping Whippet gene models, the mapping is retained only when exactly one Whippet gene matches the Step 3 Ensembl `gene_id`. Ambiguous cases are preserved in `06_whippet_duplicate_mappings.tsv`.

Whippet `Probability` is retained as Whippet output and is not treated as a p-value or FDR.

## 5. Analyze temporal patterns

```bash
python analysis/07_analyze_temporal_patterns.py \
  --timecourse analysis_results/02_ctrl_vs_timepoint_deltaPSI.tsv \
  --whippet analysis_results/06_whippet_candidate_concordance.tsv \
  --gene-annotation analysis_results/03_microexon_gene_annotation.tsv \
  --outdir analysis_results
```

This constructs the 75 × 7 ΔPSI matrix for Whippet-matched candidates and generates the temporal heatmaps. Missing values remain missing in the reported matrix and heatmaps; interpolation is used only for exploratory clustering.

```bash
python analysis/07_finalize_temporal_biology.py \
  --timecourse analysis_results/02_ctrl_vs_timepoint_deltaPSI.tsv \
  --step5 analysis_results/02_strong_responses_complete_separation.tsv \
  --whippet analysis_results/06_whippet_candidate_concordance.tsv \
  --gene-annotation analysis_results/03_microexon_gene_annotation.tsv \
  --outdir analysis_results
```

Descriptive onset bins are early (30 min or 3 h), intermediate (12 h or 24 h), and late (48 h, 72 h, or 7 d). The script also reports the stricter 3/3-replicate sensitivity analysis and selects data-driven examples for the PSI trajectory figure.

## 6. Functional enrichment

```bash
python analysis/07_functional_enrichment.py \
  --patterns analysis_results/07_temporal_pattern_labels.tsv \
  --gene-annotation analysis_results/03_microexon_gene_annotation.tsv \
  --outdir analysis_results
```

g:Profiler is queried for GO Biological Process and Reactome terms with FDR < 0.05. The custom statistical background is the set of genes represented among the 656 high-confidence microexons, rather than all human genes. The raw response and database version metadata are saved in `analysis_results/`.

## Key result files

```text
01_high_quality_PSI_matrix.tsv
02_ctrl_vs_timepoint_deltaPSI.tsv
02_strong_responses_complete_separation.tsv
03_microexon_gene_annotation.tsv
06_whippet_candidate_concordance.tsv
06_whippet_concordance_summary.tsv
07_temporal_pattern_labels.tsv
07_replicate_coverage_sensitivity.tsv
07_gprofiler_enrichment.tsv
07_biological_relevance_report.md
```

Presentation-ready figures:

```text
06_whippet_validation_scatter.png
06_whippet_direction_concordance_by_timepoint.png
07_temporal_heatmap_by_onset.png
07_onset_phase_direction_counts.png
07_representative_candidate_PSI_trajectories.png
07_gprofiler_top_terms.png
```
