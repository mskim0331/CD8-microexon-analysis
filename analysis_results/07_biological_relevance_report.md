# Step 7 biological relevance interpretation

## Scope

This document interprets the temporal microexon results after independent Whippet matching. It distinguishes three levels of evidence:

1. **Observed in this dataset:** MicroExonator PSI/DeltaPSI trajectories, Step 5 reproducibility criteria, and Whippet direction concordance.
2. **Functional enrichment:** gene-set enrichment using g:Profiler with a custom statistical background.
3. **External biological context:** published evidence about the host gene or pathway.

Published gene function is **not** treated as proof that the detected microexon is itself functional.

## Temporal result

Among the 75 Whippet-matched candidate microexons:

- 27 first met the Step 5 response criteria in the **early** phase (30 min or 3 h).
- 29 first met the criteria in the **intermediate** phase (12 h or 24 h).
- 19 first met the criteria in the **late** phase (48 h, 72 h, or 7 d).
- 26 had an increase as their consistent strong-response direction.
- 47 had a decrease as their consistent strong-response direction.
- 2 showed strong responses in both directions across the time course and are therefore retained as **mixed/direction-switching** rather than forced into an increase/decrease class.
- 42 still met the Step 5 response criteria at 7 d. This is described as **retained at 7 d** and does not imply that every intervening timepoint was strong.

The unsupervised trajectory clustering selected two clusters (mean silhouette 0.321), which mainly separated increasing from decreasing trajectories. Because this clustering did not resolve clearly distinct early/transient/late biological subtypes, timing interpretation is based primarily on the observed first strong-response timepoint rather than on post-hoc cluster names.

## Replicate-coverage sensitivity

The primary Step 5 rule allowed at least 2 of 3 non-missing biological replicates in both control and activation groups. Requiring all 3 of 3 replicates as a sensitivity analysis retained:

- 229 / 279 original Step 5 microexon-timepoint responses (82.1%)
- 65 / 83 original Step 5 unique microexons
- 202 Whippet-mapped strict responses
- 58 Whippet-mapped unique microexons

This indicates that most response-level calls persist under the stricter coverage rule, while a meaningful subset depends on allowing 2/3 valid replicates. The main analysis therefore remains exploratory and reports replicate coverage explicitly.

## Functional enrichment

g:Profiler was run with:

- organism: human
- sources: GO Biological Process and Reactome
- multiple-testing correction: FDR
- threshold: adjusted p < 0.05
- custom background: the 559 unique genes represented by the 656 high-confidence microexons
- g:Profiler database: Ensembl 114 / GRCh38.p14 snapshot recorded in `07_gprofiler_data_versions.json`

No significant terms were found for the full matched set, the increase-only set, the decrease-only set, the early-onset set, or the late-onset set.

The **intermediate-onset** set (28 unique genes) showed 13 significant GO Biological Process terms. The strongest result was:

- **regulation of actin filament bundle assembly** — FDR-adjusted p = 0.02465
  - CLASP1
  - PPFIA1
  - ROCK2
  - FHOD1
  - FLNA

A broader **actin filament organization** term was also enriched (adjusted p = 0.04234) and additionally included ARHGAP12.

A distinct enriched term, **positive regulation of cation transmembrane transport** (adjusted p = 0.04234), contained STIM2, PPP3CC, and FLNA.

These results support a phase-specific association between intermediate microexon responses and cytoskeletal / ion-signaling biology. They do **not** support a claim that all 75 candidate genes are globally enriched for a single T-cell pathway.

## Literature-supported candidates

The file `07_literature_context_candidates.tsv` lists candidates with particularly interpretable host-gene biology. Examples include:

- **STIM2** — T-cell calcium/NFAT signaling; intermediate decrease.
- **MAP4K4** — T-cell activation/MAPK signaling; intermediate increase.
- **FLNA** — T-cell immunological synapse and cytoskeletal signaling; intermediate decrease and part of the enriched actin-related terms.
- **NLRC5** — MHC-I transcriptional control, directly relevant to CD8 antigen-presentation biology; early decrease.
- **HNRNPC** and **SF3B3** — RNA-processing/spliceosome components, providing a plausible connection to regulation of the splicing program itself.
- **NCOR1** — T-cell transcriptional-state control; late increase.
- **ARHGEF1** — RhoA/actin-dependent lymphocyte migration; late decrease.

These candidates are useful for biological interpretation and follow-up, but the current RNA-seq analysis alone cannot establish the protein-level or cellular consequence of any individual microexon.

## Presentation-level interpretation

A defensible summary is:

> Activation-associated microexon responses were temporally heterogeneous, with the largest onset group appearing at 12–24 h. Functional enrichment was not significant across the entire matched set, but the intermediate-onset group was enriched for actin/cytoskeletal organization, with additional literature-supported candidates connected to T-cell calcium signaling, immune regulation, RNA processing, and lymphocyte migration.

