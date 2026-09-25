# CD8 Microexon Analysis

Microexon discovery and paired-end splicing analysis of a public human CD8+ T-cell activation RNA-seq time course.

The analysis uses **MicroExonator** for microexon discovery and microexon PSI estimation, then uses a microexon-augmented **Whippet** index for paired-end event-level splicing quantification and control-versus-activation differential analysis.

## Study design

- Dataset: **GSE247647** (24-sample activation time-course subset of the full 51-sample GEO series)
- Material: human CD8+ T-cell bulk RNA-seq
- Design: control plus seven activation time points
- Time points: 30 min, 3 h, 12 h, 24 h, 48 h, 72 h, and 7 d
- Biological replicates: 3 per condition
- Sequencing: paired-end
- Total biological samples: 24
- Total FASTQ files: 48

## Analysis overview

```text
paired-end FASTQ
        |
        v
MicroExonator discovery
        |
        +--> high-confidence microexon set
        |
        +--> MicroExonator microexon PSI
        |
        v
microexon-augmented GTF
        |
        v
Whippet index
        |
        v
pair-aware removal of read pairs containing ambiguous N bases
        |
        v
paired-end Whippet quantification
(one PSI file per biological sample)
        |
        v
Ctrl n=3 vs activation n=3
        |
        v
7 Whippet differential comparisons
```

## Repository contents

- `protocol.qmd` — complete reproducible protocol
- `previous-analysis-audit.qmd` — optional, separate methodological audit of the earlier analysis
- `scripts/` — executable preprocessing, paired-end quantification, and differential-analysis scripts
- `config/` — portable configuration templates and sample definitions
- `patches/` — minimal compatibility patch applied to the pinned MicroExonator source
- `provenance/` — software, reference, and run provenance retained from the analysis
- `docs/` — rendered Quarto HTML output after `quarto render`

## Analysis note

The paired-end Whippet `.psi.gz` and `.diff.gz` files are **general Whippet splicing-event tables**, not microexon-only tables. Microexon-specific downstream interpretation must therefore intersect or map Whippet events to the high-confidence MicroExonator microexon set.

Whippet `Probability` is not an adjusted p-value or FDR.

## Full protocol

See [protocol.qmd](protocol.qmd). After rendering with Quarto, the HTML version is written to `docs/protocol.html`.

## Status

The `protocol-draft` branch contains the working version of the protocol and associated analysis scripts.
