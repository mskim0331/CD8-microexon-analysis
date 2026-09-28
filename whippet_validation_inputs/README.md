# Whippet validation inputs

This directory contains the archived inputs used for the independent Whippet cross-validation in Step 6.

## Files

- `whippet.jls.exons.tab.gz`: exon/node mapping table generated with the Whippet index used for the paired-end analysis.
- `Delta/`: seven Ctrl-vs-activation Whippet differential-splicing tables for 30 min, 3 h, 12 h, 24 h, 48 h, 72 h, and 7 d.

The differential comparisons use Ctrl as group A and activation as group B. Whippet therefore reports raw `DeltaPsi = Psi_A - Psi_B`. The Step 6 analysis negates this value so that the repository-wide convention is `activation - control`.

## Integrity and provenance

The original SHA-256 checksums from the server bundle are preserved in:

```text
provenance/whippet_step6_SHA256SUMS.txt
```

The Step 6 GitHub Actions workflow independently checks the exon table and all seven differential files against those hashes and tests gzip integrity before analysis.

Whippet source provenance is recorded in:

```text
provenance/whippet_git_head.txt
provenance/whippet_git_status.txt
provenance/whippet_modified_files.txt
provenance/whippet_diff_stat.txt
patches/whippet_local_changes.patch
```

These files document the pinned Whippet commit and the local code changes present in the analysis environment.
