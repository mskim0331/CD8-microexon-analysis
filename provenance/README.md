# Provenance

This directory retains compact records needed to identify the computational environment and source state used for the analysis.

Recommended public files:

- `microexonator_commit.txt` — exact MicroExonator commit
- `system_and_software_versions.txt` — operating system and principal software versions
- `reference_manifest.txt` — reference file names and original file sizes
- `analysis_manifest.txt` — counts and names of principal outputs

Files that expose server-specific working paths or large internal inventories can be removed after the protocol has been finalized.

## Reproducibility caveat

The server snapshot reports Whippet v1.6.2 at git commit `f9cf2e4270af96361670b8b74b13476b3fe18529`, but the Whippet working tree contained local modifications to `src/quant.jl` and `test/test_index.jls`. The exact Whippet diff was not included in the current archive. This must be resolved before claiming a bit-for-bit Whippet source specification.
