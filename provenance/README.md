# Provenance

This directory retains compact, public-facing records needed to identify the computational state used for the analysis.

- `microexonator_commit.txt` — exact MicroExonator source commit
- `software_versions.txt` — principal operating-system and software versions
- `reference_resources.txt` — reference filenames, release identifiers, and recorded file sizes
- `analysis_summary.txt` — final input/output counts and comparison design

The original server archive contained additional machine-specific manifests, absolute paths, logs, and temporary source snapshots. Those materials were used to prepare the protocol but are not required in the final public repository.

## Whippet source-state caveat

The archived server state reported Whippet v1.6.2 at commit
`f9cf2e4270af96361670b8b74b13476b3fe18529`, but the Whippet working tree
contained local modifications to `src/quant.jl` and `test/test_index.jls`.
The exact diff was not captured in the uploaded archive. The final release
should either recover that diff or state this limitation explicitly.
