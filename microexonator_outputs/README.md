# MicroExonator outputs used downstream

This directory contains the two compact MicroExonator tables used as direct inputs to the Python analysis.

- `out.high_quality.txt`: 656 high-confidence microexons plus header.
- `out_filtered_ME.PSI.txt`: long-format PSI table for the filtered MicroExonator candidates across the 24 biological samples.

The much larger `out.high_quality.gtf` and intermediate MicroExonator working files are not stored in this repository.

The server-side run configuration, source commit, compatibility patch, reference manifest, and software versions are archived under `config/`, `source_snapshot/`, `patches/`, and `provenance/`.
