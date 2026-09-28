# Configuration and sample registry

These files preserve the configuration used on the analysis server.

- `local_samples.tsv`: 48 paired FASTQ entries (24 biological samples) with the SRA accession embedded in each original file path.
- `Paired_samples.tsv`: the 24 R1/R2 biological pairs used by MicroExonator/Whippet.
- `config.yaml`: MicroExonator configuration from the completed server run.

The paths in `config.yaml` and `local_samples.tsv` are the original absolute server paths and are kept for provenance. They must be changed before rerunning the workflow on another machine.

Key MicroExonator settings are `ME_len: 30`, `min_number_files_detected: 3`, and `min_reads_PSI: 5`.
