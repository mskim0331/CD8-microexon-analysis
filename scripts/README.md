# Server-side Whippet scripts

These scripts document the paired-end Whippet workflow used for the independent validation.

1. `filter_whippet_paired_noN.py` reads R1/R2 together and discards the entire pair if either mate contains `N`; bases are not substituted.
2. `run_whippet_paired.py` runs Whippet with both mates in one command and produces one `.psi.gz` file per biological sample.
3. `run_paired_delta.sh` compares the three Ctrl samples with the three biological replicates at each activation time point, producing seven Ctrl-vs-activation differential files.

The scripts contain the original absolute server paths to Julia and Whippet. They are retained as a record of the executed environment and should be edited for another machine.

The full Whippet index and 24 quantification outputs are not stored here. The exon/node table and seven differential files required for Step 6 are archived in `whippet_validation_inputs/`.
