#!/usr/bin/env bash

cd /home/pth1612/CD8_microexon_reproduction/03_MicroExonator_official || exit 1

echo "Waiting for paired FASTQ filtering..."

while true
do
    DONE_PAIRS=$(grep -c '^DONE ' paired_filter.log 2>/dev/null || true)

    FASTQ_COUNT=$(find FASTQ_whippet_paired_clean \
        -name "*.fastq.gz" 2>/dev/null | wc -l)

    echo "$(date) filtering: ${DONE_PAIRS}/24 pairs, ${FASTQ_COUNT}/48 FASTQs"

    if [ "$DONE_PAIRS" -eq 24 ] && [ "$FASTQ_COUNT" -eq 48 ]; then
        echo "Paired filtering complete."
        break
    fi

    if ! pgrep -f "filter_whippet_paired_noN.py" >/dev/null; then
        echo "ERROR: filtering process stopped before completion."
        exit 1
    fi

    sleep 60
done

echo "Starting paired-end Whippet quantification..."

python3 run_whippet_paired.py 2>&1 | tee whippet_paired_run.log

echo "Paired-end Whippet quantification finished."
