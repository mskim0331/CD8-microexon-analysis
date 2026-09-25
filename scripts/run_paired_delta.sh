#!/usr/bin/env bash
set -euo pipefail

: "${JULIA:?Set JULIA to the Julia executable}"
: "${WHIPPET_DELTA:?Set WHIPPET_DELTA to whippet-delta.jl}"

Q="${QUANT_DIR:-results/Whippet_Paired/Quant}"
O="${OUTPUT_DIR:-results/Whippet_Paired/Delta}"

mkdir -p "$O"

CTRL="$Q/CD8_ctrl_1.psi.gz,$Q/CD8_ctrl_2.psi.gz,$Q/CD8_ctrl_3.psi.gz"

run_delta () {
    local name="$1"
    local group_b="$2"

    echo "Running: $name"

    "$JULIA" "$WHIPPET_DELTA"         -a "$CTRL"         -b "$group_b"         -o "$O/$name"
}

run_delta     "CD8_Ctrl_vs_CD8_30min_activation"     "$Q/CD8_30min_activation_1.psi.gz,$Q/CD8_30min_activation_2.psi.gz,$Q/CD8_30min_activation_3.psi.gz"

run_delta     "CD8_Ctrl_vs_CD8_3hr_activation"     "$Q/CD8_3hr_activation_1.psi.gz,$Q/CD8_3hr_activation_2.psi.gz,$Q/CD8_3hr_activation_3.psi.gz"

run_delta     "CD8_Ctrl_vs_CD8_12hr_activation"     "$Q/CD8_12hr_activation_1.psi.gz,$Q/CD8_12hr_activation_2.psi.gz,$Q/CD8_12hr_activation_3.psi.gz"

run_delta     "CD8_Ctrl_vs_CD8_24hr_activation"     "$Q/CD8_24hr_activation_1.psi.gz,$Q/CD8_24hr_activation_2.psi.gz,$Q/CD8_24hr_activation_3.psi.gz"

run_delta     "CD8_Ctrl_vs_CD8_48hr_activation"     "$Q/CD8_48hr_activation_1.psi.gz,$Q/CD8_48hr_activation_2.psi.gz,$Q/CD8_48hr_activation_3.psi.gz"

run_delta     "CD8_Ctrl_vs_CD8_72hr_activation"     "$Q/CD8_72hr_activation_1.psi.gz,$Q/CD8_72hr_activation_2.psi.gz,$Q/CD8_72hr_activation_3.psi.gz"

run_delta     "CD8_Ctrl_vs_CD8_7day_activation"     "$Q/CD8_7day_activation_1.psi.gz,$Q/CD8_7day_activation_2.psi.gz,$Q/CD8_7day_activation_3.psi.gz"

echo "ALL 7 PAIRED-END DIFFERENTIAL COMPARISONS COMPLETE"
