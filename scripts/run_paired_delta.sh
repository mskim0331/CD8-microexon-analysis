#!/usr/bin/env bash
set -euo pipefail

JULIA="/home/pth1612/anaconda3/envs/smakemake2/bin/julia"
DELTA="/home/pth1612/Whippet.jl/bin/whippet-delta.jl"

Q="Whippet_Paired/Quant"
O="Whippet_Paired/Delta"

CTRL="${Q}/CD8_ctrl_1.psi.gz,${Q}/CD8_ctrl_2.psi.gz,${Q}/CD8_ctrl_3.psi.gz"

run_delta () {
    NAME="$1"
    B="$2"

    echo "======================================"
    echo "Running: ${NAME}"
    echo "======================================"

    "${JULIA}" "${DELTA}" \
        -a "${CTRL}" \
        -b "${B}" \
        -o "${O}/${NAME}"
}

run_delta \
"CD8_Ctrl_vs_CD8_30min_activation" \
"${Q}/CD8_30min_activation_1.psi.gz,${Q}/CD8_30min_activation_2.psi.gz,${Q}/CD8_30min_activation_3.psi.gz"

run_delta \
"CD8_Ctrl_vs_CD8_3hr_activation" \
"${Q}/CD8_3hr_activation_1.psi.gz,${Q}/CD8_3hr_activation_2.psi.gz,${Q}/CD8_3hr_activation_3.psi.gz"

run_delta \
"CD8_Ctrl_vs_CD8_12hr_activation" \
"${Q}/CD8_12hr_activation_1.psi.gz,${Q}/CD8_12hr_activation_2.psi.gz,${Q}/CD8_12hr_activation_3.psi.gz"

run_delta \
"CD8_Ctrl_vs_CD8_24hr_activation" \
"${Q}/CD8_24hr_activation_1.psi.gz,${Q}/CD8_24hr_activation_2.psi.gz,${Q}/CD8_24hr_activation_3.psi.gz"

run_delta \
"CD8_Ctrl_vs_CD8_48hr_activation" \
"${Q}/CD8_48hr_activation_1.psi.gz,${Q}/CD8_48hr_activation_2.psi.gz,${Q}/CD8_48hr_activation_3.psi.gz"

run_delta \
"CD8_Ctrl_vs_CD8_72hr_activation" \
"${Q}/CD8_72hr_activation_1.psi.gz,${Q}/CD8_72hr_activation_2.psi.gz,${Q}/CD8_72hr_activation_3.psi.gz"

run_delta \
"CD8_Ctrl_vs_CD8_7day_activation" \
"${Q}/CD8_7day_activation_1.psi.gz,${Q}/CD8_7day_activation_2.psi.gz,${Q}/CD8_7day_activation_3.psi.gz"

echo "ALL 7 PAIRED-END DIFFERENTIAL COMPARISONS COMPLETE"
