#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
source batch/config.sh
mkdir -p "$MADMAX_LOGS/pbs"

qsub_bin=${QSUB:-qsub}
common=(-q "$PBS_QUEUE" -o "$MADMAX_LOGS/pbs/")
envvars="MADMAX_ROOT=$MADMAX_ROOT,MADMAX_OUT=$MADMAX_OUT,JULIA=$JULIA,PBS_NCPUS=$PBS_NCPUS,PBS_MEM=$PBS_MEM"
if [[ -n "${JULIA_BIN:-}" ]]; then
    envvars="$envvars,JULIA_BIN=$JULIA_BIN"
fi
resources="select=1:ncpus=$PBS_NCPUS:mem=$PBS_MEM"
base="$MADMAX_ROOT/batch/jobs_paper_base.txt"
scan="$MADMAX_ROOT/batch/jobs_paper_scan.txt"
base_end=$(wc -l < "$base")
scan_end=$(wc -l < "$scan")
limit=${PBS_PAPER_CONCURRENCY:-3}

if (( $# == 0 )); then
    example=$("$qsub_bin" "${common[@]}" -N madmax_fig14 -l "$resources" -l "walltime=$PBS_WALLTIME" \
        -v "$envvars,JOBLIST=$base" batch/pbs_array.sh)
    base_job=$("$qsub_bin" "${common[@]}" -N madmax_base -l "$resources" -l "walltime=$PBS_WALLTIME" \
        -J "2-$base_end%$limit" -v "$envvars,JOBLIST=$base" batch/pbs_array.sh)
    scan_job=
elif [[ "$1" == "--continue" ]] && (( $# == 3 || $# == 4 )); then
    example=$2
    base_job=$3
    scan_job=${4:-}
else
    echo "usage: $0 [--continue figure14_job base_array [scan_array]]" >&2
    exit 2
fi
echo "Figure 14 study: $example"
echo "Independent figure studies: $base_job"

if [[ -z "$scan_job" ]]; then
    prereq=()
    if [[ ! -s "$MADMAX_DATA/tauscan.jls" || ! -s "$MADMAX_DATA/tau0.jls" ]]; then
        prereq=(-W "depend=afterok:$base_job")
    fi
    scan_job=$("$qsub_bin" "${common[@]}" -N madmax_scan -l "$resources" -l "walltime=$PBS_WALLTIME" \
        -J "1-$scan_end%$limit" "${prereq[@]}" -v "$envvars,JOBLIST=$scan" batch/pbs_array.sh)
fi
echo "Dependent scans: $scan_job"

finish=$("$qsub_bin" "${common[@]}" -l "select=1:ncpus=4:mem=$PBS_MEM" -l walltime=02:00:00 \
    -W "depend=afterok:$example:$base_job:$scan_job" -v "$envvars" batch/pbs_paper_finish.sh)
echo "Figure export and notebook: $finish"
