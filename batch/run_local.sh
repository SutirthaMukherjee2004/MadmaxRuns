#!/bin/bash
set -eo pipefail
cd "$(dirname "$0")/.."
source batch/config.sh
T=${1:-auto}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
trap 'echo "batch stopped at line $LINENO, status $?" >&2' ERR
"$JULIA" --project=. -t "$T" scripts/selftest.jl "$MADMAX_OUT/smoke" 2>&1 | tee -a "$MADMAX_LOGS/selftest.log"
for list in batch/jobs_stage1.txt batch/jobs_stage2.txt; do
    grep -v '^\s*#' "$list" | grep -v '^\s*$' | while read -r line; do
        read -r -a args <<< "$line"
        "$JULIA" --project=. -t "$T" scripts/run_study.jl "${args[@]}" output_root="$MADMAX_OUT" resume=true
    done
done
"$JULIA" --project=. scripts/merge_tables.jl "$MADMAX_DATA" 2>&1 | tee -a "$MADMAX_LOGS/tables.log"
"$JULIA" --project=plots plots/scripts/make_figures.jl --require-all
