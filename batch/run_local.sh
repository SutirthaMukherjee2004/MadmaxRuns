#!/bin/bash
set -e
cd "$(dirname "$0")/.."
T=${1:-auto}; OUT=${MADMAX_OUT:-results}; J=${JULIA:-julia}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
$J --project=. -t $T scripts/selftest.jl "$OUT/smoke" | tee selftest.log
grep -q "ALL TESTS PASSED" selftest.log || { echo "self-test failed, stopping"; exit 1; }
for list in batch/jobs_stage1.txt batch/jobs_stage2.txt; do
    grep -v '^\s*#' "$list" | grep -v '^\s*$' | while read -r line; do
        $J --project=. -t $T scripts/run_study.jl $line outdir="$OUT"
    done
done
$J --project=. scripts/merge_tables.jl "$OUT"
$J --project=plots plots/scripts/make_figures.jl "$OUT" "$OUT/figs"
