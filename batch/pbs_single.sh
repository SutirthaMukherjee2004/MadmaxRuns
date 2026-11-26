#!/bin/bash
#PBS -N madmax1
#PBS -j oe
set -eo pipefail
source "${MADMAX_ROOT:-$PBS_O_WORKDIR}/batch/config.sh"
cd "$MADMAX_ROOT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=${NCPUS:-$PBS_NCPUS}
args=(--project="$MADMAX_ROOT" scripts/selftest.jl "$MADMAX_OUT/smoke")
if [[ -n "${CMD:-}" ]]; then
    read -r -a args <<< "$CMD"
fi
echo "Job on $(hostname), $JULIA_NUM_THREADS threads, $(date): ${args[*]}"
status=0
"$JULIA" -t "$JULIA_NUM_THREADS" "${args[@]}" \
    2>&1 | tee -a "$MADMAX_LOGS/selftest.log" || status=$?
echo "exit $status at $(date)"
exit $status
