#!/bin/bash
#PBS -N madmax
#PBS -j oe
set -eo pipefail
source "${MADMAX_ROOT:-$PBS_O_WORKDIR}/batch/config.sh"
cd "$MADMAX_ROOT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=${NCPUS:-$PBS_NCPUS}
INDEX=${PBS_ARRAY_INDEX:-${MADMAX_JOB_INDEX:-1}}
LINE=$(grep -v '^\s*#' "$JOBLIST" | grep -v '^\s*$' | sed -n "${INDEX}p")
[ -z "$LINE" ] && { echo "no job $INDEX in $JOBLIST"; exit 1; }
echo "sub-job $INDEX (${PBS_JOBID}) on $(hostname): $LINE; $JULIA_NUM_THREADS threads; $(date)"
status=0
read -r -a args <<< "$LINE"
"$JULIA" --project="$MADMAX_ROOT" -t "$JULIA_NUM_THREADS" --heap-size-hint="${PBS_MEM%gb}G" \
      scripts/run_study.jl "${args[@]}" output_root="$MADMAX_OUT" resume=true || status=$?
echo "exit $status at $(date)"
exit $status
