#!/bin/bash
#PBS -N madmax
#PBS -j oe
source "${MADMAX_ROOT:-$PBS_O_WORKDIR}/batch/config.sh"
cd "$MADMAX_ROOT"
export PATH=$JULIA_BIN:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=${NCPUS:-$PBS_NCPUS}
LINE=$(grep -v '^\s*#' "$JOBLIST" | grep -v '^\s*$' | sed -n "${PBS_ARRAY_INDEX:-1}p")
[ -z "$LINE" ] && { echo "no job $PBS_ARRAY_INDEX in $JOBLIST"; exit 1; }
echo "sub-job ${PBS_ARRAY_INDEX} (${PBS_JOBID}) on $(hostname): $LINE; $JULIA_NUM_THREADS threads; $(date)"
julia --project="$MADMAX_ROOT" -t "$JULIA_NUM_THREADS" --heap-size-hint=${PBS_MEM%gb}G \
      scripts/run_study.jl $LINE outdir="$MADMAX_OUT"
status=$?
echo "exit $status at $(date)"
exit $status
