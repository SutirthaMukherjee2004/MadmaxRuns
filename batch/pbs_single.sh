#!/bin/bash
#PBS -N madmax1
#PBS -j oe
source "${MADMAX_ROOT:-$PBS_O_WORKDIR}/batch/config.sh"
cd "$MADMAX_ROOT"
export PATH=$JULIA_BIN:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=${NCPUS:-$PBS_NCPUS}
echo "$CMD on $(hostname), $JULIA_NUM_THREADS threads, $(date)"
julia -t "$JULIA_NUM_THREADS" $CMD
status=$?
echo "exit $status at $(date)"
exit $status
