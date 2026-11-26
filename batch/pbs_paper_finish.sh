#!/bin/bash
#PBS -N madmax_paper
#PBS -j oe
set -eo pipefail
source "${MADMAX_ROOT:-$PBS_O_WORKDIR}/batch/config.sh"
cd "$MADMAX_ROOT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$JULIA" --project="$MADMAX_ROOT/plots" plots/scripts/run_paper_notebook.jl
