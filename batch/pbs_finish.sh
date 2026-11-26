#!/bin/bash
#PBS -N madmax_figs
#PBS -j oe
set -eo pipefail
source "${MADMAX_ROOT:-$PBS_O_WORKDIR}/batch/config.sh"
cd "$MADMAX_ROOT"
trap 'echo "final stage stopped at line $LINENO, status $?" >&2' ERR
"$JULIA" --project="$MADMAX_ROOT" scripts/merge_tables.jl "$MADMAX_DATA" 2>&1 | tee -a "$MADMAX_LOGS/tables.log"
"$JULIA" --project="$MADMAX_ROOT/plots" plots/scripts/make_figures.jl --require-all
