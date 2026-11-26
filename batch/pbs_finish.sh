#!/bin/bash
#PBS -N madmax_figs
#PBS -j oe
source "${MADMAX_ROOT:-$PBS_O_WORKDIR}/batch/config.sh"
cd "$MADMAX_ROOT"
export PATH=$JULIA_BIN:$PATH
julia --project="$MADMAX_ROOT" scripts/merge_tables.jl "$MADMAX_OUT" || exit 1
julia --project="$MADMAX_ROOT/plots" plots/scripts/make_figures.jl "$MADMAX_OUT" "$MADMAX_OUT/figs"
