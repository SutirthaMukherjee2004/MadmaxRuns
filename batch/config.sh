export MADMAX_ROOT=${MADMAX_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
export MADMAX_OUT=${MADMAX_OUT:-$MADMAX_ROOT/results}
[[ "$MADMAX_OUT" = /* ]] || export MADMAX_OUT="$MADMAX_ROOT/$MADMAX_OUT"
export MADMAX_LOGS="$MADMAX_OUT/logs"
export MADMAX_DATA="$MADMAX_OUT/data"
export MADMAX_PLOTS="$MADMAX_OUT/paperplots"
if [[ -n "${JULIA_BIN:-}" ]]; then
    export PATH="$JULIA_BIN:$PATH"
fi
export JULIA=${JULIA:-julia}
export PBS_QUEUE=${PBS_QUEUE:-workq}
export PBS_NCPUS=${PBS_NCPUS:-32}
export PBS_MEM=${PBS_MEM:-16gb}
export PBS_WALLTIME=${PBS_WALLTIME:-24:00:00}
mkdir -p "$MADMAX_LOGS" "$MADMAX_DATA" "$MADMAX_PLOTS"
