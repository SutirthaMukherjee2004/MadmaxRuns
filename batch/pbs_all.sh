#!/bin/bash
set -e
cd "$(dirname "$0")/.."
export MADMAX_ROOT=${MADMAX_ROOT:-$PWD}
source batch/config.sh
mkdir -p "$MADMAX_OUT/pbs_logs"
count() { grep -v '^\s*#' "$1" | grep -v '^\s*$' | wc -l; }
RES="select=1:ncpus=$PBS_NCPUS:mem=$PBS_MEM"
COMMON="-q $PBS_QUEUE -o $MADMAX_OUT/pbs_logs/ -V"

T=$(qsub $COMMON -l $RES -l walltime=01:00:00 \
      -v CMD="--project=$MADMAX_ROOT scripts/selftest.jl $MADMAX_OUT/smoke",MADMAX_ROOT="$MADMAX_ROOT" batch/pbs_single.sh)
echo "self-test: $T"
N1=$(count batch/jobs_stage1.txt)
S1=$(qsub $COMMON -l $RES -l walltime=$PBS_WALLTIME -J 1-$N1 -W depend=afterok:$T \
      -v JOBLIST="$MADMAX_ROOT/batch/jobs_stage1.txt",MADMAX_ROOT="$MADMAX_ROOT" batch/pbs_array.sh)
echo "stage 1 ($N1 sub-jobs): $S1"
N2=$(count batch/jobs_stage2.txt)
S2=$(qsub $COMMON -l $RES -l walltime=$PBS_WALLTIME -J 1-$N2 -W depend=afterok:$S1 \
      -v JOBLIST="$MADMAX_ROOT/batch/jobs_stage2.txt",MADMAX_ROOT="$MADMAX_ROOT" batch/pbs_array.sh)
echo "stage 2 ($N2 sub-jobs): $S2"
F=$(qsub $COMMON -l select=1:ncpus=4:mem=$PBS_MEM -l walltime=02:00:00 -W depend=afterok:$S2 \
      -v MADMAX_ROOT="$MADMAX_ROOT" batch/pbs_finish.sh)
echo "merge and figures: $F"
