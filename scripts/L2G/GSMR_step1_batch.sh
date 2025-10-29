#!/bin/bash
set +e

# ------------------------------------------------------------------------
#  Input
# ------------------------------------------------------------------------
CONFIG=$1
SCRIPT_DIR=`yq .script.path "${CONFIG}"`
BATCH_i=${SLURM_ARRAY_TASK_ID}
BATCH_SIZE=10 # online: 100 * 111

START=$(( (BATCH_i - 1) * BATCH_SIZE + 1 ))
END=$(( BATCH_i * BATCH_SIZE ))
for ((i=START; i<=END; i++))
do
    SLURM_ARRAY_TASK_ID=$i $SCRIPT_DIR/L2G/GSMR_step1.sh "$CONFIG"
done
