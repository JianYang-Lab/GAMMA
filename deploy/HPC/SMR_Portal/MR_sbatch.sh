#!/usr/bin/env bash

CONFIG=$1

trait_name=`yq .input.trait "${CONFIG}"`
OUTPUT=`yq .input.output "${CONFIG}"`
SCRIPT_DIR=`yq .script.path "${CONFIG}"`
WORK_DIR=`yq .script.work_path "${CONFIG}"`


cd ${WORK_DIR}/L2G

mkdir -p ./out_log/mr
mkdir -p ./error_log/mr

MR_jid=$(/opt/slurm/bin/sbatch --parsable \
  -J MR_analysis \
  -c 2 \
  -p intel-sc3,amd-ep2,amd-ep2-short \
  -q normal \
  -a 1-106 \
  --ntasks-per-node 1 \
  --mem 8G \
  -o ./out_log/mr/${trait_name}_MR_comparison_%A_%a_out.txt \
  -e ./error_log/mr/${trait_name}_MR_comparison_%A_%a_error.txt \
  ${SCRIPT_DIR}/L2G/MR/MR_comparison.sh ${CONFIG})

echo "MR analsis: $MR_jid"

