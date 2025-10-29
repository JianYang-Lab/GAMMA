#!/usr/bin/env bash

CONFIG=$1

trait_name=`yq .input.trait "${CONFIG}"`
SCRIPT_DIR=`yq .script.path "${CONFIG}"`
WORK_DIR=`yq .script.work_path "${CONFIG}"`


# -----------------------------------------------------------------------
#  L2G analysis
# ------------------------------------------------------------------------
mkdir -p ${WORK_DIR}/L2G
cd ${WORK_DIR}/L2G


# SMR analysis ---------
mkdir -p ./out_log/smr
mkdir -p ./error_log/smr




QTL_list=`yq .magic.QTL_list "${CONFIG}"`
QTL_num=`cat ${QTL_list} | wc -l`
# here QTL_num == 418

SMR_jid=$(/opt/slurm/bin/sbatch --parsable \
  -J SMR_analysis \
  -c 5 \
  -p intel-sc3,amd-ep2,amd-ep2-short \
  -q huge \
  -a 1-${QTL_num} \
  --ntasks-per-node 1 \
  --mem 36G \
  -o ./out_log/smr/${trait_name}_SMR_%A_%a_out.txt \
  -e ./error_log/smr/${trait_name}_SMR_%A_%a_error.txt \
  ${SCRIPT_DIR}/L2G/SMR.sh ${CONFIG})

echo "SMR analsis: $SMR_jid"

