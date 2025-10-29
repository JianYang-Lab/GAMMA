#!/usr/bin/env bash

CONFIG=$1

trait_name=`yq .input.trait "${CONFIG}"`
GWAS_DATA=`yq .input.gwas "${CONFIG}"`
SCRIPT_DIR=`yq .script.path "${CONFIG}"`
WORK_DIR=`yq .script.work_path "${CONFIG}"`


cd ${WORK_DIR}
mkdir -p ./out_log/clumping
mkdir -p ./error_log/clumping
mkdir -p ./out_log/smr
mkdir -p ./error_log/smr
mkdir -p ./out_log/magic
mkdir -p ./error_log/magic


# ------------------------------------------------------------------------
# Clumping analysis
# ------------------------------------------------------------------------

Clumping_jid=$(/opt/slurm/bin/sbatch --parsable \
  -J Clumping_analysis \
  -c 2 \
  -p intel-sc3,amd-ep2,amd-ep2-short \
  -q normal \
  -a 1 \
  --ntasks-per-node 1 \
  --mem 16G \
  -o ./out_log/clumping/${trait_name}_Clumping_%A_%a_out.txt \
  -e ./error_log/clumping/${trait_name}_Clumping_%A_%a_error.txt \
  ${SCRIPT_DIR}/Clumping.sh ${CONFIG})

echo "Clumping analysis: $Clumping_jid"

# ------------------------------------------------------------------------
# SMR analysis
# ------------------------------------------------------------------------

QTL_list=`yq .magic.QTL_list "${CONFIG}"`
QTL_num=`cat ${QTL_list} | wc -l`
# here QTL_num == 418

SMR_jid=$(/opt/slurm/bin/sbatch --parsable \
  -J SMR_analysis \
  -c 5 \
  -p intel-sc3,amd-ep2,amd-ep2-short \
  -q normal \
  -a 1-${QTL_num} \
  --ntasks-per-node 1 \
  --mem 36G \
  -o ./out_log/smr/${trait_name}_SMR_%A_%a_out.txt \
  -e ./error_log/smr/${trait_name}_SMR_%A_%a_error.txt \
  ${SCRIPT_DIR}/SMR.sh ${CONFIG})

echo "SMR analysis: $SMR_jid"

# ------------------------------------------------------------------------
# MAGIC analysis
# ------------------------------------------------------------------------

MAGIC_jid=$(/opt/slurm/bin/sbatch --parsable \
  -d afterok:${Clumping_jid}:${SMR_jid} \
  -J MAGIC_analysis \
  -c 4 \
  -p intel-sc3,amd-ep2,amd-ep2-short \
  -q normal \
  -a 1 \
  --ntasks-per-node 1 \
  --mem 20G \
  -o ./out_log/magic/MAGIC_${trait_name}_%A_%a_out.txt \
  -e ./error_log/magic/MAGIC_${trait_name}_%A_%a_error.txt \
  ${SCRIPT_DIR}/MAGIC.sh ${CONFIG})

echo "MAGIC analysis: $MAGIC_jid"

