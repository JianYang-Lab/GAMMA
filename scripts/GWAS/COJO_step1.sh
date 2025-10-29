#!/bin/bash
set -e

# ------------------------------------------------------------------------
#  Input
# ------------------------------------------------------------------------
CONFIG=$1
SCRIPT_DIR=`yq .script.path "${CONFIG}"`
GWAS_DATA=`yq .input.gwas "${CONFIG}"`
trait_name=`yq .input.trait "${CONFIG}"`
OUTPUT=`yq .input.output "${CONFIG}"`

mkdir -p ${OUTPUT}/COJO/detail
mkdir -p ${OUTPUT}/COJO/summary

# ------------------------------------------------------------------------
#  COJO analysis
# ------------------------------------------------------------------------
GCTA=`yq .software.gcta "${CONFIG}"`
REFERENCE=`yq .reference.reference_bfile "${CONFIG}"`

# ----
i=${SLURM_ARRAY_TASK_ID}
${GCTA} \
  --cojo-file ${GWAS_DATA} \
  --bfile ${REFERENCE}_chr${i} \
  --chr ${i} \
  --maf 0.01 \
  --cojo-slct \
  --cojo-p 5e-8 \
  --diff-freq 1 \
  --cojo-wind 10000 \
  --thread-num 10 \
  --out ${OUTPUT}/COJO/detail/${trait_name}_chr${i}


