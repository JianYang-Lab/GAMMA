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
clump_field=`yq .clumping.field "${CONFIG}"`
if [ "$clump_field" = "null" ]; then
  clump_field="P"
fi

mkdir -p ${OUTPUT}/Clumping/detail
mkdir -p ${OUTPUT}/Clumping/summary

# ------------------------------------------------------------------------
#  Clumping analysis
# ------------------------------------------------------------------------
plink1_9=`yq .software.plink1_9 "${CONFIG}"`
REFERENCE=`yq .reference.reference_bfile "${CONFIG}"`

# ----
i=${SLURM_ARRAY_TASK_ID}
${plink1_9} \
    --bfile ${REFERENCE} \
    --chr ${i} \
    --maf 0.01 \
    --clump ${GWAS_DATA} \
    --clump-p1 5e-8 \
    --clump-p2 5e-8 \
    --clump-r2 0.05 \
    --clump-kb 1000 \
    --clump-snp-field SNP \
    --clump-field ${clump_field} \
    --out ${OUTPUT}/Clumping/detail/${trait_name}_chr${i}


