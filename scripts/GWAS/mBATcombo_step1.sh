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


mkdir -p ${OUTPUT}/mBATcombo/detail
mkdir -p ${OUTPUT}/mBATcombo/summary

# ------------------------------------------------------------------------
#  mBATcombo analysis
# ------------------------------------------------------------------------
GCTA=`yq .software.gcta "${CONFIG}"`
REFERENCE=`yq .reference.reference_bfile "${CONFIG}"`
gene_list=`yq .mBAT.gene_list "${CONFIG}"`

# ----
i=${SLURM_ARRAY_TASK_ID}
${GCTA} --bfile ${REFERENCE}_chr${i} \
	--mBAT-combo ${GWAS_DATA} \
	--mBAT-gene-list ${gene_list} \
	--mBAT-print-all-p \
	--chr ${i} \
	--out ${OUTPUT}/mBATcombo/detail/${trait_name}_mBATcombo_chr${i}

