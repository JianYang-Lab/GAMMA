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


mkdir -p ${OUTPUT}/CARMA/LD
mkdir -p ${OUTPUT}/CARMA/Result
mkdir -p ${OUTPUT}/CARMA/summary
mkdir -p ${OUTPUT}/CARMA/tmp

# ------------------------------------------------------------------------
#  CARMA analysis
# ------------------------------------------------------------------------
REFERENCE=`yq .reference.reference_bfile "${CONFIG}"`
plink1_9=`yq .software.plink1_9 "${CONFIG}"`
plink2=`yq .software.plink2 "${CONFIG}"`

env=`yq .environment.R_421 "${CONFIG}"`
source activate $env

# ----
# here Clumping or COJO:
# we take Clumping to do CARMA analysis
# 1-100

# locus_num_array=`awk 'NR>1{print}' ${OUTPUT}/Clumping/summary/${trait_name}.locus | wc -l`

# for ((i=${SLURM_ARRAY_TASK_ID}; i<=${locus_num_array}; i=i+100)); do

i=${SLURM_ARRAY_TASK_ID}

	Rscript ${SCRIPT_DIR}/V2G/CARMA_step1.R \
		${GWAS_DATA} \
		${trait_name} \
		${i} \
		${OUTPUT} \
		${OUTPUT}/COJO/summary/${trait_name}.locus \
		${REFERENCE} \
		${plink1_9} \
		${plink2}

# done
