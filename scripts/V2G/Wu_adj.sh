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


mkdir -p ${OUTPUT}/Wu_adj/LD
mkdir -p ${OUTPUT}/Wu_adj/Result
mkdir -p ${OUTPUT}/Wu_adj/summary


# ------------------------------------------------------------------------
# Wu_adj analysis
# ------------------------------------------------------------------------
reference_bfile=`yq .reference.reference_bfile "${CONFIG}"`
reference_all_bim=`yq .reference.reference_all_bim "${CONFIG}"`
reference_freq=`yq .reference.reference_freq "${CONFIG}"`

plink1_9=`yq .software.plink1_9 "${CONFIG}"`
plink2=`yq .software.plink2 "${CONFIG}"`

env=`yq .environment.R_421 "${CONFIG}"`
source activate $env

# --------------------------
# No matter Clumping or COJO, we only take the CHR/POS/SNP information.
# if input is COJO:

# COJO_or_Clumping_file="${OUTPUT}/Clumping/summary/${trait_name}.clumped"
# COJO_or_Clumping_file="${OUTPUT}/COJO/summary/${trait_name}.jma.cojo"
locus_path="${OUTPUT}/Clumping/summary/${trait_name}.locus"

Rscript ${SCRIPT_DIR}/V2G/Wu_adj.R \
	${GWAS_DATA} \
	${trait_name} \
	${OUTPUT} \
	${locus_path} \
	${plink1_9} \
	${reference_all_bim} \
	${reference_freq} \
	${reference_bfile}
	
