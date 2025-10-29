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

mkdir -p ${OUTPUT}/MAGIC/plot
mkdir -p ${OUTPUT}/MAGIC/summary
mkdir -p ${OUTPUT}/MAGIC/gwas
mkdir -p ${OUTPUT}/MAGIC/results

# ------------------------------------------------------------------------
#  MAGIC analysis
# ------------------------------------------------------------------------
magic_functions_file=`yq .magic.R_functions "${CONFIG}"`
gencode_file=`yq .gene.gencode "${CONFIG}"`
CpG_link_file=`yq .magic.CpG_link "${CONFIG}"`
hQTL_link_file=`yq .magic.hQTL_link "${CONFIG}"`
caQTL_link_file=`yq .magic.caQTL_link "${CONFIG}"`
reference_bim_file=`yq .reference.reference_all_bim "${CONFIG}"`
QTL_name_list_file=`yq .magic.QTL_name_list "${CONFIG}"`


# awk 'NR==1 || FNR>1' ${OUTPUT}/SMR/summary/${trait_name}_*_chrALL.msmr > ${OUTPUT}/MAGIC/results/${trait_name}_ALL.msmr

env=`yq .environment.R_421 "${CONFIG}"`
source activate ${env}

Rscript ${SCRIPT_DIR}/L2G/MAGIC.R \
    ${trait_name} \
    ${OUTPUT} \
    ${magic_functions_file} \
    ${gencode_file} \
    ${CpG_link_file} \
    ${hQTL_link_file} \
    ${caQTL_link_file} \
    ${GWAS_DATA} \
    ${reference_bim_file} \
    ${QTL_name_list_file}


Rscript ${SCRIPT_DIR}/L2G/MAGIC_HEIDI.R \
    ${trait_name} \
    ${OUTPUT} \
    ${magic_functions_file} \
    ${gencode_file} \
    ${CpG_link_file} \
    ${hQTL_link_file} \
    ${caQTL_link_file} \
    ${GWAS_DATA} \
    ${reference_bim_file} \
    ${QTL_name_list_file}




