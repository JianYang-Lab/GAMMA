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

mkdir -p ${OUTPUT}/GSMR/outcome_data
mkdir -p ${OUTPUT}/GSMR/exposure_data
mkdir -p ${OUTPUT}/GSMR/new_HEIDI
mkdir -p ${OUTPUT}/GSMR/no_HEIDI
mkdir -p ${OUTPUT}/GSMR/summary


# ------------------------------------------------------------------------
#  GSMR analysis
#  step0: outcome file
# ------------------------------------------------------------------------
outcome_file=${OUTPUT}/GSMR/outcome_data/${trait_name}_gsmr.txt
awk '!seen[$1]++{print}' ${GWAS_DATA} > ${OUTPUT}/GSMR/outcome_data/${trait_name}_cojo_uniq.txt 
echo -e "${trait_name}\t${OUTPUT}/GSMR/outcome_data/${trait_name}_cojo_uniq.txt" > ${outcome_file}

