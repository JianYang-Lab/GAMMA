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


# ------------------------------------------------------------------------
#  GSMR results
# ------------------------------------------------------------------------
awk 'NR==1 || FNR>1' ${OUTPUT}/GSMR/new_HEIDI/${trait_name}_*.gsmr > ${OUTPUT}/GSMR/summary/${trait_name}.gsmr
rm ${OUTPUT}/GSMR/new_HEIDI/${trait_name}_*.gsmr

