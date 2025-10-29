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
#  mBATcombo results
# ------------------------------------------------------------------------
awk 'NR==1 || FNR>1' ${OUTPUT}/mBATcombo/detail/${trait_name}_mBATcombo_chr*.gene.assoc.mbat > ${OUTPUT}/mBATcombo/summary/${trait_name}.gene.assoc.mbat
