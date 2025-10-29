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
#  CARMA results
# ------------------------------------------------------------------------
env=`yq .environment.R_421 "${CONFIG}"`
source activate $env

Rscript ${SCRIPT_DIR}/V2G/CARMA_step2_results.R \
	${trait_name} \
	${OUTPUT}



