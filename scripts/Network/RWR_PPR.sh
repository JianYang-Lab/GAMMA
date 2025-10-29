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

mkdir -p ${OUTPUT}/RWR_PPR/summary

# ------------------------------------------------------------------------
#  RWR and PPR analysis
# ------------------------------------------------------------------------
env=`yq .environment.python_rwr_ppr "${CONFIG}"`
source activate $env
ppi_file=`yq .rwr_ppr.ppi_file "${CONFIG}"`


python ${SCRIPT_DIR}/Network/RWR_PPR.py \
	${trait_name} \
	${OUTPUT} \
	${ppi_file}