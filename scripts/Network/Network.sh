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

mkdir -p ${OUTPUT}/Network/summary
mkdir -p ${OUTPUT}/Network/score
mkdir -p ${OUTPUT}/Network/feature
mkdir -p ${OUTPUT}/Network/plot


# ------------------------------------------------------------------------
#  Network analysis
# ------------------------------------------------------------------------
R_functions=`yq .software.R_functions "${CONFIG}"`
gamma_gene=`yq .gene.gamma_gene "${CONFIG}"`

env=`yq .environment.R_421 "${CONFIG}"`
source activate $env
Rscript  ${SCRIPT_DIR}/Network/Network_summary.R \
	${trait_name} \
	${OUTPUT} \
	${gamma_gene} \
	${R_functions}