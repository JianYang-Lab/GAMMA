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
#  COJO results
# ------------------------------------------------------------------------
gencode=`yq .gene.gencode "${CONFIG}"`
env=`yq .environment.R_421 "${CONFIG}"`
source activate $env

# ----
awk 'NR==1 || FNR>1'  ${OUTPUT}/COJO/detail/${trait_name}_chr*.jma.cojo >  ${OUTPUT}/COJO/summary/${trait_name}.jma.cojo

if [ -f "${OUTPUT}/COJO/summary/${trait_name}.jma.cojo" ]; then

	Rscript ${SCRIPT_DIR}/GWAS/COJO_step2_results.R ${trait_name} ${OUTPUT}

	Rscript ${SCRIPT_DIR}/GWAS/COJO_step2_results_locus_window.R ${trait_name} ${OUTPUT} 1000000 ${gencode}
	Rscript ${SCRIPT_DIR}/GWAS/COJO_step2_results_locus_window.R ${trait_name} ${OUTPUT} 750000 ${gencode}
	Rscript ${SCRIPT_DIR}/GWAS/COJO_step2_results_locus_window.R ${trait_name} ${OUTPUT} 500000 ${gencode}
	Rscript ${SCRIPT_DIR}/GWAS/COJO_step2_results_locus_window.R ${trait_name} ${OUTPUT} 250000 ${gencode}

fi

