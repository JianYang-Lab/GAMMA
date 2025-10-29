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

mkdir -p ${OUTPUT}/GAMMA/score
mkdir -p ${OUTPUT}/GAMMA/feature
mkdir -p ${OUTPUT}/GAMMA/plot
mkdir -p ${OUTPUT}/GWAS/manhattan_plot


# ------------------------------------------------------------------------
#  Network analysis
# ------------------------------------------------------------------------
# now MeSH_id = ""
R_functions_file=`yq .software.R_functions "${CONFIG}"`
gamma_gene_file=`yq .gene.gamma_gene "${CONFIG}"`
Pharmaprojects_data_file=`yq .gamma.pharmaprojects "${CONFIG}"`
reference_all_bim_file=`yq .reference.reference_all_bim "${CONFIG}"`

# MeSH id need to be updated
# GWAS trait MeSH id list.
gwas_mesh=`yq .mesh.gwas_mesh "${CONFIG}"` 
MeSH_id=`cat ${gwas_mesh} | awk -F "\t" -v trait_name="${trait_name}" '{if ($1 == trait_name) print $3}'`
# MeSH_id="NA"

# T2D
# MeSH_id="D003924"  # TODO: get from CONFIG
echo ${trait_name}
echo ${OUTPUT}
echo ${Pharmaprojects_data_file}
echo ${gamma_gene_file}
echo ${R_functions_file}

env=`yq .environment.R_421 "${CONFIG}"`
source activate ${env}
echo ${MeSH_id}

echo "------------"
Rscript ${SCRIPT_DIR}/GAMMA/GAMMA_summary.R \
	${trait_name} \
	${OUTPUT} \
	${Pharmaprojects_data_file} \
	${gamma_gene_file} \
	${R_functions_file} \
	${MeSH_id}


# Rscript ${SCRIPT_DIR}/GAMMA/GAMMA_summary_text.R \
	# ${trait_name} \
	# ${OUTPUT}


Rscript ${SCRIPT_DIR}/GAMMA/manhattan_plot.R \
	${trait_name} \
	${GWAS_DATA} \
	${OUTPUT} \
	${reference_all_bim_file}





