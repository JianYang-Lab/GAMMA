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


mkdir -p ${OUTPUT}/PoPS/MAGMA/annotation
mkdir -p ${OUTPUT}/PoPS/MAGMA/results
mkdir -p ${OUTPUT}/PoPS/PoPS_score

# ------------------------------------------------------------------------
#  PoPS analysis
# ------------------------------------------------------------------------
# step1 MAGMA analysis
# Do not need to run, just choose from MAGMA results.

MAGMA=`yq .software.magma "${CONFIG}"`
REFERENCE=`yq .reference.reference_bfile_GRCh37 "${CONFIG}"`
gene_annot=`yq .magma.gene_annot "${CONFIG}"`

# ----
i=${SLURM_ARRAY_TASK_ID}

${MAGMA} \
	--bfile ${REFERENCE}_chr${i} \
	--gene-annot ${gene_annot} \
	--pval ${GWAS_DATA} ncol=N \
	--gene-model snp-wise=mean \
	--out ${OUTPUT}/MAGMA/detail/${trait_name}_chr${i}


