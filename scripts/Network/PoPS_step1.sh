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
# ----
# MAGMA results
# head -n 1 ${OUTPUT}/PoPS/MAGMA/results/${trait_name}_chr1.genes.out > ${OUTPUT}/PoPS/MAGMA/results/${trait_name}.chrALL.genes.out 
# for chr in {1..22}
# do
#	tail -n +2 ${OUTPUT}/PoPS/MAGMA/results/${trait_name}_chr${chr}.genes.out >> ${OUTPUT}/PoPS/MAGMA/results/${trait_name}.chrALL.genes.out
# done

# head -n 2 ${OUTPUT}/PoPS/MAGMA/results/${trait_name}_chr1.genes.raw > ${OUTPUT}/PoPS/MAGMA/results/${trait_name}.chrALL.genes.raw 
# for chr in {1..22}
# do
#	tail -n +3 ${OUTPUT}/PoPS/MAGMA/results/${trait_name}_chr${chr}.genes.raw >> ${OUTPUT}/PoPS/MAGMA/results/${trait_name}.chrALL.genes.raw
# done


# step2 PoPS analysis
PoPS=`yq .software.pops "${CONFIG}"`
gene_annot=`yq .pops.gene_annot "${CONFIG}"`
pops_feature=`yq .pops.pops_feature "${CONFIG}"`
control_feature=`yq .pops.control_feature "${CONFIG}"`

env=`yq .environment.python_pops "${CONFIG}"`
source activate $env

python ${PoPS} \
	--verbose \
	--gene_annot_path ${gene_annot} \
	--feature_mat_prefix ${pops_feature} \
	--num_feature_chunks 116 \
	--magma_prefix ${OUTPUT}/MAGMA/summary/${trait_name}_chrALL \
	--control_features_path ${control_feature} \
	--out_prefix ${OUTPUT}/PoPS/PoPS_score/${trait_name}_PoPS
