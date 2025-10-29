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
#  MAGMA results
# ------------------------------------------------------------------------
head -n 1 ${OUTPUT}/MAGMA/detail/${trait_name}_chr1.genes.out > ${OUTPUT}/MAGMA/summary/${trait_name}_chrALL.genes.out
for chr in {1..22}
do
	tail -n +2 ${OUTPUT}/MAGMA/detail/${trait_name}_chr${chr}.genes.out >> ${OUTPUT}/MAGMA/summary/${trait_name}_chrALL.genes.out 
done

head -n 2 ${OUTPUT}/MAGMA/detail/${trait_name}_chr1.genes.raw > ${OUTPUT}/MAGMA/summary/${trait_name}_chrALL.genes.raw
for chr in {1..22}
do
	tail -n +3 ${OUTPUT}/MAGMA/detail/${trait_name}_chr${chr}.genes.raw >> ${OUTPUT}/MAGMA/summary/${trait_name}_chrALL.genes.raw 
done

