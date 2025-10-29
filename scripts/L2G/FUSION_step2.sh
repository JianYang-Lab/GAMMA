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

mkdir -p ${OUTPUT}/FUSION/GWAS_munged
mkdir -p ${OUTPUT}/FUSION/detail
mkdir -p ${OUTPUT}/FUSION/summary

# ------------------------------------------------------------------------
#  FUSION analysis
# ------------------------------------------------------------------------
env=`yq .environment.R_421 "${CONFIG}"`
FUSION=`yq .software.fusion "${CONFIG}"`
LDSC_resources_dir=`yq .fusion.ldsc_resources_dir "${CONFIG}"`
LDREF=`yq .fusion.LDREF "${CONFIG}"`
QTL_list=`yq .fusion.QTL_list "${CONFIG}"`
source activate ${env}


# ----
qtl_i=${SLURM_ARRAY_TASK_ID}
GTEx_tissue=`head -n ${qtl_i} ${QTL_list} | tail -n1 | awk -F "\t" '{print $1}'`


for i in $(seq 1 22); do

SECONDS=0

Rscript ${FUSION} \
	--sumstats ${OUTPUT}/FUSION/GWAS_munged/${trait_name}.sumstats_new \
	--weights ${LDSC_resources_dir}/FUSION_WEIGHT/GTExv8.ALL.${GTEx_tissue}.pos \
	--weights_dir ${LDSC_resources_dir}/FUSION_WEIGHT \
	--ref_ld_chr ${LDREF}/1000G.EUR. \
	--chr ${i} \
	--out ${OUTPUT}/FUSION/detail/${trait_name}.GTExv8.ALL.${GTEx_tissue}.${i}.dat

elapsed_time=$SECONDS
result_file="${OUTPUT}/FUSION/detail/${trait_name}.GTExv8.ALL.${GTEx_tissue}.${i}.dat"
elapsed_time=$SECONDS

echo "FUSION analysis results have been written into [$result_file]"
echo "FUSION analysis completed: $current_time"
echo "FUSION analysis computational time: $(($elapsed_time / 3600)):$((($elapsed_time / 60) % 60)):$(($elapsed_time % 60))"

done

awk 'NR==1 || FNR>1' ${OUTPUT}/FUSION/detail/${trait_name}.GTExv8.ALL.${GTEx_tissue}.*.dat >  ${OUTPUT}/FUSION/summary/${trait_name}.GTExv8.ALL.${GTEx_tissue}.chrALL.dat
rm  ${OUTPUT}/FUSION/detail/${trait_name}.GTExv8.ALL.${GTEx_tissue}.*.dat

