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


mkdir -p ${OUTPUT}/COLOC/detail
mkdir -p ${OUTPUT}/COLOC/summary
mkdir -p ${OUTPUT}/COLOC/tmp


# ------------------------------------------------------------------------
#  COLOC analysis
# ------------------------------------------------------------------------
reference_freq=`yq .reference.reference_freq "${CONFIG}"`
QTL_list=`yq .smr.QTL_list "${CONFIG}"`

qtl_i=$2
qtl_name=`head -n ${qtl_i} ${QTL_list} | tail -n1 | awk -F "\t" '{print $1}'`
qtl_data=`head -n ${qtl_i} ${QTL_list} | tail -n1 | awk -F "\t" '{print $2}'`
qtl_chr=`head -n ${qtl_i} ${QTL_list} | tail -n1 | awk -F "\t" '{print $3}'`
qtl_n=`head -n ${qtl_i} ${QTL_list} | tail -n1 | awk -F "\t" '{print $4}'`

COLOC=`yq .software.coloc "${CONFIG}"`
R_functions=`yq .software.R_functions "${CONFIG}"`
QTL_dir=`yq .coloc.QTL_dir "${CONFIG}"`

# ----
# for i in $(seq 1 22); do
i=${SLURM_ARRAY_TASK_ID}
SECONDS=0

if [ -f "${QTL_dir}/QTL_data/${qtl_name}_chr${i}.txt" ]; then
Rscript ${COLOC} \
	--gwas ${GWAS_DATA} \
	--reference_freq ${reference_freq} \
	--qtl_query  ${QTL_dir}/QTL_data/${qtl_name}_chr${i}.txt \
	--qtl_number ${qtl_n} \
	--R_functions ${R_functions} \
	--out ${OUTPUT}/COLOC/detail/${trait_name}_${qtl_name}_chr${i}.coloc
fi

elapsed_time=$SECONDS
result_file="${OUTPUT}/COLOC/detail/${trait_name}_${qtl_name}_chr${i}.coloc"
current_time=$(date "+%Y-%m-%d %H:%M:%S")

echo "COLOC analysis results have been written into [$result_file]"
echo "COLOC analysis completed: $current_time"
echo "COLOC analysis computational time: $(($elapsed_time / 3600)):$((($elapsed_time / 60) % 60)):$(($elapsed_time % 60))"

# done

# awk 'NR==1 || FNR>1' ${OUTPUT}/COLOC/detail/${trait_name}_${qtl_name}_chr*.coloc > ${OUTPUT}/COLOC/summary/${trait_name}_${qtl_name}_chrALL.coloc 

# # rm ${OUTPUT}/COLOC/detail/${trait_name}_${QTLName}_chr*.coloc