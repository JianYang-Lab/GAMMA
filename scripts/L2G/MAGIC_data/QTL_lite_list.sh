#!/bin/bash


CONFIG="/storage/yangjianLab/guoyazhou/github/T2D_gene_prioritization/yaml_file/T2D_EUR.yaml"
OUTPUT="/storage/yangjianLab/guoyazhou/GAMMA_github/gamma-script/scripts/L2G/MAGIC_data/QTL_lite_besd"

QTL_list=`yq .magic.QTL_list "${CONFIG}"`
QTL_list_num=`cat ${QTL_list} | wc -l`

rm ${OUTPUT}/QTL_list_lite.txt

for ((i=1;i<=${QTL_list_num};i++))
do

qtl_i=${i}
qtl_name=`awk -F "\t" -v row=$qtl_i 'NR==row {print $1}' $QTL_list`
# qtl_data=`awk -F "\t" -v row=$qtl_i 'NR==row {print $2}' $QTL_list`
qtl_chr=`awk -F "\t" -v row=$qtl_i 'NR==row {print $3}' $QTL_list`
qtl_data="${OUTPUT}/besd/${qtl_name}"

echo -e ${qtl_name}"\t"${qtl_data}"\t"${qtl_chr} >> ${OUTPUT}/QTL_list_lite.txt

done