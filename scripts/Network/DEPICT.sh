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

mkdir -p ${OUTPUT}/DEPICT/GWAS_input
mkdir -p ${OUTPUT}/DEPICT/cfg
mkdir -p ${OUTPUT}/DEPICT/output


# ------------------------------------------------------------------------
# DEPICT analysis
# ------------------------------------------------------------------------
DEPICT=`yq .software.depict "${CONFIG}"`
snp_loc_GRCh37=`yq .depict.snp_loc_GRCh37 "${CONFIG}"`
cfg_demo=`yq .depict.cfg_demo "${CONFIG}"`
env=`yq .environment.python_depict "${CONFIG}"`
source activate ${env}
plink1_9=`yq .software.plink1_9 "${CONFIG}"`


awk -v OFS='\t' 'NR == FNR {value[$1] = $2 OFS $3; next} {print $0, value[$1]}' ${snp_loc_GRCh37} ${GWAS_DATA} > ${OUTPUT}/DEPICT/GWAS_input/${trait_name}_depict.txt
awk -v OFS='\t' '{if($9 != "Y" && $9 != "X" && $9 != "" && $10 != "" && $7 != "NA" && $1 != "NA") print $1, $9,$10,$7}' ${OUTPUT}/DEPICT/GWAS_input/${trait_name}_depict.txt > ${OUTPUT}/DEPICT/GWAS_input/${trait_name}_depict.clean.txt
sed -i $'1 i \tSNP\tChr\tPos\tp' ${OUTPUT}/DEPICT/GWAS_input/${trait_name}_depict.clean.txt


cp ${cfg_demo} ${OUTPUT}/DEPICT/cfg/${trait_name}.cfg 
cfg_file=${OUTPUT}/DEPICT/cfg/${trait_name}.cfg
sed -i "7c analysis_path: ${OUTPUT}/DEPICT/output" ${cfg_file}
sed -i "40c plink_executable: ${plink1_9}" ${cfg_file}
sed -i "13c gwas_summary_statistics_file: ${OUTPUT}/DEPICT/GWAS_input/${trait_name}_depict.clean.txt" ${cfg_file}
sed -i "19c label_for_output_files: ${trait_name}_clean" ${cfg_file}


python ${DEPICT} ${cfg_file}

