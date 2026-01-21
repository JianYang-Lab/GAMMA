#!/bin/bash

# for i in `seq 1 245`
# do
# trait_name=`head -n $i ${GAMMA_HOME}/guoyazhou/GAMMA_git/GWAS_trait_list/GWAS_latest_list.txt | tail -n1 | cut -f1`
# GWAS_DATA=`head -n $i ${GAMMA_HOME}/guoyazhou/GAMMA_git/GWAS_trait_list/GWAS_latest_list.txt | tail -n1 | cut -f2`
# CONFIG="${GAMMA_HOME}/guoyazhou/GAMMA_git/yaml_file/${trait_name}.yaml"
# echo $i;echo $trait_name; echo $GWAS_DATA; echo $CONFIG
# ${GAMMA_HOME}/guoyazhou/GAMMA_github/gamma-script/deploy/HPC/hpc_update_config.sh ${CONFIG} ${trait_name} ${GWAS_DATA}
# done

# cat ${GAMMA_HOME}/guoyazhou/GAMMA_git_data/SMR/supplementary_file/SMR_Portal_GWAS_list_new.txt | while read line
# do
# trait_name=`echo $line | cut -f 1`
# GWAS_DATA=`echo $line | cut -f 2`
# echo $trait_name
# ${GAMMA_HOME}/guoyazhou/GAMMA_github/gamma-script/deploy/HPC/hpc_update_config.sh   ${GAMMA_HOME}/guoyazhou/GAMMA_git/yaml_file/${trait_name}.yaml ${trait_name} ${GWAS_DATA}
# done

# ls ${GAMMA_HOME}/guoyazhou/GAMMA_git/yaml_file/*yaml | while read line
# do
# CONFIG=${line}
# trait_name=`yq .input.trait "${CONFIG}"`
# GWAS_DATA=`yq .input.gwas_raw "${CONFIG}"`
# echo $CONFIG
# echo $trait_name
# echo $GWAS_DATA

# ${GAMMA_HOME}/guoyazhou/GAMMA_github/gamma-script/deploy/HPC/hpc_update_config.sh ${CONFIG} ${trait_name} ${GWAS_DATA}

# done


CONFIG=$1
trait_name=$2
GWAS_DATA=$3

# trait_name=`yq .input.trait "${CONFIG}"`
# GWAS_DATA=`yq .input.gwas_raw "${CONFIG}"`
# OUTPUT=`yq .input.output "${CONFIG}"`
# SCRIPT_DIR=`yq .script.path "${CONFIG}"`
# WORK_DIR=`yq .script.work_path "${CONFIG}"`

# CONFIG_template=`yq .yaml.template "${CONFIG_template}"`
CONFIG_template="${GAMMA_HOME}/guoyazhou/GAMMA_github/gamma-script/deploy/HPC/GAMMA/template.yaml"

cp ${CONFIG_template} ${CONFIG}

yq -i ".input.trait = \"$trait_name\"" "$CONFIG"
yq -i ".input.gwas_raw = \"$GWAS_DATA\"" "$CONFIG"
yq -i ".input.gwas = \"$GWAS_DATA\"" "$CONFIG"

clump_field=`awk 'NR==1 {print $7}' ${GWAS_DATA}`
yq -i ".clumping.field = \"$clump_field\"" "$CONFIG"

# yq -i ".input.output = \"$OUTPUT\"" "$CONFIG"

# yq -i ".script.work_path = \"$WORK_DIR\"" "$CONFIG"
# yq -i ".script.path = \"$SCRIPT_DIR\"" "$CONFIG"

yq -i ".yaml.config = \"$CONFIG\"" "$CONFIG"



