#!/bin/bash
# CONFIG_template="/storage/yangjianLab/guoyazhou/GAMMA_github/gamma-script/deploy/HPC/GAMMA/template.yaml"
CONFIG_template=$1
trait_name=$2
GWAS_DATA=$3
WORK_DIR=${4:-"/storage/yangjianLab/guoyazhou/GAMMA_git"}
OUTPUT=${5:-"/storage/yangjianLab/guoyazhou/GAMMA_git_data"}
SCRIPT_DIR=${6:-"/storage/yangjianLab/guoyazhou/GAMMA_github/gamma-script/scripts"}
report_dir=$7
MESH_ID=$8

mkdir -p ${OUTPUT}/GWAS/COJO_format
COJO_file=${OUTPUT}/GWAS/COJO_format/${trait_name}.txt

# CONFIG_template=`yq .yaml.template "${CONFIG_template}"`
mkdir -p ${WORK_DIR}/yaml_file 
CONFIG=${WORK_DIR}/yaml_file/${trait_name}.yaml 
cp ${CONFIG_template} ${CONFIG}



yq -i ".input.trait = \"$trait_name\"" "$CONFIG"
yq -i ".input.gwas_raw = \"$GWAS_DATA\"" "$CONFIG"
yq -i ".input.gwas = \"$COJO_file\"" "$CONFIG"
yq -i ".input.output = \"$OUTPUT\"" "$CONFIG"
yq -i ".input.mesh_id = \"$MESH_ID\"" "$CONFIG"

yq -i ".script.work_path = \"$WORK_DIR\"" "$CONFIG"
yq -i ".script.path = \"$SCRIPT_DIR\"" "$CONFIG"

yq -i ".yaml.config = \"$CONFIG\"" "$CONFIG"
yq -i ".input.report_dir = \"$report_dir\"" "$CONFIG"
