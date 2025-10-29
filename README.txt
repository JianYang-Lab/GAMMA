# GAMMA pipeline

## Input
CONFIG_template=$1
trait_name=$2
GWAS_DATA=$3
WORK_DIR=${4:-"/storage/yangjianLab/guoyazhou/GAMMA_git"}
OUTPUT=${5:-"/storage/yangjianLab/guoyazhou/GAMMA_git_data"}
SCRIPT_DIR=${6:-"/storage/yangjianLab/guoyazhou/GAMMA_github/gamma-script/scripts"}

## Config first
./scripts/config.sh \
    ./deploy/HPC/template.yaml \
    ${trait_name} \
    ${GWAS_DATA} \
    ${WORK_DIR} \
    ${OUTPUT} \
    ${SCRIPT_DIR}

## Check GWAS data 
./scripts/gwas_format.sh ${CONFIG}

## GAMMA pipeline
CONFIG=${WORK_DIR}/yaml_file/${trait_name}.yaml 
./deploy/HPC/GAMMA_sbatch.sh ${CONFIG}
