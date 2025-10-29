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
# step1: Munge GWAS summary data
env=`yq .environment.python_ldsc "${CONFIG}"`
source activate ${env}

munge_sumstats=`yq .fusion.munge_sumstats "${CONFIG}"`
LDSC_resources_dir=`yq .fusion.ldsc_resources_dir "${CONFIG}"`

# ----
python ${munge_sumstats} \
	--sumstats ${GWAS_DATA} \
	--out ${OUTPUT}/FUSION/GWAS_munged/${trait_name} \
	--chunksize 500000 \
	--merge-alleles ${LDSC_resources_dir}/hapmap3_snps/w_hm3.snplist


if [ -f "${OUTPUT}/FUSION/GWAS_munged/${trait_name}.sumstats" ]; then
	rm ${OUTPUT}/FUSION/GWAS_munged/${trait_name}.sumstats
fi

gunzip ${OUTPUT}/FUSION/GWAS_munged/${trait_name}.sumstats.gz
awk '$4 != ""' ${OUTPUT}/FUSION/GWAS_munged/${trait_name}.sumstats > ${OUTPUT}/FUSION/GWAS_munged/${trait_name}.sumstats_new

