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

mkdir -p ${OUTPUT}/GSMR/outcome_data
mkdir -p ${OUTPUT}/GSMR/exposure_data
mkdir -p ${OUTPUT}/GSMR/new_HEIDI
mkdir -p ${OUTPUT}/GSMR/no_HEIDI
mkdir -p ${OUTPUT}/GSMR/summary


# ------------------------------------------------------------------------
#  GSMR analysis
# ------------------------------------------------------------------------
outcome_file=${OUTPUT}/GSMR/outcome_data/${trait_name}_gsmr.txt

GCTA=`yq .software.gcta "${CONFIG}"`
reffile=`yq .gsmr.reffile "${CONFIG}"`
expofile=`yq .gsmr.expofile "${CONFIG}"`
expo_dir=`yq .gsmr.expo_dir "${CONFIG}"`

# ----
i=${SLURM_ARRAY_TASK_ID}
prot=`sed -n ${i}'p' ${expofile} | awk '{print $1}'`
if [ -z "$prot" ]; then
    echo "no prot in line $i, skip"
    exit 0
fi
expo=${expo_dir}/${prot}.txt

# step 1. run GSMR analysis
${GCTA} \
  --mbfile ${reffile} \
  --maf 0.01 \
  --diff-freq 0.2 \
  --gsmr2-beta \
  --gsmr-file ${expo} ${outcome_file} \
  --gsmr-direction 0 \
  --threads 1 \
  --clump-r2 0.05 \
  --gsmr-snp-min 2 \
  --effect-plot \
  --out ${OUTPUT}/GSMR/new_HEIDI/${trait_name}_${prot}


