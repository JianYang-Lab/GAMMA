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


mkdir -p ${OUTPUT}/CARMA/LD
mkdir -p ${OUTPUT}/CARMA/Result
mkdir -p ${OUTPUT}/CARMA/summary
mkdir -p ${OUTPUT}/CARMA/tmp

# ------------------------------------------------------------------------
#  CARMA analysis
# ------------------------------------------------------------------------
REFERENCE=`yq .reference.reference_bfile "${CONFIG}"`
plink1_9=`yq .software.plink1_9 "${CONFIG}"`
plink2=`yq .software.plink2 "${CONFIG}"`

env=`yq .environment.R_421 "${CONFIG}"`
source activate $env

# ----
# here Clumping or COJO:
# we take Clumping to do CARMA analysis
# 1-${locus_num_array}

# locus_file="${OUTPUT}/Clumping/summary/${trait_name}.locus"
locus_file="${OUTPUT}/COJO/summary/${trait_name}.locus"
locus_num_array=`awk 'NR>1{print}' ${locus_file} | wc -l`

for ((i=${SLURM_ARRAY_TASK_ID}; i<=${locus_num_array}; i=i+${locus_num_array})); do

	# ------------------------
	# Check if the file exists
	echo $i
	locus_line=$((i + 1))
	locus=$(awk -v locus_line="$locus_line" 'NR == locus_line {print $5}' ${locus_file})
	output_file=${OUTPUT}/CARMA/Result/${trait_name}_${locus}.txt.gz
	if [ -f $output_file ]; then
		echo "----------------- File already exists --------------: $output_file"
	else

		echo "----------------- Running CARMA --------------"
		Rscript ${SCRIPT_DIR}/V2G/CARMA_step1.R \
			${GWAS_DATA} \
			${trait_name} \
			${i} \
			${OUTPUT} \
			${locus_file} \
			${REFERENCE} \
			${plink1_9} \
			${plink2}

	fi

done
