#!/bin/bash
#SBATCH --job-name=generate_QTL_lite_besd
#SBATCH -c 2
#SBATCH -o ./QTL_lite_besd/out_log/QTL_lite_%A_%a_out.txt
#SBATCH -e ./QTL_lite_besd/error_log/QTL_lite_%A_%a_error.txt
#SBATCH -p amd-ep2,intel-sc3
#SBATCH --array=451
####SBATCH --array=1-441
#SBATCH --ntasks-per-node=1
#SBATCH --mem=120G
#SBATCH --qos=normal

# ------------------------------------------------------------------------
#  Input
# ------------------------------------------------------------------------
CONFIG="/storage/yangjianLab/guoyazhou/github/T2D_gene_prioritization/yaml_file/T2D_EUR.yaml"

SMR=`yq .software.smr "${CONFIG}"`
# SMR="/storage/yangjianLab/guoyazhou/bin/smr"
SMR="/storage/yangjianLab/guoyazhou/software/SMR/smr_new/smr-1.4.0-linux-x86_64/smr"

QTL_list=`yq .magic.QTL_list "${CONFIG}"`
qtl_i=${SLURM_ARRAY_TASK_ID}
qtl_name=`awk -F "\t" -v row=$qtl_i 'NR==row {print $1}' $QTL_list`
qtl_data=`awk -F "\t" -v row=$qtl_i 'NR==row {print $2}' $QTL_list`
qtl_chr=`awk -F "\t" -v row=$qtl_i 'NR==row {print $3}' $QTL_list`


OUTPUT="/storage/yangjianLab/guoyazhou/GAMMA_github/gamma-script/scripts/L2G/MAGIC_data/QTL_lite_besd"
mkdir -p ${OUTPUT}/probe
mkdir -p ${OUTPUT}/query
mkdir -p ${OUTPUT}/besd


if [ "$qtl_chr" = "TRUE" ]; then
	for i in $(seq 1 22); do
		QTL_data="${qtl_data}${i}"
		${SMR} --beqtl-summary ${QTL_data} \
			--query 5e-8 \
			--out ${OUTPUT}/query/${qtl_name}_query_chr${i}

		cat ${OUTPUT}/query/${qtl_name}_query_chr${i}.txt | awk -F "\t" 'NR>1 {print $7}' | uniq > ${OUTPUT}/probe/${qtl_name}_probe_chr${i}.txt

		${SMR} --beqtl-summary ${QTL_data} \
			--extract-probe ${OUTPUT}/probe/${qtl_name}_probe_chr${i}.txt \
			--make-besd --out ${OUTPUT}/besd/${qtl_name}_chr${i}
		
	done

else
    QTL_data="${qtl_data}"

	${SMR} --beqtl-summary ${QTL_data} \
		--query 5e-8 \
		--out ${OUTPUT}/query/${qtl_name}_query 

	cat ${OUTPUT}/query/${qtl_name}_query.txt | awk -F "\t" 'NR>1 {print $7}' | uniq  > ${OUTPUT}/probe/${qtl_name}_probe.txt

	${SMR} --beqtl-summary ${QTL_data} \
		--extract-probe ${OUTPUT}/probe/${qtl_name}_probe.txt \
		--make-besd --out ${OUTPUT}/besd/${qtl_name}

fi
