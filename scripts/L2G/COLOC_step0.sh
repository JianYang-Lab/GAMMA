#!/bin/bash
set -e

# ------------------------------------------------------------------------
#  Input
# ------------------------------------------------------------------------
CONFIG=$1
SMR=`yq .software.smr "${CONFIG}"`
# QTL_list=`yq .smr.QTL_list "${CONFIG}"`
QTL_list=`yq .magic.QTL_list "${CONFIG}"`
QTL_dir=`yq .coloc.QTL_dir "${CONFIG}"`


qtl_i=${SLURM_ARRAY_TASK_ID}
# **** Here we take qtl as a unit.
# 1-78: eQTL
# 79-129: sQTL
# 130-131: mQTL
# 132-136: pQTL

qtl_name=`head -n ${qtl_i} ${QTL_list} | tail -n1 | awk -F "\t" '{print $1}'`
qtl_data=`head -n ${qtl_i} ${QTL_list} | tail -n1 | awk -F "\t" '{print $2}'`
qtl_chr=`head -n ${qtl_i} ${QTL_list} | tail -n1 | awk -F "\t" '{print $3}'`
qtl_n=`head -n ${qtl_i} ${QTL_list} | tail -n1 | awk -F "\t" '{print $4}'`



# ------------------------------------------------------------------------
#  QTL data
# ------------------------------------------------------------------------

if [ ! -f "${QTL_dir}/QTL_data/${qtl_name}_chr${i}.txt" ]; then

for i in $(seq 1 22); do

    if [ "$qtl_chr" = "TRUE" ]; then
        
		QTL_data="${qtl_data}${i}"
		
		${SMR} --beqtl-summary ${QTL_data} \
			--probe-chr ${i} \
			--descriptive-cis \
			--peqtl-cis 5e-5 \
			--out ${QTL_dir}/QTL_top/${qtl_name}_chr${i}

		awk -F "\t" 'NR>1 {print $1}' ${QTL_dir}/QTL_top/${qtl_name}_chr${i}.cis.summary.txt > ${QTL_dir}/QTL_top/${qtl_name}_chr${i}.cis.probe.list 

    else
        QTL_data="${qtl_data}"
		
		${SMR} --beqtl-summary ${QTL_data} \
			--descriptive-cis \
			--peqtl-cis 5e-5 \
			--out ${QTL_dir}/QTL_top/${qtl_name}

		awk -F "\t" -v chr=${i} 'NR > 1 && $2 == chr {print}' ${QTL_dir}/QTL_top/${qtl_name}.cis.summary.txt > ${QTL_dir}/QTL_top/${qtl_name}_chr${i}.cis.probe.list  
		
    fi


	${SMR} --beqtl-summary ${QTL_data} \
		--extract-probe ${QTL_dir}/QTL_top/${qtl_name}_chr${i}.cis.probe.list  \
		--query 1 \
		--out ${QTL_dir}/QTL_data/${qtl_name}_chr${i}

done

else
    echo "File already exists."
fi
