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

mkdir -p ${OUTPUT}/V2G/summary
mkdir -p ${OUTPUT}/V2G/score
mkdir -p ${OUTPUT}/V2G/feature
mkdir -p ${OUTPUT}/V2G/plot

# ------------------------ CARMA results summary --------------------------------------
CARMA_bed_file="${OUTPUT}/CARMA/summary/${trait_name}_CARMA.bed"
CARMA_CS_file="${OUTPUT}/CARMA/summary/${trait_name}_CARMA.CS"

# Wu_adj_bed_file="${OUTPUT}/Wu_adj/summary/${trait_name}.bed"
# Wu_adj_CS_file="${OUTPUT}/Wu_adj/summary/${trait_name}.CS"


bed_file=${CARMA_bed_file}
CS_file=${CARMA_CS_file}
# bed_file=${Wu_adj_bed_file}
# CS_file=${Wu_adj_CS_file}

GWAS_LOCUS_file="${OUTPUT}/COJO/summary/${trait_name}.locus"
# ------------------------------------------------------------------------
#  V2G analysis
# ------------------------------------------------------------------------
Exon_bed=`yq .v2g.Exon_bed "${CONFIG}"`
ABC_bed=`yq .v2g.ABC_bed "${CONFIG}"`
EpiMap_bed=`yq .v2g.EpiMap_bed "${CONFIG}"`
RoadMap_bed=`yq .v2g.RoadMap_bed "${CONFIG}"`
PCHiC_bed=`yq .v2g.PCHiC_bed "${CONFIG}"`
bedtools=`yq .software.bedtools "${CONFIG}"`

R_functions_file=`yq .software.R_functions "${CONFIG}"`

gencode=`yq .gene.gencode "${CONFIG}"`
gencode_tss=`yq .gene.gencode_tss "${CONFIG}"`
gamma_gene_file=`yq .gene.gamma_gene "${CONFIG}"`

env=`yq .environment.R_421 "${CONFIG}"`
source activate $env

# ----
Rscript ${SCRIPT_DIR}/V2G/V2G_TSS.R \
	${trait_name} \
	${OUTPUT} \
	${gencode_tss} \
	${bed_file}

${bedtools} \
	intersect -wa -wb \
	-a ${bed_file} \
	-b ${Exon_bed} > ${OUTPUT}/V2G/summary/${trait_name}_Exon.annot.summary	

${bedtools} \
	intersect -wa -wb \
	-a ${bed_file} \
	-b ${ABC_bed} > ${OUTPUT}/V2G/summary/${trait_name}_ABC.annot.summary

${bedtools} \
	intersect -wa -wb \
	-a ${bed_file} \
	-b ${EpiMap_bed} > ${OUTPUT}/V2G/summary/${trait_name}_EpiMap.annot.summary

${bedtools} \
	intersect -wa -wb \
	-a ${bed_file} \
	-b ${RoadMap_bed} > ${OUTPUT}/V2G/summary/${trait_name}_RoadMap.annot.summary

${bedtools} \
	intersect -wa -wb \
	-a ${bed_file} \
	-b ${PCHiC_bed} > ${OUTPUT}/V2G/summary/${trait_name}_PCHiC.annot.summary

Rscript ${SCRIPT_DIR}/V2G/V2G_summary.R \
	${trait_name} \
	${OUTPUT} \
	${gencode} \
	${gamma_gene_file} \
	${R_functions_file} \
	${CS_file} \
	${GWAS_LOCUS_file}


