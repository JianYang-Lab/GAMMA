#!/usr/bin/env bash

CONFIG=$1

trait_name=`yq .input.trait "${CONFIG}"`
SCRIPT_DIR=`yq .script.path "${CONFIG}"`
WORK_DIR=`yq .script.work_path "${CONFIG}"`

#### ------------------------------------------------------------------------- 
#  GWAS analysis
##### ------------------------------------------------------------------------

mkdir -p ${WORK_DIR}/GWAS
cd ${WORK_DIR}/GWAS

#  Clumping analysis --------
mkdir -p ./out_log/clumping
mkdir -p ./error_log/clumping

clumping_jid1=$(/opt/slurm/bin/sbatch --parsable \
	-J Clumping_analysis \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q normal \
	-a 1-22 \
	--ntasks-per-node 1 \
	--mem 12G \
	-o ./out_log/clumping/${trait_name}_clumping_step1_%A_%a_out.txt \
	-e ./error_log/clumping/${trait_name}_clumping_step1_%A_%a_error.txt \
	${SCRIPT_DIR}/GWAS/Clumping_step1.sh ${CONFIG}) 

clumping_jid2=$(/opt/slurm/bin/sbatch --parsable \
	-d afterok:${clumping_jid1} \
	-J Clumping_results \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q normal \
	-a 1 \
	--ntasks-per-node 1 \
	--mem 6G \
	-o ./out_log/clumping/Clumping_step2_${trait_name}_%A_%a_out.txt \
	-e ./error_log/clumping/Clumping_step2_${trait_name}_%A_%a_error.txt \
	${SCRIPT_DIR}/GWAS/Clumping_step2_results.sh ${CONFIG})

echo "Clumping analsis: $clumping_jid1"
echo "Clumping analsis: $clumping_jid2"

#  COJO analysis --------
mkdir -p ./out_log/cojo
mkdir -p ./error_log/cojo

COJO_jid1=$(/opt/slurm/bin/sbatch --parsable \
	-J COJO_analysis \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q normal \
	-a 1-22 \
	--ntasks-per-node 1 \
	--mem 20G \
	-o ./out_log/cojo/${trait_name}_cojo_step1_%A_%a_out.txt \
	-e ./error_log/cojo/${trait_name}_cojo_step1_%A_%a_error.txt \
	${SCRIPT_DIR}/GWAS/COJO_step1.sh ${CONFIG})


COJO_jid2=$(/opt/slurm/bin/sbatch --parsable \
	-d afterok:${COJO_jid1} \
	-J COJO_results \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q normal \
	-a 1 \
	--ntasks-per-node 1 \
	--mem 6G \
	-o ./out_log/cojo/COJO_step2_${trait_name}_%A_%a_out.txt \
	-e ./error_log/cojo/COJO_step2_${trait_name}_%A_%a_error.txt \
	${SCRIPT_DIR}/GWAS/COJO_step2_results.sh ${CONFIG})

echo "COJO analsis: $COJO_jid1"
echo "COJO analsis: $COJO_jid2"

# MAGMA analysis ---------

mkdir -p ./out_log/MAGMA
mkdir -p ./error_log/MAGMA

MAGMA_jid1=$(/opt/slurm/bin/sbatch --parsable \
	-J MAGMA_analysis \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q huge \
	-a 1-22 \
	--ntasks-per-node 1 \
	--mem 30G \
	-o ./out_log/MAGMA/${trait_name}_MAGMA_step1_%A_%a_out.txt \
	-e ./error_log/MAGMA/${trait_name}_MAGMA_step1_%A_%a_error.txt \
	${SCRIPT_DIR}/GWAS/MAGMA_step1.sh ${CONFIG})
	
MAGMA_jid2=$(/opt/slurm/bin/sbatch --parsable \
	-d afterok:${MAGMA_jid1} \
	-J MAGMA_results \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q normal \
	-a 1 \
	--ntasks-per-node 1 \
	--mem 6G \
	-o ./out_log/MAGMA/MAGMA_step2_${trait_name}_%A_%a_out.txt \
	-e ./error_log/MAGMA/MAGMA_step2_${trait_name}_%A_%a_error.txt \
	${SCRIPT_DIR}/GWAS/MAGMA_step2_results.sh ${CONFIG})

echo "MAGMA analsis: $MAGMA_jid1"
echo "MAGMA analsis: $MAGMA_jid2"

# mBATcombo analysis ---------

mkdir -p ./out_log/mBATcombo
mkdir -p ./error_log/mBATcombo

mBATcombo_jid1=$(/opt/slurm/bin/sbatch --parsable \
	-J mBAT_analysis \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q huge \
	-a 1-22 \
	--ntasks-per-node 1 \
	--mem 20G \
	-o ./out_log/mBATcombo/${trait_name}_mBATcombo_step1_%A_%a_out.txt \
	-e ./error_log/mBATcombo/${trait_name}_mBATcombo_step1_%A_%a_error.txt \
	${SCRIPT_DIR}/GWAS/mBATcombo_step1.sh ${CONFIG})

mBATcombo_jid2=$(/opt/slurm/bin/sbatch --parsable \
	-d afterok:${mBATcombo_jid1} \
	-J mBAT_results \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q normal \
	-a 1 \
	--ntasks-per-node 1 \
	--mem 6G \
	-o ./out_log/mBATcombo/mBATcombo_step2_${trait_name}_%A_%a_out.txt \
	-e ./error_log/mBATcombo/mBATcombo_step2_${trait_name}_%A_%a_error.txt \
	${SCRIPT_DIR}/GWAS/mBATcombo_step2_results.sh ${CONFIG})

echo "mBATcombo analsis: $mBATcombo_jid1"
echo "mBATcombo analsis: $mBATcombo_jid2"

#### ------------------------------------------------------------------------
#  V2G analysis
#### ------------------------------------------------------------------------

mkdir -p ${WORK_DIR}/V2G
cd ${WORK_DIR}/V2G

# CARMA analysis ---------

mkdir -p ./out_log/CARMA
mkdir -p ./error_log/CARMA

# Depends on using clumping or COJO. Here we take Clumping.
# ---- clumping_jid2
# here you need to notice that clumping is already finished.
#### **** So it best to finish clumping first.

CARMA_jid1=$(/opt/slurm/bin/sbatch --parsable \
	-d afterok:${clumping_jid2} \
	-J CARMA_analysis \
	-c 10 \
	-p intel-sc3,amd-ep2 \
	-q huge \
	-a 1-100 \
	--ntasks-per-node 1 \
	--mem 80G \
	-o ./out_log/CARMA/${trait_name}_CARMA_%A_%a_out.txt \
	-e ./error_log/CARMA/${trait_name}_CARMA_%A_%a_error.txt \
	${SCRIPT_DIR}/V2G/CARMA_step1.sh ${CONFIG})

CARMA_jid2=$(/opt/slurm/bin/sbatch --parsable \
	-d afterok:${CARMA_jid1} \
	-J CARMA_results \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q normal \
	-a 1 \
	--ntasks-per-node 1 \
	--mem 10G \
	-o ./out_log/CARMA/CARMA_results_${trait_name}_%A_%a_out.txt \
	-e ./error_log/CARMA/CARMA_results_${trait_name}_%A_%a_error.txt \
	${SCRIPT_DIR}/V2G/CARMA_step2_results.sh ${CONFIG})

echo "CARMA analsis: $CARMA_jid1"
echo "CARMA analsis: $CARMA_jid2"

# Wu adj analysis ----------


mkdir -p ./out_log/Wu_adj
mkdir -p ./error_log/Wu_adj


#  ---------
#  Wu_adj analysis
#  ---------

Wu_adj_jid=$(/opt/slurm/bin/sbatch --parsable \
	-d afterok:${clumping_jid2} \
	-J Wu_adj_analysis \
	-c 2 \
	-p intel-sc3,amd-ep2 \
	-q normal \
	-a 1 \
	--ntasks-per-node 1 \
	--mem 30G \
	-o ./out_log/Wu_adj/${trait_name}_Wu_adj_%A_%a_out.txt \
	-e ./error_log/Wu_adj/${trait_name}_Wu_adj_%A_%a_error.txt \
	${SCRIPT_DIR}/V2G/Wu_adj.sh ${CONFIG})

echo "Wu adj analsis: $Wu_adj_jid"

# V2G analysis ---------

mkdir -p ./out_log/V2G
mkdir -p ./error_log/V2G

V2G_jid=$(/opt/slurm/bin/sbatch  --parsable \
	-d afterok:${CARMA_jid2} \
	-J V2G_analysis \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q normal \
	-a 1 \
	--ntasks-per-node 1 \
	--mem 30G \
	-o ./out_log/V2G/${trait_name}_v2g_%A_%a_out.txt \
	-e ./error_log/V2G/${trait_name}_v2g_%A_%a_error.txt \
	${SCRIPT_DIR}/V2G/V2G.sh ${CONFIG})

echo "V2G analsis: $V2G_jid"


# ------------------------------------------------------------------------
#  L2G analysis
# ------------------------------------------------------------------------
mkdir -p ${WORK_DIR}/L2G
cd ${WORK_DIR}/L2G


# SMR analysis ---------
mkdir -p ./out_log/smr
mkdir -p ./error_log/smr

QTL_list=`yq .magic.QTL_list "${CONFIG}"`
QTL_num=`cat ${QTL_list} | wc -l`
# here QTL_num == 418

SMR_jid=$(/opt/slurm/bin/sbatch --parsable \
  -J SMR_analysis \
  -c 5 \
  -p intel-sc3,amd-ep2,amd-ep2-short \
  -q huge \
  -a 1-${QTL_num} \
  --ntasks-per-node 1 \
  --mem 36G \
  -o ./out_log/smr/${trait_name}_SMR_%A_%a_out.txt \
  -e ./error_log/smr/${trait_name}_SMR_%A_%a_error.txt \
  ${SCRIPT_DIR}/L2G/SMR.sh ${CONFIG})

echo "SMR analsis: $SMR_jid"

# ------------------------------------------------------------------------
# MAGIC analysis
# ------------------------------------------------------------------------
mkdir -p ./out_log/magic
mkdir -p ./error_log/magic

MAGIC_jid=$(/opt/slurm/bin/sbatch --parsable \
  -d afterok:${clumping_jid2}:${SMR_jid} \
  -J MAGIC_analysis \
  -c 4 \
  -p intel-sc3,amd-ep2,amd-ep2-short \
  -q normal \
  -a 1 \
  --ntasks-per-node 1 \
  --mem 20G \
  -o ./out_log/magic/MAGIC_${trait_name}_%A_%a_out.txt \
  -e ./error_log/magic/MAGIC_${trait_name}_%A_%a_error.txt \
  ${SCRIPT_DIR}/MAGIC.sh ${CONFIG})

echo "MAGIC analysis: $MAGIC_jid"


# COLOC analysis ---------
mkdir -p ./out_log/coloc
mkdir -p ./error_log/coloc

# now selec 5e-5 probes, see COLOC_step0.sh
COLOC_jid=$(/opt/slurm/bin/sbatch --parsable \
	-J COLOC_analysis \
	-c 4 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q huge \
	-a 1-136 \
	--ntasks-per-node 1 \
	--mem 20G \
	-o ./out_log/coloc/${trait_name}_coloc_%A_%a_out.txt \
	-e ./error_log/coloc/${trait_name}_coloc_%A_%a_error.txt \
	${SCRIPT_DIR}/L2G/COLOC_step1.sh ${CONFIG})

echo "COLOC analsis: $COLOC_jid"


# FUSION analysis ---------
mkdir -p ./out_log/fusion
mkdir -p ./error_log/fusion

FUSION_jid1=$(/opt/slurm/bin/sbatch --parsable \
	-J FUSION_analysis \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q normal \
	-a 1 \
	--ntasks-per-node 1 \
	--mem 16G \
	-o ./out_log/fusion/FUSION_step1_${trait_name}_%A_%a_out.txt \
	-e ./error_log/fusion/FUSION_step1_${trait_name}_%A_%a_error.txt \
	${SCRIPT_DIR}/L2G/FUSION_step1.sh ${CONFIG}) 

FUSION_jid2=$(/opt/slurm/bin/sbatch --parsable \
	-d afterok:${FUSION_jid1} \
	-J FUSION_analysis \
	-c 4 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q huge \
	-a 1-49 \
	--ntasks-per-node 1 \
	--mem 16G \
	-o ./out_log/fusion/${trait_name}_FUSION_step2_%A_%a_out.txt \
	-e ./error_log/fusion/${trait_name}_FUSION_step2_%A_%a_error.txt \
	${SCRIPT_DIR}/L2G/FUSION_step2.sh ${CONFIG})
	
echo "FUSION analsis: $FUSION_jid2"
echo "FUSION analsis: $FUSION_jid1"

# ------------------------------------------------------------------------
# GSMR bash scrit ---------- 

mkdir -p ./out_log/gsmr
mkdir -p ./error_log/gsmr

# # Step 0: outcome format
# GSMR_jid0=$(/opt/slurm/bin/sbatch --parsable \
#     -J GSMR_step0 \
#     -c 2 \
#     -p intel-sc3,amd-ep2,amd-ep2-short \
#     -q normal \
#     -a 1 \
#     --ntasks-per-node 1 \
#     --mem 6G \
#     -o ./out_log/gsmr/GSMR_step0_${trait_name}_%A_%a_out.txt \
#     -e ./error_log/gsmr/GSMR_step0_${trait_name}_%A_%a_error.txt \
#     ${SCRIPT_DIR}/L2G/GSMR_step0.sh ${CONFIG})

# # Step 1: GSMR analysis
# GSMR_jid1=$(/opt/slurm/bin/sbatch --parsable \
#     -d afterok:$jid_0 \
#     -J GSMR_analysis \
#     -c 2 \
#     -p intel-sc3,amd-ep2,amd-ep2-short \
#     -q huge \
#     -a 1-11043 \
#     --ntasks-per-node 1 \
#     --mem 16G \
#     -o ./out_log/gsmr/${trait_name}_gsmr_%A_%a_out.txt \
#     -e ./error_log/gsmr/${trait_name}_gsmr_%A_%a_error.txt \
#     ${SCRIPT_DIR}/L2G/GSMR_step1.sh ${CONFIG})

# # Step 2: GSMR results
# GSMR_jid2=$(/opt/slurm/bin/sbatch --parsable \
# 	-d afterok:$jid_1 \
#     -J GSMR_results \
#     -c 2 \
#     -p intel-sc3,amd-ep2,amd-ep2-short \
#     -q normal \
#     -a 1 \
#     --ntasks-per-node 1 \
#     --mem 12G \
#     -o ./out_log/gsmr/GSMR_step2_${trait_name}_%A_%a_out.txt \
#     -e ./error_log/gsmr/GSMR_step2_${trait_name}_%A_%a_error.txt \
#     ${SCRIPT_DIR}/L2G/GSMR_step2_results.sh ${CONFIG})

# echo "GSMR analsis: $GSMR_jid0"
# echo "GSMR analsis: $GSMR_jid1"
# echo "GSMR analsis: $GSMR_jid2"

# ------------------------------------------------------------------------
# GSMR R script
mkdir -p ./out_log/mr_gsmr
mkdir -p ./error_log/mr_gsmr

GSMR_R_jid=$(/opt/slurm/bin/sbatch --parsable \
  -J MR_GSMR_analysis \
  -c 2 \
  -p intel-sc3,amd-ep2,amd-ep2-short \
  -q normal \
  -a 1-6 \
  --ntasks-per-node 1 \
  --mem 40G \
  -o ./out_log/mr_gsmr/${trait_name}_MR_comparison_%A_%a_out.txt \
  -e ./error_log/mr_gsmr/${trait_name}_MR_comparison_%A_%a_error.txt \
  ${SCRIPT_DIR}/L2G/MR/GSMR_MR.sh ${CONFIG})



# L2G analysis ---------- 
mkdir -p ./out_log/L2G
mkdir -p ./error_log/L2G

L2G_jid=$(/opt/slurm/bin/sbatch --parsable \
  -d afterok:${SMR_jid}:${COLOC_jid}:${FUSION_jid2}:${GSMR_R_jid} \
  -J L2G_analysis \
  -c 2 \
  -p intel-sc3,amd-ep2,amd-ep2-short \
  -q normal \
  -a 1 \
  --ntasks-per-node 1 \
  --mem 30G \
  -o ./out_log/L2G/${trait_name}_l2g_%A_%a_out.txt \
  -e ./error_log/L2G/${trait_name}_l2g_%A_%a_error.txt \
  ${SCRIPT_DIR}/L2G/L2G.sh ${CONFIG})

echo "L2G analsis: $L2G_jid"

# ------------------------------------------------------------------------
#  Network analysis
# ------------------------------------------------------------------------
mkdir -p ${WORK_DIR}/Network
cd ${WORK_DIR}/Network

# PoPS analysis --------------

mkdir -p ./out_log/pops
mkdir -p ./error_log/pops

# ---- step1: PoPS analysis

PoPS_jid=$(/opt/slurm/bin/sbatch --parsable \
	-d afterok:${MAGMA_jid2} \
	-J PoPS_analysis \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q normal \
	-a 1 \
	--ntasks-per-node 1 \
	--mem 30G \
	-o ./out_log/pops/${trait_name}_PoPS_%A_%a_out.txt \
	-e ./error_log/pops/${trait_name}_PoPS_%A_%a_error.txt \
	${SCRIPT_DIR}/Network/PoPS_step1.sh ${CONFIG})

echo "PoPS analsis: $PoPS_jid"

# DEPICT analysis --------------
mkdir -p ./out_log/depict
mkdir -p ./error_log/depict

DEPCIT_jid=$(/opt/slurm/bin/sbatch --parsable \
	-J DEPICT_analysis \
	-c 10 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q normal \
	-a 1 \
	--ntasks-per-node 1 \
	--mem 240G \
	-o ./out_log/depict/${trait_name}_DEPICT_%A_%a_out.txt \
	-e ./error_log/depict/${trait_name}_DEPICT_%A_%a_error.txt \
	${SCRIPT_DIR}/Network/DEPICT.sh ${CONFIG})

echo "DEPICT analsis: $DEPCIT_jid"

# RWR and PPR analysis --------------
mkdir -p ./out_log/rwr_ppr
mkdir -p ./error_log/rwr_ppr

RWR_PPR_jid=$(/opt/slurm/bin/sbatch --parsable \
	-d afterok:${V2G_jid}:${L2G_jid} \
	-J RWR_PPR_analysis \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q normal \
	-a 1 \
	--ntasks-per-node 1 \
	--mem 30G \
	-o ./out_log/rwr_ppr/RWR_PPR_step1_${trait_name}_%A_%a_out.txt \
	-e ./error_log/rwr_ppr/RWR_PPR_step1_${trait_name}_%A_%a_error.txt \
	${SCRIPT_DIR}/Network/RWR_PPR.sh ${CONFIG})

echo "RWR and PPR analsis: $RWR_PPR_jid"

# Network analysis --------------
mkdir -p ./out_log/network
mkdir -p ./error_log/network

Network_jid=$(/opt/slurm/bin/sbatch --parsable \
  -J Network_analysis \
  -c 2 \
  -p intel-sc3,amd-ep2,amd-ep2-short \
  -q normal \
  -a 1 \
  --ntasks-per-node 1 \
  --mem 12G \
  -o ./out_log/network/${trait_name}_network_%A_%a_out.txt \
  -e ./error_log/network/${trait_name}_network_%A_%a_error.txt \
  ${SCRIPT_DIR}/Network/Network.sh ${CONFIG})

echo "Network analsis: $Network_jid"

# ------------------------------------------------------------------------
#  GAMMA analysis
# ------------------------------------------------------------------------
mkdir -p ${WORK_DIR}/GAMMA 
cd ${WORK_DIR}/GAMMA

mkdir -p ./out_log/GAMMA
mkdir -p ./error_log/GAMMA

GAMMA_jid=$(/opt/slurm/bin/sbatch --parsable \
  -J GAMMA_analysis \
  -c 2 \
  -p intel-sc3,amd-ep2,amd-ep2-short \
  -q normal \
  -a 1 \
  --ntasks-per-node 1 \
  --mem 32G \
  -o ./out_log/GAMMA/${trait_name}_gamma_%A_%a_out.txt \
  -e ./error_log/GAMMA/${trait_name}_gamma_%A_%a_error.txt \
  ${SCRIPT_DIR}/GAMMA/GAMMA.sh ${CONFIG})

echo "GAMMA analsis: $GAMMA_jid"



GAMMA_ML_jid=$(/opt/slurm/bin/sbatch --parsable \
  -d afterok:${GAMMA_jid} \
  -J GAMMA_ML_analysis \
  -c 2 \
  -p intel-sc3,amd-ep2,amd-ep2-short \
  -q normal \
  -a 1 \
  --ntasks-per-node 1 \
  --mem 10G \
  -o ./out_log/GAMMA/${trait_name}_gamma_ml_%A_%a_out.txt \
  -e ./error_log/GAMMA/${trait_name}_gamma_ml_%A_%a_error.txt \
  ${SCRIPT_DIR}/GAMMA_ML/GAMMA_ML.sh ${CONFIG})

echo "GAMMA ML analsis: $GAMMA_ML_jid"


