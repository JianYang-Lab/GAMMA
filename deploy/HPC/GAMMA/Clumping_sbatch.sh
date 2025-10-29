#!/usr/bin/env bash

CONFIG=$1

script_dir="/storage/yangjianLab/guoyazhou/GAMMA_github/gamma-script/scripts"
trait_name=`yq .input.trait "${CONFIG}"`

cd /storage/yangjianLab/guoyazhou/GAMMA_git/GWAS

mkdir -p ./out_log/clumping
mkdir -p ./error_log/clumping


jid=$(/opt/slurm/bin/sbatch --parsable \
	-J Clumping_analysis \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q normal \
	-a 1-22 \
	--ntasks-per-node 1 \
	--mem 12G \
	-o ./out_log/clumping/${trait_name}_clumping_step1_%A_%a_out.txt \
	-e ./error_log/clumping/${trait_name}_clumping_step1_%A_%a_error.txt \
	${script_dir}/GWAS/Clumping_step1.sh ${CONFIG}) \
	&& \
	/opt/slurm/bin/sbatch -d afterok:${jid} \
	-J Clumping_results \
	-c 2 \
	-p intel-sc3,amd-ep2,amd-ep2-short \
	-q normal \
	-a 1 \
	--ntasks-per-node 1 \
	--mem 12G \
	-o ./out_log/clumping/${trait_name}_clumping_step2_%A_%a_out.txt \
	-e ./error_log/clumping/${trait_name}_clumping_step2_%A_%a_error.txt \
	${script_dir}/GWAS/Clumping_step2_results.sh ${CONFIG}

