#!/bin/bash
set -e

# ------------------------------------------------------------------------
#  Input
# ------------------------------------------------------------------------
CONFIG=$1
omim_path=`yq .model_data.omim "${CONFIG}"`
clinvar_path=`yq .model_data.clinvar "${CONFIG}"`
mgi_path=`yq .model_data.mgi "${CONFIG}"`
gene_features=`yq .model_data.gene_features "${CONFIG}"`
trait_name=`yq .input.trait "${CONFIG}"`
SCRIPT_DIR=`yq .script.path "${CONFIG}"`
MESH_ID=D003924  # TODO
pharmap_path=`yq .gamma.pharmaprojects_0507_shufeng "${CONFIG}"`
model_path=`yq .model_model.model "${CONFIG}"`
feature_path=`yq .model_model.feature "${CONFIG}"`
OUTPUT=`yq .input.output "${CONFIG}"`

mkdir -p ${OUTPUT}/GAMMA/score

env=`yq .environment.python_ml "${CONFIG}"`
source activate $env

python ${SCRIPT_DIR}/GAMMA_ML/0.0_backend.py \
    --gamma ${OUTPUT}/GAMMA/feature/${trait_name}_GAMMA.feature \
    --mesh_id ${MESH_ID} \
    --output ${OUTPUT}/GAMMA/score/AI_score.csv \
    --omim_path $omim_path \
    --clinvar_path $clinvar_path \
    --mgi_path $mgi_path \
    --gene_features $gene_features \
    --model_path $model_path \
    --feature_path $feature_path \
    --pharmap_path $pharmap_path