# %%
import os
import pandas as pd
import numpy as np
import argparse

# %%
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gamma", help = "input gamma feature file")
parser.add_argument("-i", "--mesh_id", nargs='?', const='NO_MESH_INPUT', default='NO_MESH_INPUT', help = "input disease mesh id")
parser.add_argument("-o", "--output", help = "output path")
parser.add_argument("-omim", "--omim_path", help = "omim data path")
parser.add_argument("-clinvar", "--clinvar_path", help = "clinvar data path")
parser.add_argument("-mgi", "--mgi_path", help = "mgi data path")
parser.add_argument("-gene", "--gene_features", help = "gene features path")
parser.add_argument("-uniprot", "--uniprot_path", help = "uniprot protein family path")
parser.add_argument("-pharmap", "--pharmap_path", help = "pharmap current clinical results")

args = parser.parse_args()

gamma_path = args.gamma
mesh_id = args.mesh_id
out_path = args.output
omim_path = args.omim_path
clinvar_path = args.clinvar_path
mgi_path = args.mgi_path
gene_feature_path = args.gene_features
uniprot_path = args.uniprot_path
pharmap_path = args.pharmap_path

# %%
# get dir of out_path
out_dir = os.path.dirname(out_path)
os.makedirs(os.path.join(out_dir, 'gene_relations'), exist_ok=True)

# %%
df_gamma = pd.read_csv(gamma_path, sep='\t')
df_output = pd.read_csv(out_path)

omim = pd.read_csv(omim_path)
clinvar = pd.read_csv(clinvar_path)
mgi = pd.read_csv(mgi_path)

gene_features = pd.read_csv(gene_feature_path)
dummy_family = pd.read_csv(uniprot_path)
pharmap = pd.read_csv(pharmap_path)

clinvar = clinvar.rename(columns={'clinvar_MeSH_id': 'MeSH_id', 'clinvar_Ensembl':'gene_id'})

# %%
top_genes = df_output.sort_values(by='yhat_approve', ascending=False).reset_index(drop=True).iloc[:100]

df_info = df_gamma[['gene_id', 'gene_name']+list([i for i in df_gamma.columns if i.startswith('GAMMA')])]
df_info = df_info.merge(df_output[['entrez_id_single', 'gene_id', 'gene_name']+list([i for i in df_output.columns if i.startswith('yhat_')])], 
                        on=['gene_id', 'gene_name'])

if mesh_id != 'NO_MESH_INPUT':
    pharmap = pharmap.loc[pharmap['MeSH_id'] == mesh_id]
    df_info = df_info.merge(pharmap[['entrez_id_single', 'Highest Status Reached Value']], 
                            how='left', on=['entrez_id_single'])
    df_info.loc[:, 'Highest Status Reached Value'] = df_info['Highest Status Reached Value'].fillna(-1)
    
    if len(pharmap) == 0:
        mesh_id = 'MESH_ID_NOT_IN_PHARMAP'

def get_pathway_highest(gene_set):
    gene_set_info = df_info.loc[df_info['gene_id'].isin(gene_set)]
    GAMMA_highest = gene_set_info['GAMMA'].max()
    if GAMMA_highest >0:
        GAMMA_highest_genes = gene_set_info['gene_name'].loc[gene_set_info['GAMMA'] == GAMMA_highest].tolist()
        GAMMA_highest_genes = ', '.join(GAMMA_highest_genes)
    else:
        GAMMA_highest_genes = 'No gene with GAMMA > 0'
    GAMMA_ML_highest = gene_set_info['yhat_approve'].max()
    GAMMA_ML_highest_genes = gene_set_info['gene_name'].loc[gene_set_info['yhat_approve'] == GAMMA_ML_highest].tolist()
    GAMMA_ML_highest_genes = ', '.join(GAMMA_ML_highest_genes)
    if mesh_id == 'NO_MESH_INPUT':
        return [GAMMA_highest, GAMMA_highest_genes, GAMMA_ML_highest, GAMMA_ML_highest_genes, 'No MeSH input', 'No MeSH input']
    elif mesh_id == 'MESH_ID_NOT_IN_PHARMAP':
        return [GAMMA_highest, GAMMA_highest_genes, GAMMA_ML_highest, GAMMA_ML_highest_genes, 'No record found for MeSH', 'No record found for MeSH']
    else:
        pharmap_highest = gene_set_info['Highest Status Reached Value'].max()
        if pharmap_highest == -1:
            pharmap_highest_genes = 'No gene entered clinical status'
        else:
            pharmap_highest_genes = gene_set_info['gene_name'].loc[gene_set_info['Highest Status Reached Value'] == pharmap_highest].tolist()
            pharmap_highest_genes = ', '.join(pharmap_highest_genes)
        return [GAMMA_highest, GAMMA_highest_genes, GAMMA_ML_highest, GAMMA_ML_highest_genes, pharmap_highest, pharmap_highest_genes]

def get_gene_info(prediction):
    entrez_id = prediction['entrez_id_single']
    ensembl_id = prediction['gene_id']
    gene_symbol = prediction['gene_name']
    
    gene_relation_dir = os.path.join(out_dir, 'gene_relations', f"{ensembl_id}_{gene_symbol}")
    os.makedirs(gene_relation_dir, exist_ok=True)
    
    if mesh_id == 'NO_MESH_INPUT':
        disease_database_info = pd.DataFrame([prediction[['entrez_id_single', 'gene_id', 'gene_name']].values], columns=['entrez_id_single', 'gene_id', 'gene_name'])
        disease_database_info.to_csv(os.path.join(gene_relation_dir, f"disease_database_info.csv"), index=False)
    elif mesh_id != 'NO_MESH_INPUT':
        omim_info = omim.loc[(omim['MeSH_id'] == mesh_id) & (omim['gene_id'] == ensembl_id)]
        clinvar_info = clinvar.loc[(clinvar['MeSH_id'] == mesh_id) & (clinvar['gene_id'] == ensembl_id)]
        mgi_info = mgi.loc[(mgi['entrez_id_single'] == entrez_id) & (mgi['MeSH_id'] == mesh_id)]
        
        disease_database_info = pd.DataFrame([prediction[['entrez_id_single', 'gene_id', 'gene_name']].values], columns=['entrez_id_single', 'gene_id', 'gene_name'])
        disease_database_info = disease_database_info.merge(omim_info, how='left', on='gene_id')
        disease_database_info = disease_database_info.merge(clinvar_info, how='left', on=['gene_id', 'MeSH_id'])
        disease_database_info = disease_database_info.merge(mgi_info, how='left', on=['entrez_id_single', 'MeSH_id'])
        
        disease_database_info.to_csv(os.path.join(gene_relation_dir, f"disease_database_info.csv"), index=False)
    
    gene_feature_info = gene_features.loc[gene_features['gene_id'] == ensembl_id]
    
    basic_info = gene_feature_info[['gene_id', 'gene_name', 'Intracellular', 'Membrane', 'Secreted']]
    
    related_pathways1 = gene_feature_info[list(gene_feature_info.columns[59:127]) + list(gene_feature_info.columns[260:3026])]
    related_pathways1 = related_pathways1.columns[(related_pathways1 > 0).any()].tolist()

    related_genes1 = [list(gene_features['gene_id'].loc[gene_features[pathway_i]>0]) for pathway_i in related_pathways1]

    related_pathways2 = dummy_family.loc[dummy_family['gene_id'] == ensembl_id]
    related_pathways2 = related_pathways2.drop(columns=['uniprot_family_nan', 'gene_id'])
    related_pathways2 = related_pathways2.columns[(related_pathways2 > 0).any()].tolist()
    
    related_genes2 = [list(dummy_family['gene_id'].loc[dummy_family[family_i]>0]) for family_i in related_pathways2]
    
    # collecet related_genes3 from PPI.csv
    
    related_pathways = related_pathways1 + related_pathways2
    related_genes = related_genes1 + related_genes2
    
    related_outputs = [get_pathway_highest(gene_set) for gene_set in related_genes]

    output_df = pd.DataFrame(related_outputs, columns=['GAMMA_highest', 'GAMMA_highest_genes', 'GAMMA_ML_highest', 'GAMMA_ML_highest_genes', 'pharmap_highest', 'pharmap_highest_genes'])
    output_df.loc[:, 'pathway'] = related_pathways
    output_df = output_df.sort_values(by='GAMMA_highest', ascending=False).reset_index(drop=True)

    output_df = output_df[['pathway', 'GAMMA_highest', 'GAMMA_highest_genes', 'GAMMA_ML_highest', 'GAMMA_ML_highest_genes', 'pharmap_highest', 'pharmap_highest_genes']]
    output_df = output_df.rename(columns={
        'pathway': 'Pathway',
        'GAMMA_highest': 'Pathway Highest GAMMA Score',
        'GAMMA_highest_genes': 'Pathway Genes with Highest GAMMA Score',
        'GAMMA_ML_highest': 'Pathway Highest GAMMA-ML Score',
        'GAMMA_ML_highest_genes': 'Pathway Genes with Highest GAMMA-ML Score',
        'pharmap_highest': 'Pathway Highest Clinical Status',
        'pharmap_highest_genes': 'Pathway Genes with Highest Clinical Status'
    })
    
    output_df.to_csv(os.path.join(gene_relation_dir, f"related_pathways_info.csv"), index=False)
    basic_info.to_csv(os.path.join(gene_relation_dir, f"basic_info.csv"), index=False)
    
    return basic_info

# %%
top_genes.apply(get_gene_info, axis=1)

# %%
