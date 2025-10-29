# %%
from sklearn.ensemble import RandomForestClassifier
import joblib
import pandas as pd
import numpy as np
import argparse

# %%
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gamma", help = "input gamma feature file")
parser.add_argument("-i", "--mesh_id", help = "input disease mesh id")
parser.add_argument("-o", "--output", help = "output path")
parser.add_argument("-m", "--model_path", help = "input model path")
parser.add_argument("-f", "--feature_path", help = "input feature path")
parser.add_argument("-omim", "--omim_path", help = "omim data path")
parser.add_argument("-clinvar", "--clinvar_path", help = "clinvar data path")
parser.add_argument("-mgi", "--mgi_path", help = "mgi data path")
parser.add_argument("-gene", "--gene_features", help = "gene features path")
parser.add_argument("-pharmap", "--pharmap_path", help = "pharmap current clinical results")

args = parser.parse_args()

# %%
gamma_path = args.gamma
mesh_id = args.mesh_id
model_path = args.model_path
feature_path = args.feature_path
omim_path = args.omim_path
clinvar_path = args.clinvar_path
mgi_path = args.mgi_path
gene_feature_path = args.gene_features
out_path = args.output
pharmap_path = args.pharmap_path

# %%
# gamma_path = '/storage/yangjianLab/guoyazhou/GAMMA_git_data/GAMMA/feature/T2D_GAMMA.feature'
# mesh_id = 'D003924'
# model_path = '/storage/yangjianLab/sunshufeng/gamma_v4/0.0.1_valid_feature_model.joblib'
# feature_path = '/storage/yangjianLab/sunshufeng/gamma_v4/0.0.2_feature_list.txt'
# omim_path = '/storage/yangjianLab/sunshufeng/gamma_v4/data/OMIM_processed.csv'
# clinvar_path = '/storage/yangjianLab/sunshufeng/gamma_v4/data/clinvar.csv'
# mgi_path = '/storage/yangjianLab/sunshufeng/gamma_v4/data/MGI_processed.csv'
# gene_feature_path = '/storage/yangjianLab/sunshufeng/gamma_v4/0.0.3_gene_features.csv'
# out_path = '/storage/yangjianLab/sunshufeng/gamma_v4/0.0.4_tmp_test_out.csv'
# pharmap_path = '/storage/yangjianLab/sunshufeng/gamma_v4/data/pharmap_processed.csv'

# %%
clf = joblib.load(model_path)

# %%
df = pd.read_csv(gamma_path, sep='\t')

# %%
# Collect SMR
smr_test_df = df[[i for i in df.columns if i.startswith('p_SMR')]]
smr_test = ((smr_test_df * (~smr_test_df.isna()).sum()) < 5e-2).values
heidi_test = (df[[i for i in df.columns if (i.startswith('p_HEIDI')) and (not i.endswith('OPERA'))]] > 5e-2).values
result_array = np.logical_and(smr_test, heidi_test)
df.loc[:, 'GAMMA_SMR'] = result_array.sum(axis=1)

# Collect COLOC
df.loc[:, 'GAMMA_COLOC'] = (df[[i for i in df.columns if i.startswith('pp4')]] > 0.7).sum(axis=1)

# Collect FUSION
fusion_test_df = df[[i for i in df.columns if i.startswith('p_FUSION')]]
df.loc[:, 'GAMMA_FUSION'] = ((fusion_test_df * ((~fusion_test_df.isna()).sum())) < 5e-2).sum(axis=1)

# %%
omim = pd.read_csv(omim_path)
clinvar = pd.read_csv(clinvar_path)
mgi = pd.read_csv(mgi_path)

# %%
clinvar = clinvar.rename(columns={'clinvar_Ensembl':'gene_id', 
                                  'clinvar_MeSH_id':'MeSH_id'})
df = df.rename(columns={'Gene_ID':'entrez_id_single'})

# %%
omim = omim.loc[omim['MeSH_id']==mesh_id]
clinvar = clinvar.loc[clinvar['MeSH_id']==mesh_id]
mgi = mgi.loc[mgi['MeSH_id']==mesh_id]

df = pd.merge(df, omim.drop('MeSH_id', axis=1), on=['gene_id'], how='left')
df = pd.merge(df, clinvar.drop('MeSH_id', axis=1), on=['gene_id'], how='left')
df = pd.merge(df, mgi.drop('MeSH_id', axis=1), on=['entrez_id_single'], how='left')

df['OMIM'] = df['OMIM'].fillna(0)
df[clinvar.columns[2:]] = df[clinvar.columns[2:]].fillna(0)
df['MGI_score'] = df['MGI_score'].fillna(0)

# %%
df2 = pd.read_csv(gene_feature_path)

# %%
df = df.merge(df2.drop('gene_name', axis=1))

# %%
df['PCHiC'] = df['PCHiC'].replace([np.inf, -np.inf], 500)
df['PoPS'] = df['PoPS'].fillna(0)

columns_w_na=[]
for i in df.columns[18:]:
    if df[i].isna().sum() > 0:
        if i.startswith('p_') or i.startswith('P_'):
            df[i] = df[i].fillna(1)
        elif i.startswith('z_'):
            df[i] = df[i].fillna(0)
        elif i.startswith('pp4'):
            df[i] = df[i].fillna(0)
        elif i.endswith('MAGIC'):
            df[i] = df[i].fillna(1)
        elif i.endswith('SMR'):
            df[i] = df[i].fillna(1)
        elif i.endswith('COLOC'):
            df[i] = df[i].fillna(0)
        elif i.endswith('FUSION'):
            df[i] = df[i].fillna(1)
        elif i=='MAGMA':
            df[i] = df[i].fillna(1)
        elif i=='mBATcombo':
            df[i] = df[i].fillna(1)
        elif i=='DistanceTSS':
            df[i] = df[i].fillna(df[i].max())
        else:
            print(i)
            # df[i] = df[i].fillna(df[i].mean())

# %%
with open(feature_path) as file:
    features = [line.rstrip() for line in file]
X_test = np.array(df[features])

y_pred = clf.predict_proba(X_test)

# %%
out_df = df[df.columns[:10]]

phase_list = ['pre_cli', 'phase_1', 'phase_2', 'phase_3', 'approve']
for i in range(len(phase_list)):
    out_df['yhat_'+phase_list[i]] = y_pred[i][:,1]

# %%
pharmap = pd.read_csv(pharmap_path)
pharmap = pharmap.loc[pharmap['MeSH_id']==mesh_id]
out_df = pd.merge(out_df, pharmap[['entrez_id_single', 'Highest Status Reached Value']], how='left')
out_df['Highest Status Reached Value'] = out_df['Highest Status Reached Value'].fillna(-1)

# %%
out_df.to_csv(out_path, index=False)

# %%
