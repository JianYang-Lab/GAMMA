# %%
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import joblib
import pandas as pd
import numpy as np
import argparse

# %%
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gamma", help = "input gamma feature file")
parser.add_argument("-i", "--mesh_id", nargs='?', const='NO_MESH_INPUT', default='NO_MESH_INPUT', help = "input disease mesh id")
parser.add_argument("-o", "--output", help = "output path")
parser.add_argument("-m", "--model_path", help = "input model path")
parser.add_argument("-s", "--scaler_path", help = "input scaler path")
# parser.add_argument("-f", "--feature_path", help = "input feature path")
parser.add_argument("-omim", "--omim_path", help = "omim data path")
parser.add_argument("-clinvar", "--clinvar_path", help = "clinvar data path")
parser.add_argument("-mgi", "--mgi_path", help = "mgi data path")
parser.add_argument("-gene", "--gene_features", help = "gene features path")
parser.add_argument("-pharmap", "--pharmap_path", help = "pharmap current clinical results")

args = parser.parse_args()

gamma_path = args.gamma
mesh_id = args.mesh_id
out_path = args.output
model_path = args.model_path
scaler_path = args.scaler_path
# feature_path = args.feature_path
omim_path = args.omim_path
clinvar_path = args.clinvar_path
mgi_path = args.mgi_path
gene_feature_path = args.gene_features
pharmap_path = args.pharmap_path

# %%
# Device configuration
DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f'Predicting on: {DEVICE}')

# Batch size for prediction
batch_size = 2048

##########################
### MODEL DEFINITION (must match training)
##########################
class MLP(torch.nn.Module):
    def __init__(self, in_features, num_classes, num_hidden_1=128, num_hidden_2=32):
        super().__init__()
        self.my_network = torch.nn.Sequential(
            torch.nn.Linear(in_features, num_hidden_1, bias=False),
            torch.nn.LeakyReLU(),
            torch.nn.Dropout(0.2),
            torch.nn.BatchNorm1d(num_hidden_1),
            torch.nn.Linear(num_hidden_1, num_hidden_2, bias=False),
            torch.nn.LeakyReLU(),
            torch.nn.Dropout(0.2),
            torch.nn.BatchNorm1d(num_hidden_2),
        )
        self.fc = torch.nn.Linear(num_hidden_2, num_classes)

    def forward(self, x):
        x = self.my_network(x)
        logits = self.fc(x)
        probas = torch.sigmoid(logits)
        return logits, probas


##########################
### DATASET CLASS
##########################
class PredictDataset(Dataset):
    def __init__(self, feature_array):
        self.features = feature_array.astype(np.float32)
    def __len__(self):
        return len(self.features)
    def __getitem__(self, index):
        return self.features[index]


##########################
### LOAD MODEL AND SCALER
##########################

print("Loading model and scaler...")

# Load checkpoint
checkpoint = torch.load(model_path, map_location=DEVICE)
feature_names = checkpoint['feature_names']
num_features = len(feature_names)
num_classes = 5

print(f"Model expects {num_features} features")
print(f"Hyperparameters: {checkpoint['hyperparameters']}")

# Initialize model
model = MLP(in_features=num_features, num_classes=num_classes)
model.load_state_dict(checkpoint['model_state_dict'])
model.to(DEVICE)
model.eval()

# Load scaler
scaler = joblib.load(scaler_path)

print("Model and scaler loaded successfully!")

# %%
df = pd.read_csv(gamma_path, sep='\t')
df = df.rename(columns={'Gene_ID':'entrez_id_single'})
df.loc[:, 'GWAS'] = (~df['Lead_SNP'].isna()).astype(int)

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

clinvar = clinvar.rename(columns={'clinvar_Ensembl':'gene_id', 
                                  'clinvar_MeSH_id':'MeSH_id'})

omim = omim.loc[omim['MeSH_id']==mesh_id]
clinvar = clinvar.loc[clinvar['MeSH_id']==mesh_id]
mgi = mgi.loc[mgi['MeSH_id']==mesh_id]

df = pd.merge(df, omim.drop('MeSH_id', axis=1), on=['gene_id'], how='left')
df = pd.merge(df, clinvar.drop('MeSH_id', axis=1), on=['gene_id'], how='left')
df = pd.merge(df, mgi.drop('MeSH_id', axis=1), on=['entrez_id_single'], how='left')

df['OMIM'] = df['OMIM'].fillna(0)
df[clinvar.columns[2:]] = df[clinvar.columns[2:]].fillna(0)
df['MGI_human'] = df['MGI_human'].fillna(0)
df['MGI_mouse'] = df['MGI_mouse'].fillna(0)

# %%
df2 = pd.read_csv(gene_feature_path)

# %%
df = df.merge(df2.drop('gene_name', axis=1))

# %%
##########################
### VERIFY FEATURES
##########################
# Check if all required features are present
missing_features = set(feature_names) - set(df.columns)
if missing_features:
    raise ValueError(f"Missing features in new data: {missing_features}")

# Extract features in the correct order
df_features = df[feature_names]

print(f"Features extracted: {df_features.shape}")

# %%
##########################
### PREPROCESS DATA
##########################
print("Preprocessing data...")

df_features['PCHiC'] = df_features['PCHiC'].replace([np.inf, -np.inf], 500)
df_features['PoPS'] = df_features['PoPS'].fillna(0)

columns_w_na=[]
for i in df_features.columns:
    if df_features[i].isna().sum() > 0:
        if i.startswith('p_') or i.startswith('P_'):
            df_features[i] = df_features[i].fillna(1)
        elif i.startswith('z_'):
            df_features[i] = df_features[i].fillna(0)
        elif i.startswith('pp4'):
            df_features[i] = df_features[i].fillna(0)
        elif i.endswith('MAGIC'):
            df_features[i] = df_features[i].fillna(1)
        elif i.endswith('SMR'):
            df_features[i] = df_features[i].fillna(1)
        elif i.endswith('COLOC'):
            df_features[i] = df_features[i].fillna(0)
        elif i.endswith('FUSION'):
            df_features[i] = df_features[i].fillna(1)
        elif i=='MAGMA':
            df_features[i] = df_features[i].fillna(1)
        elif i=='mBATcombo':
            df_features[i] = df_features[i].fillna(1)
        elif i=='DistanceTSS':
            df_features[i] = df_features[i].fillna(df[i].max())
        else:
            print(i)
            # df_features[i] = df_features[i].fillna(df_features[i].mean())

features_scaled = scaler.transform(df_features)
print(f"Data preprocessed: {features_scaled.shape}")

# %%
##########################
### MAKE PREDICTIONS
##########################

print("Making predictions...")

# Create dataset and dataloader
predict_dataset = PredictDataset(features_scaled)
predict_loader = DataLoader(predict_dataset, batch_size=batch_size, shuffle=False)

# Collect predictions
all_probas = []
all_predicted_labels = []

with torch.no_grad():
    for batch_idx, features_batch in enumerate(predict_loader):
        features_batch = features_batch.to(DEVICE)
        logits, probas = model(features_batch)
        predicted_labels = (probas > 0.5).float()
        
        all_probas.append(probas.cpu())
        all_predicted_labels.append(predicted_labels.cpu())
        
        if batch_idx % 50 == 0:
            print(f"Processed batch {batch_idx}/{len(predict_loader)}")

# Concatenate all predictions
all_probas = torch.cat(all_probas, dim=0)
all_predicted_labels = torch.cat(all_predicted_labels, dim=0)

print("Predictions completed!")

# %%
##########################
### CREATE OUTPUT DATAFRAME
##########################

print("Creating output dataframe...")

# Create probability dataframe
df_probs = pd.DataFrame(
    all_probas.numpy(), 
    columns=["yhat_pre_cli", "yhat_phase_1", "yhat_phase_2", "yhat_phase_3", "yhat_approve"]
)

# Create binary predictions dataframe
df_predictions = pd.DataFrame(
    all_predicted_labels.numpy(), 
    columns=["pred_pre_cli", "pred_phase_1", "pred_phase_2", "pred_phase_3", "pred_approve"]
)

# Combine with original data
df_output = pd.concat([df[['entrez_id_single', 'gene_id', 'gene_name']], df_probs, df_predictions], axis=1)

# %%
##########################
### SAVE PREDICTIONS
##########################
phase_list = ['pre_cli', 'phase_1', 'phase_2', 'phase_3', 'approve']

print(f"Saving predictions to: {out_path}")
df_output.to_csv(out_path, index=False)

print("Done!")
print(f"\nSummary:")
print(f"  - Total samples: {len(df_output)}")
print(f"  - Predictions saved to: {out_path}")
print(f"\nPrediction statistics:")
for phase in phase_list:
    prob_col = f"yhat_{phase}"
    pred_col = f"pred_{phase}"
    mean_prob = df_output[prob_col].mean()
    positive_pct = df_output[pred_col].sum() / len(df_output) * 100
    print(f"  {phase}: Mean probability = {mean_prob:.4f}, Predicted positive = {positive_pct:.2f}%")

# %%
# Clear memory
del model
torch.cuda.empty_cache()