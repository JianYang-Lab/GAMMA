# %%
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import GroupKFold
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.utils.class_weight import compute_class_weight

# %%
import os
import re
import sys
import sklearn.metrics
from sklearn.model_selection import GroupKFold, GridSearchCV, GroupShuffleSplit, RandomizedSearchCV

phase_list = ['pre_cli', 'phase_1', 'phase_2', 'phase_3', 'approve']
y_dict = {
    -1: [0,0,0,0,0],
    0 : [1,0,0,0,0],
    1 : [1,1,0,0,0],
    2 : [1,1,1,0,0],
    3 : [1,1,1,1,0],
    4 : [1,1,1,1,1]
}

# %%
# FULL feature data: 5e-07_2e+06_df_OC1_RP_MAGIC_KeReGo_Ufamily_GHexp.parquet
data_path = '/storage/yangjianLab/sunshufeng/gamma_v5/data'
data_matr = 'locus2e6_disease215_detailGAMMA_OC_networks_ppi_GOembedding_cancer_property_NAmean.parquet'

df_origin = pd.read_parquet(os.path.join(data_path, data_matr))
df_origin.loc[:, 'GWAS'] = (~df_origin['Lead_SNP'].isna()).astype(int)

# %%
df_origin['PCHiC'] = df_origin['PCHiC'].replace([np.inf, -np.inf], 500)
df_origin['PoPS'] = df_origin['PoPS'].fillna(0)

columns_w_na=[]
for i in df_origin.columns[18:]:
    if df_origin[i].isna().sum() > 0:
        if i.startswith('p_') or i.startswith('P_'):
            df_origin[i] = df_origin[i].fillna(1)
        elif i.startswith('z_'):
            df_origin[i] = df_origin[i].fillna(0)
        elif i.startswith('pp4'):
            df_origin[i] = df_origin[i].fillna(0)
        elif i.endswith('MAGIC'):
            df_origin[i] = df_origin[i].fillna(1)
        elif i.endswith('SMR'):
            df_origin[i] = df_origin[i].fillna(1)
        elif i.endswith('COLOC'):
            df_origin[i] = df_origin[i].fillna(0)
        elif i.endswith('FUSION'):
            df_origin[i] = df_origin[i].fillna(1)
        elif i=='MAGMA':
            df_origin[i] = df_origin[i].fillna(1)
        elif i=='mBATcombo':
            df_origin[i] = df_origin[i].fillna(1)
        elif i=='DistanceTSS':
            df_origin[i] = df_origin[i].fillna(df_origin[i].max())
        else:
            print(i)

# %%
# Split the dataframe
locus_indices = df_origin.index[~df_origin['snp_group'].isna()]
positive_indices = df_origin.index[df_origin['pharma_240507'] >= 0]
negative_indices = df_origin.index[(~df_origin.index.isin(locus_indices)) & (~df_origin.index.isin(positive_indices))]
# Sample 4000 random rows from the dataframe where 'pharma' < 0
negative_indices_sample = np.random.choice(negative_indices, size=1000000, replace=False)
# Combine the indices
sample_indices = np.concatenate([locus_indices, positive_indices, negative_indices_sample])
sample_indices = np.unique(sample_indices)

# %%
# Get the rows corresponding to these indices
df = df_origin.loc[sample_indices]

# %%
df = df.sample(frac=1).reset_index(drop=True)

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
feature_genetics = df.columns[np.r_[17:435]]
feature_genetics = [i for i in feature_genetics if not i.startswith('p_HEIDI')] 
feature_genetics = [i for i in feature_genetics if not i.endswith('DEPICT')]
feature_genetics = [i for i in feature_genetics if not i.endswith('OPERA')]
feature_genetics = [i for i in feature_genetics if ((not i.endswith('PPR')) and (not i.endswith('RWR')))]
feature_genetics = [i for i in feature_genetics if 'pops' not in i.lower()]
feature_genetics = [i for i in feature_genetics if not i.endswith('DEPICT')]
feature_genetics = [i for i in feature_genetics if not i.endswith('Network')]

feature_specific = [i for i in df.columns[np.r_[17:435]] if i not in feature_genetics]
feature_specific = [i for i in feature_specific if not i.startswith('p_HEIDI')]
feature_specific = [i for i in feature_specific if not i.endswith('DEPICT')]

# Valid features
gamma_example = pd.read_csv('/storage/yangjianLab/guoyazhou/GAMMA_git_data/GAMMA/feature/T2D_GAMMA.feature', sep='\t')
valid_list = gamma_example.isna().sum().sort_values()[:36].index
feature_genetics = [i for i in feature_genetics if i in valid_list] + list(df.columns[np.r_[432:435]])

feature_genetics = feature_genetics + ['GWAS']

feature_omim_clinvar_mgi = df.columns[np.r_[435:448]]

features_path_net = df.columns[448:3472]

feature_pretrain = df.columns[3472:3984]

feature_ppi_pca = df.columns[3984:4484]

feature_cancer_specific = df.columns[4484:4501]
feature_cancer_genefeature = df.columns[4501:4526]

# feature_all = feature_genetics
# feature_all = feature_genetics + feature_specific + list(feature_omim_clinvar_mgi) + list(feature_cancer_specific)
# feature_all = feature_specific + list(features_path_net) + list(feature_pretrain) + list(feature_ppi_pca)
# feature_all = list(features_path_net) + list(feature_pretrain) + list(feature_ppi_pca) + list(feature_cancer_genefeature)
feature_all = feature_genetics + feature_specific + list(feature_omim_clinvar_mgi) + list(features_path_net) + list(feature_pretrain) + list(feature_ppi_pca) # + list(feature_cancer_specific) + list(feature_cancer_genefeature)

print(f'n features: {len(feature_all)}')

# %%
feature_df = df[feature_all]

# Convert labels to one-hot encoding for multi-label classification
y_multilabel = np.array([y_dict[label] for label in df['pharma_240507']])
groups = df.chr
phase_list = ['pre_cli', 'phase_1', 'phase_2', 'phase_3', 'approve']
percent_top = [0.025,0.05,0.075,0.10,0.15,0.20,0.30,0.40,0.50,0.75,1]
combi = [['no_pr', -1, 0],
         ['pr_p1', 0, 1],
         ['pr_p2', 0, 2],
         ['pr_p3', 0, 3],
         ['pr_ap', 0, 4],
         ['p1_ap', 1, 4]]

# %%
# Split the data into train, validation, and test sets
group_kfold = GroupKFold(n_splits=22)
gss = GroupShuffleSplit(n_splits=1, test_size=1/5)

all_df_out = []
for train_val_index, test_index in group_kfold.split(df, y_multilabel, groups):
    train_index, val_index = next(gss.split(df.iloc[train_val_index], y_multilabel[train_val_index], groups.iloc[train_val_index]))

    train_index = train_val_index[train_index]
    val_index = train_val_index[val_index]

    scaler = StandardScaler()
    feature_train = feature_df.iloc[train_val_index]
    feature_train = scaler.fit_transform(feature_train)
    
    # Transform the validation and test sets
    feature_val = scaler.transform(feature_df.iloc[val_index])
    feature_test = scaler.transform(feature_df.iloc[test_index])
    
    assert len(df.iloc[test_index].chr.unique()) == 1
    test_chromosome = df.iloc[test_index].chr.unique()[0]
    
    output_prefix = '/storage/yangjianLab/sunshufeng/gamma_v7/results_250723/valid_epoch3_512+32_neg1000000_nn_not_weighted_0723'
    if os.path.exists(f'{output_prefix}_chr{test_chromosome}.csv'):
        print(f'File chr{test_chromosome} exists, skip')
        continue
    
    feature_output = scaler.transform(df_origin.loc[df_origin.chr == test_chromosome, feature_all])
    
    ### %%
    ##########################
    ### SETTINGS
    ##########################

    # Hyperparameters
    random_seed = 1
    learning_rate = 0.001
    num_epochs = 3
    batch_size = 2048

    # Architecture
    NUM_CLASSES = 5  # Multi-label classification with 5 classes

    # Other
    DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print('Training on', DEVICE)

    ### %%
    # Calculate class weights for each label separately
    def compute_multilabel_class_weights(y_multilabel):
        """Compute class weights for multi-label classification."""
        class_weights = []
        for i in range(y_multilabel.shape[1]):  # For each label
            # Get unique classes and their weights
            unique_classes = np.unique(y_multilabel[:, i])
            if len(unique_classes) > 1:  # Only if both classes exist
                weights = compute_class_weight('balanced', 
                                             classes=unique_classes, 
                                             y=y_multilabel[:, i])
                # Create weight dict
                weight_dict = dict(zip(unique_classes, weights))
                # For binary classification, we need the weight for positive class
                pos_weight = weight_dict.get(1, 1.0) / weight_dict.get(0, 1.0)
            else:
                pos_weight = 1.0
            class_weights.append(pos_weight)
        return torch.tensor(class_weights, dtype=torch.float32)

    # Compute class weights for the training set
    # pos_weights = compute_multilabel_class_weights(y_multilabel[train_val_index])
    # Not weighted version of the model
    pos_weights = torch.tensor([1.0, 1.0, 1.0, 1.0, 1.0], dtype=torch.float32)
    pos_weights = pos_weights.to(DEVICE)
    print(f"Positive class weights: {pos_weights}")

    ### %%
    # Prepare the dataset class
    class MyDataset(Dataset):
        def __init__(self, feature_array, label_array, dtype=np.float32):
            self.features = feature_array.astype(np.float32)
            self.labels = label_array.astype(np.float32)  # Multi-label targets should be float
        def __getitem__(self, index):
            inputs = self.features[index]
            label = self.labels[index]
            return inputs, label
        def __len__(self):
            return self.labels.shape[0]

    class OutDataset(Dataset):
        def __init__(self, feature_array):
            self.features = feature_array.astype(np.float32)
        def __len__(self):
            return len(self.features)
        def __getitem__(self, index):
            return self.features[index]
    
    ### %%
    # Create datasets
    train_dataset = MyDataset(feature_train, y_multilabel[train_val_index])
    val_dataset = MyDataset(feature_val, y_multilabel[val_index])
    test_dataset = MyDataset(feature_test, y_multilabel[test_index])
    
    output_dataset = OutDataset(feature_output)

    # Create dataloaders
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    
    output_loader = DataLoader(output_dataset, batch_size=batch_size, shuffle=False)

    # Checking the dataset
    for inputs, labels in train_loader:  
        print('Input batch dimensions:', inputs.shape)
        print('Input label dimensions:', labels.shape)
        break

    ### %%
    # Initialize the model - Standard Multi-Label Neural Network
    class MLP(torch.nn.Module):

        def __init__(self, in_features, num_classes, num_hidden_1=512, num_hidden_2=32, num_hidden_3=16):
            super().__init__()

            self.my_network = torch.nn.Sequential(

                # 1st hidden layer
                torch.nn.Linear(in_features, num_hidden_1, bias=False),
                torch.nn.LeakyReLU(),
                torch.nn.Dropout(0.2),
                torch.nn.BatchNorm1d(num_hidden_1),

                # 2nd hidden layer
                torch.nn.Linear(num_hidden_1, num_hidden_2, bias=False),
                torch.nn.LeakyReLU(),
                torch.nn.Dropout(0.2),
                torch.nn.BatchNorm1d(num_hidden_2),

                # # 3rd hidden layer
                # torch.nn.Linear(num_hidden_2, num_hidden_3, bias=False),
                # torch.nn.LeakyReLU(),
                # torch.nn.Dropout(0.2),
                # torch.nn.BatchNorm1d(num_hidden_3),
            )

            # Standard output layer for multi-label classification
            self.fc = torch.nn.Linear(num_hidden_2, num_classes)

        def forward(self, x):
            x = self.my_network(x)
            # Get logits and probabilities
            logits = self.fc(x)
            probas = torch.sigmoid(logits)
            return logits, probas

    # Custom evaluation function for multi-label classification
    def compute_multilabel_metrics(model, data_loader, device):
        model.eval()
        with torch.no_grad():
            all_targets = []
            all_predictions = []
            all_probas = []
            
            for features, targets in data_loader:
                features = features.to(device)
                targets = targets.to(device)
                
                logits, probas = model(features)
                predictions = (probas > 0.5).float()  # Threshold at 0.5
                
                all_targets.append(targets.cpu().numpy())
                all_predictions.append(predictions.cpu().numpy())
                all_probas.append(probas.cpu().numpy())
            
            all_targets = np.concatenate(all_targets, axis=0)
            all_predictions = np.concatenate(all_predictions, axis=0)
            all_probas = np.concatenate(all_probas, axis=0)
            
            # Calculate metrics
            from sklearn.metrics import hamming_loss, accuracy_score, f1_score
            
            # Hamming loss (fraction of incorrect labels)
            hamming = hamming_loss(all_targets, all_predictions)
            
            # Exact match accuracy (all labels must be correct)
            exact_match = accuracy_score(all_targets, all_predictions)
            
            # F1 score (macro average)
            f1_macro = f1_score(all_targets, all_predictions, average='macro', zero_division=0)
            
            return hamming, exact_match, f1_macro

    torch.manual_seed(random_seed)
    model = MLP(in_features=len(feature_all), num_classes=NUM_CLASSES)
    model.to(DEVICE)

    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    
    # Use Weighted Binary Cross Entropy Loss with class weights
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weights)
    
    # Early stopping parameters
    patience = 5
    best_val_loss = float('inf')
    epochs_no_improve = 0
    best_model_state = None

    ### %%
    # Train the model
    for epoch in range(num_epochs):
        model.train()
        for batch_idx, (features, class_labels) in enumerate(train_loader):
            
            features = features.to(DEVICE)
            class_labels = class_labels.to(DEVICE)
            logits, probas = model(features)

            # Weighted BCE loss for multi-label classification
            loss = criterion(logits, class_labels)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            ### LOGGING
            if not batch_idx % 200:
                print ('Epoch: %03d/%03d | Batch %03d/%03d | Loss: %.4f' 
                    %(epoch+1, num_epochs, batch_idx, 
                        len(train_loader), loss))

        # Evaluate on validation set
        val_hamming, val_exact_match, val_f1 = compute_multilabel_metrics(model, val_loader, DEVICE)
        print(f'Validation Hamming Loss: {val_hamming:.4f}, Exact Match: {val_exact_match:.4f}, F1 Macro: {val_f1:.4f}')

        # Early Stopping (using Hamming loss - lower is better)
        if val_hamming < best_val_loss:
            best_val_loss = val_hamming
            epochs_no_improve = 0
            best_model_state = model.state_dict()
        else:
            epochs_no_improve += 1
            if epochs_no_improve >= patience:
                print('Early stopping!')
                break
        
    # Load the best model state
    if best_model_state is not None:
        model.load_state_dict(best_model_state)
        
    train_hamming, train_exact_match, train_f1 = compute_multilabel_metrics(model, train_loader, DEVICE)
    val_hamming, val_exact_match, val_f1 = compute_multilabel_metrics(model, val_loader, DEVICE)
    test_hamming, test_exact_match, test_f1 = compute_multilabel_metrics(model, test_loader, DEVICE)

    print(f'Hamming Loss (train/val/test): {train_hamming:.4f} | {val_hamming:.4f} | {test_hamming:.4f}')
    print(f'Exact Match (train/val/test): {train_exact_match:.4f} | {val_exact_match:.4f} | {test_exact_match:.4f}')
    print(f'F1 Macro (train/val/test): {train_f1:.4f} | {val_f1:.4f} | {test_f1:.4f}')

    ### %%
    all_probas = []
    all_predicted_labels = []
    model.eval()

    with torch.no_grad():
        for features_batch in output_loader:
            features_batch = features_batch.to(DEVICE)
            logits, probas = model(features_batch)
            predicted_labels = (probas > 0.5).float()  # Binary threshold
            all_probas.append(probas.cpu())
            all_predicted_labels.append(predicted_labels.cpu())

    # Concatenate all predicted labels
    all_probas = torch.cat(all_probas, dim=0)
    all_predicted_labels = torch.cat(all_predicted_labels, dim=0)

    df_probs = pd.DataFrame(all_probas.numpy(), columns=["yhat_pre_cli", "yhat_phase_1", "yhat_phase_2", "yhat_phase_3", "yhat_approve"])
    df_predictions = pd.DataFrame(all_predicted_labels.numpy(), columns=["pred_pre_cli", "pred_phase_1", "pred_phase_2", "pred_phase_3", "pred_approve"])

    ### %%
    df_out = df_origin[df_origin.columns[:18]].loc[df_origin.chr == test_chromosome].reset_index(drop=True)

    for i in range(len(phase_list)):
        df_out['y_'+phase_list[i]] = np.array(df_out['pharma_240507'].map(y_dict).tolist())[:,i]

    df_out = pd.concat([df_out, df_probs, df_predictions], axis=1)
    df_out.to_csv(f'{output_prefix}_chr{test_chromosome}.csv', index=False)

# %%
# clear the memory
del model
torch.cuda.empty_cache()
del df_origin
del df

# concatenate all the results
all_df_out = []
for i in range(1,23):
    df = pd.read_csv(f'{output_prefix}_chr{i}.csv')
    all_df_out.append(df)

df_out = pd.concat(all_df_out)
df_out.to_csv(f'{output_prefix}.csv', index=False)

for i in range(1,23):
    os.remove(f'{output_prefix}_chr{i}.csv')