import warnings
warnings.filterwarnings("ignore")
import os
import logging
import random
import torch
import sys

import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix
import scanpy as sc
import matplotlib.pyplot as plt
import json
import pandas as pd
import h5py
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report, confusion_matrix


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
DEVICE = device
dataset = 'human_gene' # 'mouse_gene'
method = 'CellPLM' # 'CellPLM_finetune', 'geneformer', 'scGPT', 'scvi'
classify_type = 'disease' # 'gender', 'age', 'organ_system'
num_epochs = 600


def set_seed(rndseed, cuda: bool = True, extreme_mode: bool = False):
    os.environ["PYTHONHASHSEED"] = str(rndseed)
    random.seed(rndseed)
    np.random.seed(rndseed)
    torch.manual_seed(rndseed)
    if cuda:
        torch.cuda.manual_seed(rndseed)
        torch.cuda.manual_seed_all(rndseed)
    if extreme_mode:
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
    # dgl.seed(rndseed)
    # dgl.random.seed(rndseed)
    logging.info(f"Setting global random seed to {rndseed}")

set_seed(42)


if dataset == 'human_gene':
    adata = sc.read_h5ad('./data/data/human_gene_v2.2_split.h5ad')

if dataset == 'mouse_gene':
    adata = sc.read_h5ad('./data/data/mouse_gene_v2.2_split_converted.h5ad')
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=200)
    metadata_train = pd.read_csv("./data/data/GEO_MOUSE_100K_trainin_V2.csv")
    metadata_test = pd.read_csv("./data/data/GEO_MOUSE_100K_test_V2.csv")
    metadata = pd.concat([metadata_train, metadata_test], ignore_index=True)
    merged_age = adata.obs.merge(metadata[['GSM_ID', 'Age_Group_New']], 
                                left_on='geo_accession', 
                                right_on='GSM_ID', 
                                how='left')
    adata.obs['age'] = merged_age['Age_Group_New'].values                            

if method == 'geneformer':
    embeddings = torch.load(f'./data/data/geneformer/human_embs.pt')
    index = pd.Index(list(range(512)))
    filtered_index = index[index.map(lambda x: isinstance(x, int))]
    embeddings = np.array(embeddings[filtered_index])
    adata.obsm['emb'] = embeddings
    


adata_train = adata[adata.obs['split'] == 'train']  
adata_test = adata[adata.obs['split'] == 'test']
    
# adata.obs['celltype'] = adata.obs['gender']
if classify_type == 'disease':
    adata = adata[adata.obs['disease'].notnull()]
    adata = adata[adata.obs['disease'] != 'No Disease Mentioned']
    disease_train = set(adata_train.obs['disease'])
    disease_test = set(adata_test.obs['disease'])
    overlap = disease_train & disease_test
    adata = adata[adata.obs['disease'].isin(overlap)]
elif classify_type == 'organ_system':
    adata = adata[adata.obs['organ_system'].notnull()]
    organ_system_train = set(adata_train.obs['organ_system'])
    organ_system_test = set(adata_test.obs['organ_system'])
    overlap = organ_system_train & organ_system_test
    adata = adata[adata.obs['organ_system'].isin(overlap)]
elif classify_type == 'gender':
    adata = adata[adata.obs['gender'].notnull()]
    adata = adata[adata.obs['gender'] != 'Unknown',:]
elif classify_type == 'age':
    adata = adata[adata.obs['age'].notnull()]
adata.var_names_make_unique()


if method == 'CellPLM':
    PRETRAIN_VERSION = "20231027_85M"
    from CellPLM.utils import set_seed
    from CellPLM.pipeline.cell_embedding import CellEmbeddingPipeline
    pipeline = CellEmbeddingPipeline(pretrain_prefix=PRETRAIN_VERSION, # Specify the pretrain checkpoint to load
                                    pretrain_directory='./data/data/CellPLM/ckpt/')


    embedding = pipeline.predict(adata, # An AnnData object
                    # covariate_fields=["species"],  
                    inference_config={"batch_size": 4000},           
                    device=DEVICE) # Specify a gpu or cpu for model inference
    
if method == 'CellPLM_finetune':
    PRETRAIN_VERSION = "20231027_85M"
    import sys
    from scipy.sparse import csr_matrix
    from CellPLM.utils import set_seed
    from CellPLM.pipeline.cell_type_annotation import CellTypeAnnotationPipeline, CellTypeAnnotationDefaultPipelineConfig, CellTypeAnnotationDefaultModelConfig
    pipeline_config = CellTypeAnnotationDefaultPipelineConfig.copy()
    model_config = CellTypeAnnotationDefaultModelConfig.copy()
    model_config['out_dim'] = adata.obs[classify_type].nunique()
    pipeline = CellTypeAnnotationPipeline(pretrain_prefix=PRETRAIN_VERSION, # Specify the pretrain checkpoint to load
                                      overwrite_config=model_config,  # This is for overwriting part of the pretrain config
                                      pretrain_directory='./data/data/CellPLM/ckpt/')
    pipeline.fit_no_valid(adata, # An AnnData object
            pipeline_config, # The config dictionary we created previously, optional
            split_field = 'split', #  Specify a column in .obs that contains split information
            train_split = 'train',
            # valid_split = 'valid',
            label_fields = [classify_type], # Specify a column in .obs that contains cell type labels
            device = device) 
    y_pred_labels = pipeline.predict(
                adata, # An AnnData object
                pipeline_config, # The config dictionary we created previously, optional
                device = device
            )
    y_pred_labels = y_pred_labels.cpu().numpy()
    y_pred_labels = y_pred_labels[adata.obs['split'] == 'test']
    label_encoder = LabelEncoder()

    y_test = label_encoder.fit_transform(adata[adata.obs['split'] == 'test'].obs[classify_type])
    y_test_original = adata[adata.obs['split'] == 'test'].obs[classify_type]

    
    # Calculate accuracy
    accuracy = accuracy_score(y_test, y_pred_labels)
    print(f"Test Accuracy for {classify_type}: {accuracy:.4f}")

    class_names = label_encoder.classes_
    y_pred_decoded = label_encoder.inverse_transform(y_pred_labels)

    print(f"\nPer-category Accuracy for {classify_type}:")
    accuracy_list = []
    for category in class_names:
        mask = y_test_original == category
        category_pred = y_pred_decoded[mask]
        category_true = y_test_original[mask]
        category_accuracy = accuracy_score(category_true, category_pred)
        category_samples = sum(mask)
        accuracy_list.append(category_accuracy)
        print(f"{category}: {category_accuracy:.4f} (n={category_samples})")
    print("mean of accuracy_list:", np.mean(accuracy_list))
    print("std of accuracy_list:", np.std(accuracy_list))

    unique_labels = y_test_original.dropna().unique()

    report = classification_report(y_test_original, y_pred_decoded, labels=unique_labels, target_names=unique_labels, output_dict=True)

    report_df = pd.DataFrame(report).transpose()

    print(report_df)
    report_df.to_csv('{}_{}_{}_classification_report.csv'.format(dataset, classify_type, method))
    print("mean of report_df['f1-score']:", report_df['f1-score'].mean())
    print("std of report_df['f1-score']:", report_df['f1-score'].std())
    print("mean of report_df['precision']:", report_df['precision'].mean())
    print("std of report_df['precision']:", report_df['precision'].std())
    print("mean of report_df['recall']:", report_df['recall'].mean())
    print("std of report_df['recall']:", report_df['recall'].std())

    print("adata.shape:", adata.shape)
    print("Number of categories in adata.obs['gender']:", adata.obs['gender'].nunique())
    print("Number of categories in adata.obs['age']:", adata.obs['age'].nunique())
    print("Number of categories in adata.obs['organ_system']:", adata.obs['organ_system'].nunique())
    print("Number of categories in adata.obs['disease']:", adata.obs['disease'].nunique())
    print("Number of categories in adata.obs['series_id']:", adata.obs['series_id'].nunique())
    sys.exit() 

elif method == 'scvi':
    import scvi
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
      adata,
      n_top_genes=2000,
      subset=True)
    
    scvi.model.SCVI.setup_anndata(adata, batch_key="series_id")
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=16, gene_likelihood="nb")
    model.train()
    embedding = model.get_latent_representation()

elif method == 'scGPT':
    adata.X = adata.X.astype('int32')
    from pathlib import Path
    import warnings

    import scanpy as sc
    import scib
    import numpy as np
    import sys

    sys.path.insert(0, "../")

    import scgpt as scg
    import matplotlib.pyplot as plt
    import anndata

    plt.style.context('default')
    warnings.simplefilter("ignore", ResourceWarning)

    model_dir = Path("./data/data/scGPT/ckpt/scGPT_human")

    if dataset == 'human_gene':
        gene_col = "symbol"
    elif dataset == 'mouse_gene':   
        gene_col = "symbol_human"
    embed_adata = scg.tasks.embed_data(
    adata,
    model_dir,
    gene_col=gene_col,
    batch_size=64,
)
    embedding = embed_adata.obsm["X_scGPT"]


import matplotlib.pyplot as plt


if method == 'scvi':
    adata.obsm['emb'] = embedding
elif method == 'CellPLM':
    adata.obsm['emb'] = embedding.cpu().numpy()
elif method == 'raw':
    adata.obsm['emb'] = adata.X
elif method == 'scGPT':
    adata.obsm['emb'] = embedding   


import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

class MLPClassifier(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(MLPClassifier, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim)
        )
    
    def forward(self, x):
        return self.model(x)
latent_embeddings = adata.obsm['emb']

label_encoder = LabelEncoder()
labels = adata.obs[classify_type]
labels_encoded = label_encoder.fit_transform(labels)

# Split data into training and testing sets
X_train = adata[adata.obs['split'] == 'train'].obsm['emb']
X_test = adata[adata.obs['split'] == 'test'].obsm['emb']
y_train = label_encoder.fit_transform(adata[adata.obs['split'] == 'train'].obs[classify_type])
y_test = label_encoder.fit_transform(adata[adata.obs['split'] == 'test'].obs[classify_type])
y_test_original = adata[adata.obs['split'] == 'test'].obs[classify_type]

 
X_train_tensor = torch.tensor(X_train, dtype=torch.float32).to(device)
X_test_tensor = torch.tensor(X_test, dtype=torch.float32).to(device)
y_train_tensor = torch.tensor(y_train, dtype=torch.long).to(device)
y_test_tensor = torch.tensor(y_test, dtype=torch.long).to(device)

input_dim = latent_embeddings.shape[1]
hidden_dim = 256  # You can adjust this
output_dim = len(np.unique(labels_encoded))

model = MLPClassifier(input_dim, hidden_dim, output_dim).to(device)
criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

batch_size = 32

for epoch in range(num_epochs):
    model.train()
    permutation = torch.randperm(X_train_tensor.size(0))
    
    epoch_loss = 0.0
    for i in range(0, X_train_tensor.size(0), batch_size):
        indices = permutation[i:i + batch_size]
        batch_X, batch_y = X_train_tensor[indices], y_train_tensor[indices]

        # Forward pass
        outputs = model(batch_X)
        loss = criterion(outputs, batch_y)

        # Backward pass and optimization
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        epoch_loss += loss.item()
    
    print(f"Epoch [{epoch+1}/{num_epochs}], Loss: {epoch_loss:.4f}")

# Evaluate the model
model.eval()
with torch.no_grad():
    y_pred = model(X_test_tensor)
    y_pred_labels = torch.argmax(y_pred, axis=1).cpu().numpy()

accuracy = accuracy_score(y_test, y_pred_labels)
print(f"Test Accuracy for {classify_type}: {accuracy:.4f}")

class_names = label_encoder.classes_
y_pred_decoded = label_encoder.inverse_transform(y_pred_labels)

print(f"\nPer-category Accuracy for {classify_type}:")
accuracy_list = []
for category in class_names:
    mask = y_test_original == category
    category_pred = y_pred_decoded[mask]
    category_true = y_test_original[mask]
    category_accuracy = accuracy_score(category_true, category_pred)
    category_samples = sum(mask)
    accuracy_list.append(category_accuracy)
    print(f"{category}: {category_accuracy:.4f} (n={category_samples})")
print("mean of accuracy_list:", np.mean(accuracy_list))
print("std of accuracy_list:", np.std(accuracy_list))



unique_labels = y_test_original.dropna().unique()

report = classification_report(y_test_original, y_pred_decoded, labels=unique_labels, target_names=unique_labels, output_dict=True)

report_df = pd.DataFrame(report).transpose()

print(report_df)
report_df.to_csv('./data/results/{}_{}_{}_classification_report.csv'.format(dataset, classify_type, method))
print("mean of report_df['f1-score']:", report_df['f1-score'].mean())
print("std of report_df['f1-score']:", report_df['f1-score'].std())
print("mean of report_df['precision']:", report_df['precision'].mean())
print("std of report_df['precision']:", report_df['precision'].std())
print("mean of report_df['recall']:", report_df['recall'].mean())
print("std of report_df['recall']:", report_df['recall'].std())



print("adata.shape:", adata.shape)
print("Number of categories in adata.obs['gender']:", adata.obs['gender'].nunique())
print("Number of categories in adata.obs['age']:", adata.obs['age'].nunique())
print("Number of categories in adata.obs['organ_system']:", adata.obs['organ_system'].nunique())
print("Number of categories in adata.obs['disease']:", adata.obs['disease'].nunique())
print("Number of categories in adata.obs['series_id']:", adata.obs['series_id'].nunique())