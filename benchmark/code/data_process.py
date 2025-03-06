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

def ensembl_to_symbol(gene_list):
    import mygene
    mg = mygene.MyGeneInfo()
    return mg.querymany(gene_list, scopes='ensembl.gene', fields='symbol', as_dataframe=True,
                 species='human').reset_index().drop_duplicates(subset='query')['symbol'].fillna('0').tolist()


dataset = 'mouse_gene' # 'human_gene'

mapping_file_path = "../benchmark/data/mouse_human_mapping.json"
with open(mapping_file_path, "r") as file:
    mouse_human_mapping = json.load(file)


if dataset == 'human_gene':
    destination_file = '../benchmark/data/human_gene_v2.2.h5'
    metadata = pd.read_csv("../benchmark/data/GEO_HUMAN_100K_V2.csv")
    metadata['GSM_ID'] = [bytes(sample, 'utf-8') for sample in metadata['GSM_ID']]
    with h5py.File(destination_file, 'r') as f:
        num_samples = len(f['meta/samples/geo_accession'][:])
        ensembl_id = f['meta/genes/ensembl_gene_id'][:]
        symbol = f['meta/genes/symbol'][:]
        expression = f['data/expression'][:]
        geo_accession = f['meta/samples/geo_accession'][:]
        sample = f['meta/samples/sample'][:]
        series_id = f['meta/samples/series_id'][:]
        adata = ad.AnnData(X=expression.T, obs={'geo_accession': geo_accession, 'sample': sample, 'series_id': series_id}, 
        var={'ensembl_id': ensembl_id, 'symbol': symbol})
        adata.var['ensembl_id'] = adata.var['ensembl_id'].apply(lambda x: x.decode('utf-8'))
        adata.var['ensembl.gene'] = adata.var['ensembl_id']
        adata.var.index = adata.var['ensembl.gene']
        adata.obs['geo_accession'] = adata.obs['geo_accession'].astype('category')
        adata.obs['sample'] = adata.obs['sample'].astype('category')
        adata.obs['series_id'] = adata.obs['series_id'].astype('category')
        adata.var_names_make_unique()
        metadata_batch = metadata[metadata['GSM_ID'].isin(adata.obs['geo_accession'])]
        duplicates = metadata_batch[metadata_batch.duplicated(subset='GSM_ID', keep=False)]
        print(f"Number of duplicate GSM_IDs: {len(duplicates)}")
        metadata_batch = metadata_batch.drop_duplicates(subset='GSM_ID', keep='first')

        merged_gender = adata.obs.merge(metadata_batch[['GSM_ID', 'Gender']], 
                            left_on='geo_accession', 
                            right_on='GSM_ID', 
                            how='left')
        
        adata.obs['gender'] = merged_gender['Gender'].values

        merged_organ = adata.obs.merge(metadata_batch[['GSM_ID', 'Organ']], 
                            left_on='geo_accession', 
                            right_on='GSM_ID', 
                            how='left')
        adata.obs['organ_system'] = merged_organ['Organ'].values

        merged_disease = adata.obs.merge(metadata_batch[['GSM_ID', 'Disease']], 
                            left_on='geo_accession', 
                            right_on='GSM_ID', 
                            how='left')
        adata.obs['disease'] = merged_disease['Disease'].values

        merged_age = adata.obs.merge(metadata_batch[['GSM_ID', 'Age_Group']], 
                            left_on='geo_accession', 
                            right_on='GSM_ID', 
                            how='left')
        adata.obs['age'] = merged_age['Age_Group'].values
        adata.obs['gender'] = adata.obs['gender'].astype('category')
        adata.obs['organ_system'] = adata.obs['organ_system'].astype('category')
        adata.obs['disease'] = adata.obs['disease'].astype('category')
        adata.obs['age'] = adata.obs['age'].astype('category')
        del adata.var['ensembl.gene']
        adata = adata[adata.obs['geo_accession'].isin(metadata['GSM_ID'])]
    
        metadata_train = pd.read_csv("../benchmark/data/GEO_HUMAN_100K_training.csv")
        train_categories = set(metadata_train['GSM_ID'])
        train_categories = np.array(list(train_categories))
        adata.obs['split'] = 'test'
        adata.obs['split'] = np.where(
        adata.obs['sample'].isin(train_categories), 'train', 'test')
        adata.write('../benchmark/data/human_gene_v2.2_split.h5ad')
        adata = sc.read_h5ad('../benchmark/data/human_gene_v2.2_split.h5ad')
        adata.obs['split'] = np.where(
        adata.obs['sample'].isin(train_categories), 'train', 'test')
        adata.write('../benchmark/data/human_gene_v2.2_split.h5ad')


elif dataset == 'mouse_gene':
    destination_file = '../benchmark/data/mouse_gene_v2.2.h5'
    metadata = pd.read_csv("../benchmark/data/GEO_MOUSE_100K.csv")
    metadata['GSM_ID'] = [bytes(sample, 'utf-8') for sample in metadata['GSM_ID']]
    with h5py.File(destination_file, 'r') as f:
        num_samples = len(f['meta/samples/geo_accession'][:])
        ensembl_id = f['meta/genes/ensembl_gene_id'][:]
        symbol = f['meta/genes/symbol'][:]
        expression = f['data/expression'][:]
        geo_accession = f['meta/samples/geo_accession'][:]
        sample = f['meta/samples/sample'][:]
        series_id = f['meta/samples/series_id'][:]
        adata = ad.AnnData(X=expression.T, obs={'geo_accession': geo_accession, 'sample': sample, 'series_id': series_id}, 
        var={'ensembl_id': ensembl_id, 'symbol': symbol})
        adata.var['ensembl_id'] = adata.var['ensembl_id'].apply(lambda x: x.decode('utf-8'))
        adata.var['ensembl.gene'] = adata.var['ensembl_id']
        adata.var.index = adata.var['ensembl.gene']
        adata.obs['geo_accession'] = adata.obs['geo_accession'].astype('category')
        adata.obs['sample'] = adata.obs['sample'].astype('category')
        adata.obs['series_id'] = adata.obs['series_id'].astype('category')
        adata.var_names_make_unique()
        metadata_batch = metadata[metadata['GSM_ID'].isin(adata.obs['geo_accession'])]
        duplicates = metadata_batch[metadata_batch.duplicated(subset='GSM_ID', keep=False)]
        print(f"Number of duplicate GSM_IDs: {len(duplicates)}")
        metadata_batch = metadata_batch.drop_duplicates(subset='GSM_ID', keep='first')

        merged_gender = adata.obs.merge(metadata_batch[['GSM_ID', 'Gender']], 
                            left_on='geo_accession', 
                            right_on='GSM_ID', 
                            how='left')
        
        adata.obs['gender'] = merged_gender['Gender'].values

        merged_organ = adata.obs.merge(metadata_batch[['GSM_ID', 'Organ']], 
                            left_on='geo_accession', 
                            right_on='GSM_ID', 
                            how='left')
        adata.obs['organ_system'] = merged_organ['Organ'].values

        merged_disease = adata.obs.merge(metadata_batch[['GSM_ID', 'Disease']], 
                            left_on='geo_accession', 
                            right_on='GSM_ID', 
                            how='left')
        adata.obs['disease'] = merged_disease['Disease'].values

        merged_age = adata.obs.merge(metadata_batch[['GSM_ID', 'Age_Group']], 
                            left_on='geo_accession', 
                            right_on='GSM_ID', 
                            how='left')
        adata.obs['age'] = merged_age['Age_Group'].values
        adata.obs['gender'] = adata.obs['gender'].astype('category')
        adata.obs['organ_system'] = adata.obs['organ_system'].astype('category')
        adata.obs['disease'] = adata.obs['disease'].astype('category')
        adata.obs['age'] = adata.obs['age'].astype('category')
        del adata.var['ensembl.gene']
        adata.write('../benchmark/data/mouse_gene_v2.2.h5ad')
        adata = adata[adata.obs['geo_accession'].isin(metadata['GSM_ID'])]
    
        metadata_train = pd.read_csv("../benchmark/data/GEO_MOUSE_100K_training.csv")
        train_categories = set(metadata_train['GSM_ID'])
        train_categories = np.array(list(train_categories))
        adata.obs['split'] = 'test'
        adata.obs['split'] = np.where(
        adata.obs['sample'].isin(train_categories), 'train', 'test')
        adata.write('../benchmark/data/mouse_gene_v2.2_split.h5ad')
        adata = sc.read_h5ad('../benchmark/data/mouse_gene_v2.2_split.h5ad')
        adata.obs['split'] = np.where(
        adata.obs['sample'].isin(train_categories), 'train', 'test')
        adata.write('../benchmark/data/mouse_gene_v2.2_split.h5ad')
        adata = sc.read_h5ad('../benchmark/data/mouse_gene_v2.2_split.h5ad')
        adata = adata[:, adata.var.index.isin(mouse_human_mapping.keys())]
        adata.var.index = adata.var.index.map(mouse_human_mapping)
        adata.var_names_make_unique()
        adata.var['symbol_human'] = ensembl_to_symbol(adata.var.index.tolist())
        adata.write_h5ad('../benchmark/data/mouse_gene_v2.2_split_converted.h5ad')
