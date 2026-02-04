### Tacco
import os
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import anndata as ad
import scanpy as sc

import tacco as tc

output_path = '../NEC1/tacco'

### load st_adata
adata_st = sc.read_h5ad('../NEC1/data/cellbins_scanpy.h5ad')
adata_st.var['SYMBOL'] = adata_st.var['real_gene_name']
# duplicated_genes = adata_st.var['SYMBOL'][adata_st.var['SYMBOL'].duplicated()]

genes_to_keep = adata_st.var.index[~adata_st.var.index.isin(remove_genes)]
adata_st=adata_st[:,genes_to_keep].copy()
adata_st.var_names = adata_st.var['SYMBOL']
adata_st.var_names.name = None

from scipy.sparse import csr_matrix
adata_st.X = adata_st.X.toarray()
adata_st.X = adata_st.X.astype(float)
sc.pp.calculate_qc_metrics(adata_st, inplace=True)
adata_st.X = csr_matrix(adata_st.X)
adata_st.var['MT'] = [gene.startswith('MT-') for gene in adata_st.var['SYMBOL']]
adata_st.obs['MT_frac'] = adata_st[:, adata_st.var['MT'].tolist()].X.sum(1).A.squeeze() / adata_st.obs['total_counts']

adata_st.obs["sample"] = 'NEC1'
adata_st.obs_names = adata_st.obs["sample"] \
                  + '_' + adata_st.obs_names
adata_st.obs.index.name = 'spot_id'

adata_st.obsm['MT'] = adata_st[:, adata_st.var['MT'].values].X.toarray()
adata_st = adata_st[:, ~adata_st.var['MT'].values]

st_adata = adata_st


### load sc_adata
import scipy.sparse as sp
sc_adata = sc.read_text("../summary/scRNA-seq_reference/NEC/scRNA_counts.csv", delimiter='\t', first_column_names=True)
sc_adata.X = sp.csr_matrix(sc_adata.X)
metadata = pd.read_csv("../summary/scRNA-seq_reference/NEC/metadata.csv", sep = '\t', header = 0, index_col = 0)
sc_adata.obs = metadata
r_embedding = pd.read_csv("../summary/scRNA-seq_reference/NEC/scRNA_embedding.csv", sep = '\t', header = 0, index_col = 0)
sc_adata.obsm["X_umap"] = r_embedding.values
sc_adata.var['SYMBOL'] = sc_adata.var.index


### annotate Major_clusters
tc.tl.annotate(st_adata, sc_adata,'Major_clusters', result_key = 'Major_clusters')
st_adata.obsm['Major_clusters']['label'] = st_adata.obsm['Major_clusters'].idxmax(axis=1)
st_adata.obsm['Major_clusters']['label'].value_counts()

st_adata.obsm['Major_clusters'].to_csv(output_path + '/Tacco_result_maj.csv')