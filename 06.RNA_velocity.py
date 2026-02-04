import scvelo as scv
import pandas as pd
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")

loom_data = scv.read('data/combined.loom', cache=False)
loom_data.obs

loom_data.obs = loom_data.obs.rename(index = lambda x: 'Nec_'+x.replace('x', '').split("_")[0]+'_'+x.replace('x', '').split(":")[1]+'-1' if re.search(pattern='[nN][eE][cC]1',string=x) else 'Normal_'+x.replace('x', '').split("_")[0]+'_'+x.replace('x', '').split(":")[1]+'-1' )
loom_data.obs.head()

meta_path = "../NEC/sc_final/data"
sample_obs = pd.read_csv(os.path.join(meta_path, "neu_cellID_obs.csv"))
cell_umap= pd.read_csv(os.path.join(meta_path, "neu_cell_embeddings.csv"), header=0, names=["CellID", "UMAP_1", "UMAP_2"])
cell_celltype = pd.read_csv(os.path.join(meta_path, "neu_cell_celltype.csv"), header=0, names=["CellID", "celltype"])

sample_one = loom_data[np.isin(loom_data.obs.index, sample_obs)]
sample_one.obs.head()

annData = sc.read_h5ad('/lustre/user/taowlab/mengfj/project/NEC/sc_final/data/round2/round2_NEC.neu.harmony_singlet.h5ad')

adata = scv.utils.merge(sample_one, annData)
sc.pl.umap(adata, color='Minor_clusters', legend_loc='on data')
os.chdir("../NEC/sc_final/figure/velocyto")

scv.pp.filter_and_normalize(adata, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=60, n_neighbors=45)

scv.tl.velocity(adata)
scv.tl.velocity_graph(adata, n_jobs=8)
adata.uns['Minor_clusters_colors'] = np.array(['#C95100', '#F9B858', '#FFDB6E', '#ff7f5c', '#FFAF84', '#FF9BCE', '#ffc5c5', '#db89a2', '#e0b794', '#d1c5c6'])
scv.pl.velocity_embedding_stream(adata, basis='umap', color="Minor_clusters", fontsize=5, figsize= (10, 10), dpi=300, use_raw=True, show=True, save='neu_scvelo_embedding_stream.svg',n_neighbors = 80)