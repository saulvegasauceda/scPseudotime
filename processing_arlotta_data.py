import anndata
import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sbn
import scipy
import matplotlib.pyplot as plt

path_to_umap = "/home/gridsan/ssauceda/neuroTF_shared/mouse_cortex_development/umap/"
umap_coords_file = path_to_umap + "cluster_scDevSC.merged.umap.txt"
adata_file = path_to_umap + "reduced_counts_adata.h5ad"

adata = anndata.read_h5ad(adata_file)
umap_info = pd.read_csv(umap_coords_file, header=None, sep='\t')

#====== Adding UMAP coordinates ======#

axes = list(umap_info.iloc[0, 1:])
barcodes = list(umap_info.iloc[2:, 0])
umap_info = umap_info.iloc[2:, 1:].astype(float)
umap_info.columns = axes
umap_info.index = barcodes

# dealing with barcode discrepancy between Arlotta obj & our obj
our_barcodes_ordered = adata.obs.index
matching_barcodes = umap_info.index.isin(our_barcodes_ordered)
umap_info = umap_info[matching_barcodes]
umap_info = umap_info.reindex(our_barcodes_ordered)
umap_info = np.array(umap_info)

adata.obsm["X_umap"] = umap_info

#====== Filtering to excitatory cells only ======#

# only interested in analyzing excitatory cell types
excitatory = ['Excitatory neurons', 'Apical progenitors', 'Intermediate progenitors']
barcodes_excitatory = adata.obs['Gral_cellType'].isin(excitatory)
adata = adata[barcodes_excitatory]

#====== Running scVelo ======#

sum_spliced_features = np.array(adata.layers["spliced"].sum(axis=1))
sum_unspliced_features = np.array(adata.layers["unspliced"].sum(axis=1))

# filter out cells with less than 500 spliced or unspliced features
cells_w_enough_splice_info = np.asarray((sum_spliced_features + sum_unspliced_features) >= 500)

adata = adata[cells_w_enough_splice_info, :]

scv.pp.normalize_per_cell(adata, counts_per_cell_after=1e4, enforce=True)
sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=3000)
scv.pp.log1p(adata)

regress_out = ['nCount_RNA', 'nFeature_RNA', 'percent_mito', 'CC_Difference']
sc.pp.regress_out(adata, keys=regress_out)

scv.pp.moments(adata, n_pcs=30, n_neighbors=30, use_highly_variable=True)
scv.tl.recover_dynamics(adata, n_jobs=20)
scv.tl.velocity(adata, mode='dynamical', use_highly_variable=True)
scv.tl.velocity_graph(adata, n_jobs=20)

scv.tl.latent_time(adata)
scv.tl.velocity_pseudotime(adata)
# this is needed due to a current bug
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups="New_cellType")
