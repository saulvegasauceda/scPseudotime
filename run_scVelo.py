import anndata
import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import scipy

#====== Loading anndata objects ======#
path_to_scVelo = "/home/gridsan/ssauceda/neuroTF_shared/mouse_cortex_development/scVelo/"

# Input files
file_imputed_counts = "arlotta.imputed.high_var_only.regress_out.umap.filter_splice.h5ad"
file_norm_counts = "arlotta.norm.high_var_only.regress_out.umap.filter_splice.h5ad"

input_files = [file_imputed_counts, file_norm_counts]

#====== Running scVelo ======#
for file in input_files:
	print("Reading :", file)
	adata = anndata.read_h5ad(path_to_scVelo + file)
	output_file = file[: -4] + "scvelo." +  file[-4:]
	# Running modeling
	scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
	scv.tl.recover_dynamics(adata, n_jobs=20)
	scv.tl.velocity(adata, mode='dynamical')
	scv.tl.velocity_graph(adata, n_jobs=20)

	# Run pseudotime algorithms
	scv.tl.latent_time(adata)
	scv.tl.velocity_pseudotime(adata)
	# this is needed due to a current bug
	adata.uns['neighbors']['distances'] = adata.obsp['distances']
	adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
	scv.tl.paga(adata, groups="New_cellType")

	# Exporting h5ad
	print("Saving :", output_file)
	adata.write_h5ad(path_to_scVelo + output_file)
	print()
