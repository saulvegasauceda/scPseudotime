#!/bin/bash

#SBATCH --partition=xeon-p8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=saulv@mit.edu
#SBATCH --job-name=magic_pseudotemp
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --mem=192000
#SBATCH --output=/home/gridsan/ssauceda/neuroTF_shared/mouse_cortex_development/magic/logs/magic_pseudotemp.out

python ./run_magic.py -i ../umap/reduced_counts_adata.h5ad -o ./ -n reduced_counts_adata --n_jobs -1
