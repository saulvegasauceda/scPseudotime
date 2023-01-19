#!/bin/bash

#SBATCH --partition=xeon-g6-volta
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=saulv@mit.edu
#SBATCH --job-name=scVelo_modeling
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --mem=384000
#SBATCH --output=/home/gridsan/ssauceda/neuroTF_shared/mouse_cortex_development/scVelo/logs/scVelo_modeling.out

python ./run_scVelo.py
