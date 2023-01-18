import argparse

# Argument parser -------------------------------------------------------------
epilog = 'Runs MAGIC from https://github.com/KrishnaswamyLab/MAGIC to  ' + \
    'impute the absolute counts of the scRNA-seq data. We reccomend looking at ' + \
    'the following: https://magic.readthedocs.io/en/stable/  ' + \
    'There you will find the parameters for the magic operator under class magic.MAGIC'

parser = argparse.ArgumentParser(epilog=epilog)
parser.add_argument('-i','--input',
                    type=str,
                    dest='input',
                    help='input directory path to the h5ad file of RNA counts')
parser.add_argument('-o','--output',
                    type=str,
                    dest='output',
                    help='output directory path for the imputed RNA counts')
parser.add_argument('-n','--name',
                    type=str,
                    dest='sample_name',
                    default='sample',
                    help='sample name')
parser.add_argument('--knn',
                    type=int,
                    dest='knn',
                    default=5,
                    help='')
parser.add_argument('--knn_max',
                    type=int,
                    dest='knn_max',
                    default=None,
                    help='')
parser.add_argument('--decay',
                    type=int,
                    dest='decay',
                    default=1,
                    help='')
parser.add_argument('--t',
                    type=int,
                    dest='t',
                    default=3,
                    help='')
parser.add_argument('--n_pca',
                    type=int,
                    dest='n_pca',
                    default=100,
                    help='')
parser.add_argument('--solver',
                    type=str,
                    dest='solver',
                    default='exact',
                    help='')
parser.add_argument('--knn_dist',
                    type=str,
                    dest='knn_dist',
                    default='euclidean',
                    help='')
parser.add_argument('--n_jobs',
                    type=int,
                    dest='n_jobs',
                    default=1,
                    help='')
parser.add_argument('--random_state',
                    type=int,
                    dest='random_state',
                    default=None,
                    help='')
parser.add_argument('--verbose',
                    type=bool,
                    dest='verbose',
                    default=True,
                    help='')

args = parser.parse_args()

# Run MAGIC --------------------------------------------------------------------

import scanpy as sc
import scanpy.external as sce
import magic

print(f"Loading input file ...")
adata = sc.read_h5ad(args.input)
adata.var_names_make_unique()

print("Normalizing counts ...")
sc.pp.normalize_total(adata)
sc.pp.sqrt(adata)

print("Running MAGIC ...")
sce.pp.magic(adata = adata, knn = args.knn, knn_max = args.knn_max, decay = args.decay, 
             t = args.t, n_pca = args.n_pca, solver = args.solver, 
             knn_dist = args.knn_dist, n_jobs = args.n_jobs, 
             random_state = args.random_state, verbose = args.verbose)


adata.X = adata.X**2

print("Saving output as h5ad file...")
output_dir = args.output if args.output[-1] == "/" else (args.output + "/")
suffix = "magic.h5ad"
file_name = args.sample_name + '.' + suffix if args.sample_name else suffix
output_file = output_dir + file_name
adata.write_h5ad(output_file)

print("All done...")

