#! /bin/bash -l

# Converts a Seurat object (.rds file) into a .loom file and
# AnnData object (.h5ad file) in the same directory

# Specify paths
main='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'
seurat_object=$main'/analysis/06_cluster/seurat_object.rds'

cd $main

# activate conda environment
source activate Sauron.v1

# Convert Seurat .rds object to .loom
loom_object="${seurat_object//.rds/.loom}"
Rscript $main'/scripts/scanpy/seurat2loom.R' -i $seurat_object -o $loom_object

# change conda environment to scanpy-env
conda deactivate
source activate scanpy-env

# convert loom file to h5ad file (anndata object) for use in scanpy
ann_object="${loom_object//.loom/.h5ad}"
python $main'/scripts/scanpy/loom2anndata.py' -i $loom_object -o $ann_object




