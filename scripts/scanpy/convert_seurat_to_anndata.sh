#! /bin/bash -l

########################
### DEFINE VARIABLES ###
########################
# main='/home/jonrob/projects/d_angeletti_1910'
# main='/cephyr/users/jonrob/Hebbe/projects/d_angeletti_1910'
main='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'

cd $main


# Activate conda environment
conda activate Sauron.v1  # macOS (local)
# source activate Sauron.v1  # linux/unix (cluster)


###############################################
### CONVERT SEURAT OBJECT TO ANNDATA OBJECT ###
###############################################
mkdir -p $main'/analysis/trajectory_paga'
Rscript /Users/jonrob/Documents/NBIS/repos/single-cell-hackathon-2020/scanpy/seurat2loom.R \
-i $main'/analysis/immcantation/seurat_object_VDJannot.rds' \
-o $main'/analysis/trajectory_paga/seurat_object_VDJannot.loom'

# activate scanpy environment
conda deactivate
conda activate hackathon_scanpy

# convert loom file to h5ad file (anndata object) for use in scanpy
python /Users/jonrob/Documents/NBIS/repos/single-cell-hackathon-2020/scanpy/loom2anndata.py \
-i $main'/analysis/trajectory_paga/seurat_object_VDJannot.loom' \
-o $main'/analysis/trajectory_paga/anndata_object_VDJannot.h5ad'




