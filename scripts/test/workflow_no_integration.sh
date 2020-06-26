#! /bin/bash -l

########################
### DEFINE VARIABLES ###
########################
var_to_plot='dataset,mouse_nr,day_post_infection,organ,infection'
var_to_regress='nFeature_RNA,percent_mito,S.Score,G2M.Score'
# main='/home/jonrob/projects/d_angeletti_1910'
# script_path='/home/jonrob/repos/sauron/scripts'
main='/cephyr/users/jonrob/Hebbe/projects/d_angeletti_1910'
script_path='/cephyr/users/jonrob/Hebbe/repos/sauron/scripts'
# main='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'
# script_path='/Users/jonrob/Documents/NBIS/repos/sauron/scripts'

cd $main
mkdir analysis
mkdir log


##################################
### ACTIVATE CONDA ENVIRONMENT ###
##################################
# conda activate sauron-env  # macOS (local)
source activate Sauron.v1  # linux/unix (cluster)


# Start from Seurat object in 01_qc from original workflow


##################################################
### RUN INTEGRATION SCRIPT WITHOUT INTEGRATING ###
##################################################
mkdir $main/'analysis/02_cluster_noint'
Rscript $script_path/02_integrate.R \
--Seurat_object_path $main/'analysis/01_qc/filt_seurat_object.rds' \
--var_genes 'seurat' \
--integration_method 'NONE' \
--cluster_use 'all' \
--assay 'RNA' \
--output_path $main/'analysis/02_cluster_noint' \
2>&1 | tee $main/'analysis/02_cluster_noint/02_integrate_log.txt'



###################################################
### RUN DIMENSIONALITY REDUCTION AND CLUSTERING ###
###################################################
Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $main/'analysis/02_cluster_noint/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'var,1' \
--var_genes 'seurat' \
--dim_reduct_use 'umap' \
--dim_reduct_params 'umap, n.neighbors=30, min.dist=0.01, spread=3, n.epochs=500, learning.rate=0.5; umap10, n.neighbors=30, min.dist=0.01, n.epochs=500' \
--pre_dim_reduct 'pca' \
--cluster_use 'all' \
--cluster_method 'louvain,HC' \
--assay 'RNA' \
--output_path $main/'analysis/02_cluster_noint' \
2>&1 | tee $main/'analysis/02_cluster_noint/03_dr_and_cluster_log.txt'




conda deactivate
