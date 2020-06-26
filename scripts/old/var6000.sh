
var_to_plot='dataset,mouse_nr,day_post_infection,organ,infection'
var_to_regress='nCount_RNA,nFeature_RNA,perc_mito,CC.Diff'  # regress CC.Diff (S - G2M) instead of S and G2M separately

# main='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'
# script_path='/Users/jonrob/Documents/NBIS/repos/sauron/scripts'
main='/cephyr/users/jonrob/Hebbe/projects/d_angeletti_1910'
script_path='/cephyr/users/jonrob/Hebbe/repos/sauron/scripts'


cd $main

source activate Sauron.v1


# ###################################################
# ### RUN INTEGRATION WITH PRE-SCALING/REGRESSION ###
# ###################################################
# mkdir $main/'analysis/06_cluster_var6000'
# Rscript $main/scripts/02b_scale_integrate.R \
# --Seurat_object_path $main/'analysis/05_qc/filt_seurat_object.rds' \
# --sauron_script_path $script_path \
# --regress $var_to_regress \
# --integration_method 'mnn,dataset' \
# --var_genes 'seurat' \
# --nSeurat '6000' \
# --cluster_use 'all' \
# --assay 'RNA' \
# --output_path $main/'analysis/06_cluster_var6000' \
# 2>&1 | tee $main/'analysis/06_cluster_var6000/integrate_log.txt'


###################################################
### RUN DIMENSIONALITY REDUCTION AND CLUSTERING ###
###################################################
Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $main/'analysis/MITO/06_cluster_var6000/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'top,50' \
--var_genes 'seurat' \
--dim_reduct_use 'umap' \
--dim_reduct_params 'umap, n.neighbors=50, min.dist=0.1, spread=5, repulsion.strength=0.5, n.epochs=500, learning.rate=0.5, negative.sample.rate=7, metric="euclidean", seed.use=42;'\
'umap10, n.neighbors=50, min.dist=0.1, spread=5, repulsion.strength=0.5, n.epochs=500, learning.rate=0.5, negative.sample.rate=7, metric="euclidean", seed.use=42' \
--pre_dim_reduct 'mnn' \
--cluster_use 'all' \
--cluster_method 'louvain,HC' \
--assay 'RNA' \
--output_path $main/'analysis/MITO/06_cluster_var6000' \
2>&1 | tee $main/'analysis/MITO/06_cluster_var6000/dr_and_cluster_log.txt'






