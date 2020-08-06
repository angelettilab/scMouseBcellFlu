#! /bin/bash -l

########################
### DEFINE VARIABLES ###
########################
var_to_plot='dataset,mouse_nr,day_post_infection,organ,infection'
var_to_regress='nCount_RNA,nFeature_RNA,perc_mito,S.Score,G2M.Score'
# main='/home/jonrob/projects/d_angeletti_1910'
# script_path='/home/jonrob/repos/sauron/scripts'
main='/cephyr/users/jonrob/Hebbe/projects/d_angeletti_1910/scMouseBcellFlu'
script_path='/cephyr/users/jonrob/Hebbe/repos/sauron/scripts'
# main='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'
# script_path='/Users/jonrob/Documents/NBIS/repos/sauron/scripts'

cd $main
mkdir analysis


##################################
### ACTIVATE CONDA ENVIRONMENT ###
##################################
# conda activate sauron-env  # macOS (local)
source activate Sauron.v1  # linux/unix (cluster)


# Start from Seurat object in 05_qc_scale from worflow_scale_integrate.sh pipeline


###################################################
### RUN INTEGRATION WITH PRE-SCALING/REGRESSION ###
###################################################
mkdir $main/'analysis/06_cluster_scaleCC'
Rscript $main/scripts/02b_scale_integrate.R \
--Seurat_object_path $main/'analysis/05_qc_scale/filt_seurat_object.rds' \
--sauron_script_path $script_path \
--regress $var_to_regress \
--integration_method 'mnn,dataset' \
--var_genes 'seurat' \
--cluster_use 'all' \
--assay 'RNA' \
--output_path $main/'analysis/06_cluster_scaleCC' \
2>&1 | tee $main/'analysis/06_cluster_scaleCC/integrate_log.txt'



###################################################
### RUN DIMENSIONALITY REDUCTION AND CLUSTERING ###
###################################################
Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $main/'analysis/06_cluster_scaleCC/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'top,50' \
--var_genes 'seurat' \
--dim_reduct_use 'umap' \
--dim_reduct_params 'umap, n.neighbors=50, min.dist=0.1, spread=4, repulsion.strength=0.5, n.epochs=500, learning.rate=0.5, negative.sample.rate=7, metric="euclidean", seed.use=42;'\
'umap10, n.neighbors=50, min.dist=0.1, spread=4, repulsion.strength=0.5, n.epochs=500, learning.rate=0.5, negative.sample.rate=7, metric="euclidean", seed.use=42' \
--pre_dim_reduct 'mnn' \
--cluster_use 'all' \
--cluster_method 'louvain,HC,hdbscan' \
--assay 'RNA' \
--output_path $main/'analysis/06_cluster_scaleCC' \
2>&1 | tee $main/'analysis/06_cluster_scaleCC/dr_and_cluster_log.txt'



############################
### CELL TYPE PREDICTION ###
############################
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $main/'analysis/06_cluster_scaleCC/seurat_object.rds' \
--marker_lists $main/'data/gene_lists/main_cell_types.csv,'$main/'data/gene_lists/bcell_types.csv,'$main/'data/gene_lists/bcell_types_germsub.csv,'$main/'data/gene_lists/bcell_types_germsub_zonesub.csv' \
--clustering_use 'louvain_0.9' \
--assay 'RNA' \
--output_path $main/'analysis/06_cluster_scaleCC/cell_type_prediction' \
2>&1 | tee $main/'analysis/06_cluster_scaleCC/cell_subtype_prediction_log.txt'



###################################
### RUN DIFFERENTIAL EXPRESSION ###
###################################
Rscript $script_path/04_diff_gene_expr.R \
--Seurat_object_path $main/'analysis/06_cluster_scaleCC/seurat_object.rds' \
--clustering_use 'louvain_0.9' \
--metadata_use 'organ,infection' \
--exclude_cluster 'NONE' \
--assay 'RNA' \
--output_path $main/'analysis/07_diff_expr_scaleCC' \
2>&1 | tee $main/'analysis/07_diff_expr_scaleCC/diff_expr_log.txt'



# ############################
# ### VARY UMAP PARAMETERS ###
# ############################
# Rscript $main/scripts/plotting/vary_umap_params.R \
# --Seurat_object_path $main/'analysis/06_cluster_scale/seurat_object.rds' \
# --plot_groups 'none' \
# --plot_features 'NaiveBcell,GC_DarkZone,GC_LightZone' \
# --pre_dim_reduct 'mnn' \
# --umap_base_params 'n.neighbors=50, min.dist=0.1, spread=4, repulsion.strength=0.5, n.epochs=500, learning.rate=0.5, negative.sample.rate=7, metric="euclidean", seed.use=42' \
# --umap_vary_params 'n.neighbors,10,20,50,100,200; metric,"euclidean","correlation","cosine"; dims,10,20,30,50; n.epochs,200,300,400,500; learning.rate,0.1,0.3,0.5,1,2; min.dist,0.001,0.01,0.1,0.5; spread,0.1,0.5,1,3,5; set.op.mix.ratio,0,0.25,0.5,0.75,1; repulsion.strength,0.1,0.2,0.5,1,2,5; negative.sample.rate,1,2,5,7,10' \
# --assay 'RNA' \
# --add_metadata $main/'analysis/06_cluster_scale/cell_type_prediction/bcell_types_germsub/cell_pred_correlation_bcell_types_germsub.csv' \
# --output_path $main/'analysis/06_cluster_scale/umap_param_variation' \
# 2>&1 | tee $main/'analysis/06_cluster_scale/vary_umap_params_log.txt'




# ########################
# ### RUN VDJ ANALYSIS ###
# ########################
# # Only involves Ig, as the data is from B-cells
# # - top_TCRs - most abundant for visualization
# # - paired_only - alpha & beta, or only alpha, beta separately (heavy and light chain in our case)
# # - cdr3 is variable region within TCR or antibody - confers specificity and affinity to targeted protein
# # - only_coding_cdr3 - only select options that actually code for something
# # - same_scale - for different metadata comparisons
# Rscript $script_path/VDJ_analysis.R \
# --Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
# --VDJ_annotation_path $main/'data/VDJ_OTUs' \
# --columns_metadata 'infection,organ' \
# --top_TCRs '10' \
# --paired_only 'true' \
# --only_coding_cdr3 'true' \
# --same_scale 'true' \
# --assay 'RNA' \
# --output_path $main/'analysis/06_VDJ_analysis_paired' \
# 2>&1 | tee $main/log/'14_VDJ_analysis_paired_log.txt'

# # run again, but unpaired
# Rscript $script_path/VDJ_analysis.R \
# --Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
# --VDJ_annotation_path $main/'data/VDJ_OTUs' \
# --columns_metadata 'infection,organ' \
# --top_TCRs '10' \
# --paired_only 'false' \
# --only_coding_cdr3 'true' \
# --same_scale 'true' \
# --assay 'RNA' \
# --output_path $main/'analysis/06_VDJ_analysis_unpaired' \
# 2>&1 | tee $main/log/'15_VDJ_analysis_unpaired_log.txt'



conda deactivate
