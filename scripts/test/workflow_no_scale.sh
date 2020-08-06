#! /bin/bash -l

# NOTE: this was the original workflow formulation, which does not involve
# scaling or regressing the counts prior to integration. It is therefore
# recommended that this workflow is NOT used.

########################
### DEFINE VARIABLES ###
########################
var_to_plot='dataset,mouse_nr,day_post_infection,organ,infection'
var_to_regress='nFeature_RNA,percent_mito,S.Score,G2M.Score'
main='/home/jonrob/projects/d_angeletti_1910/scMouseBcellFlu'
script_path='/home/jonrob/repos/sauron/scripts'
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


#####################
### LOAD DATASETS ###
#####################

# - Rename all "features.tsv" to "genes.tsv"
# - Remove all zipped (.gz) versions of data
# - Use "metadata.csv" to include ALL data
# - Outputs object as "raw_seurat_object.rds"
Rscript $script_path/00_load_data.R \
--input_path $main/'data/cellranger' \
--dataset_metadata_path $main/'data/metadata.csv' \
--assay 'RNA' \
--output_path $main/'analysis/01_qc' \
2>&1 | tee $main/log/'00_load_data_log.txt'



###########################
### RUN QUALITY CONTROL ###
###########################

# - Outputs object as "filt_seurat_object.rds"
Rscript $script_path/01_qc_filter.R \
--Seurat_object_path $main/'analysis/01_qc/raw_seurat_object.rds' \
--columns_metadata $var_to_plot \
--species_use 'mmusculus' \
--remove_non_coding 'True' \
--plot_gene_family 'RPS,RPL,mito,HB' \
--keep_genes "$(cat $main/data/gene_lists/igh_genes_to_keep.txt)" \
--remove_gene_family 'mito' \
--min_gene_count '5' \
--min_gene_per_cell '200' \
--pct_mito_range '0,10' \
--pct_ribo_range '0,25' \
--assay 'RNA' \
--output_path $main/'analysis/01_qc' \
2>&1 | tee $main/log/'01_qc_log.txt'



##############################################################
### RUN DATA INTEGRATION, NORMALIZE AND GET VARIABLE GENES ###
##############################################################

# - Integration method will only affect visualization and clustering, but NOT DE analysis
# - Outputs "seurat_object.rds" with added integrated slot ('mnn') and variable genes
Rscript $script_path/02_integrate.R \
--Seurat_object_path $main/'analysis/01_qc/filt_seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--var_genes 'seurat' \
--integration_method 'mnn,dataset' \
--cluster_use 'NONE' \
--assay 'RNA' \
--output_path $main/'analysis/02_cluster' \
2>&1 | tee $main/log/'02_integrate_log.txt'



###################################################
### RUN DIMENSIONALITY REDUCTION AND CLUSTERING ###
###################################################

# using object integrated by dataset
Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $main/'analysis/02_cluster/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'var,1' \
--var_genes 'seurat' \
--dim_reduct_use 'umap' \
--dim_reduct_params 'umap, n.neighbors=30, min.dist=0.01, spread=3, n.epochs=500, learning.rate=0.5; umap10, n.neighbors=30, min.dist=0.01, n.epochs=500' \
--pre_dim_reduct 'mnn' \
--cluster_use 'none' \
--cluster_method 'HC,louvain' \
--assay 'RNA' \
--output_path $main/'analysis/02_cluster' \
2>&1 | tee $main/log/'03_dr_and_cluster_log.txt'

# At this point, we choose a clustering result yielding many small clusters
# to help remove all non-B-cells from the data


########################################
### RUN CLUSTER CORRELATION ANALYSIS ###
########################################
Rscript $script_path/05_cluster_correlation.R \
--Seurat_object_path $main/'analysis/02_cluster/seurat_object.rds' \
--clustering_use 'louvain_0.7' \
--exclude_cluster 'NONE' \
--merge_cluster '0.95,0.9,0.85,0.8,0.75,0.7' \
--output_path $main/'analysis/02_cluster/cluster_correlations' \
2>&1 | tee $main/'log/04_clust_corr.txt'



###################################
### RUN DIFFERENTIAL EXPRESSION ###
###################################
Rscript $script_path/04_diff_gene_expr.R \
--Seurat_object_path $main/'analysis/02_cluster/seurat_object.rds' \
--clustering_use 'louvain_0.7' \
--metadata_use 'organ,infection' \
--exclude_cluster 'NONE' \
--assay 'RNA' \
--output_path $main/'analysis/03_diff_expr' \
2>&1 | tee $main/'log/05_diff_expr_log.txt'



###########################
## CELL TYPE PREDICTION ###
###########################
# - "main_cell_types" contains some cells, but many more types can be added (the group may have some in mind)
# - Don't need to use DE to predict which cell types are which, we can use this annotation to define them ourselves
# - Prediction is currently based on correlation, but other packages may be implemented later
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $main/'analysis/02_cluster/seurat_object.rds' \
--marker_lists $main/'data/gene_lists/main_cell_types.csv' \
--assay 'RNA' \
--clustering_use 'louvain_0.7' \
--output_path $main/'analysis/02_cluster/cell_type_prediction' \
2>&1 | tee $main/'log/06_cell_type_prediction_log.txt'



#######################################
## REMOVE NON-B-CELLS FROM THE DATA ###
#######################################
# Cells that are NOT predicted as B-cells will be removed from the data, and the pipeline
# re-run starting from the DATA INTEGRATION step.
# The following code will remove cells that were in the specified clusters OR cells that
# were not predicted by the correlation analysis to be B-cells.
Rscript $main/scripts/remove_cells.R \
--Seurat_object_path $main/'analysis/02_cluster/seurat_object.rds' \
--remove 'louvain_0.7,0,11,12,13,15,16,17' \
--keep 'cell_pred_correlation_main_cell_types,B_cell' \
--combine_method 'union' \
--output_path $main/'analysis/04_cluster2' \
2>&1 | tee $main/'log/07_cell_type_prediction_log.txt'



##############################################################
### RUN DATA INTEGRATION, NORMALIZE AND GET VARIABLE GENES ###
##############################################################

# Re-integrate and normalize data now that non-B-cells have been removed.
Rscript $script_path/02_integrate.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--var_genes 'seurat' \
--integration_method 'mnn,dataset' \
--cluster_use 'NONE' \
--assay 'RNA' \
--output_path $main/'analysis/04_cluster' \
2>&1 | tee $main/log/'08_integrate_log.txt'


###################################################
### RUN DIMENSIONALITY REDUCTION AND CLUSTERING ###
###################################################

Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'var,1' \
--var_genes 'seurat' \
--dim_reduct_use 'umap' \
--dim_reduct_params 'umap, n.neighbors=30, min.dist=0.01, spread=3; umap10, n.neighbors=30' \
--pre_dim_reduct 'mnn' \
--cluster_use 'none' \
--cluster_method 'louvain' \
--assay 'RNA' \
--output_path $main/'analysis/04_cluster' \
2>&1 | tee $main/log/'09_dr_and_cluster_log.txt'

# analyze louvain_0.8 clustering


########################################
### RUN CLUSTER CORRELATION ANALYSIS ###
########################################
Rscript $script_path/05_cluster_correlation.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--clustering_use 'louvain_0.8' \
--exclude_cluster 'NONE' \
--merge_cluster '0.95,0.9,0.85,0.8,0.75,0.7' \
--output_path $main/'analysis/04_cluster/cluster_correlations' \
2>&1 | tee $main/'log/10_clust_corr.txt'



###################################
### RUN DIFFERENTIAL EXPRESSION ###
###################################
Rscript $script_path/04_diff_gene_expr.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--clustering_use 'louvain_0.8' \
--metadata_use 'organ,infection' \
--exclude_cluster 'NONE' \
--assay 'RNA' \
--output_path $main/'analysis/05_diff_expr' \
2>&1 | tee $main/'log/11_diff_expr_log.txt'



############################
### CELL TYPE PREDICTION ###
############################
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--marker_lists $main/'data/gene_lists/main_cell_types.csv' \
--clustering_use 'louvain_0.8' \
--assay 'RNA' \
--output_path $main/'analysis/04_cluster/cell_type_prediction' \
2>&1 | tee $main/'log/12_cell_type_prediction_log.txt'



################################
### SUB-CELL TYPE PREDICTION ###
################################
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--marker_lists $main/'data/gene_lists/bcell_types.csv,'$main/'data/gene_lists/bcell_types_germsub.csv,'$main/'data/gene_lists/bcell_types_germsub_zonesub.csv' \
--clustering_use 'louvain_0.8' \
--assay 'RNA' \
--output_path $main/'analysis/04_cluster/cell_type_prediction' \
2>&1 | tee $main/'log/13_cell_type_prediction_log.txt'



########################
### RUN VDJ ANALYSIS ###
########################
# Only involves Ig, as the data is from B-cells
# - top_TCRs - most abundant for visualization
# - paired_only - alpha & beta, or only alpha, beta separately (heavy and light chain in our case)
# - cdr3 is variable region within TCR or antibody - confers specificity and affinity to targeted protein
# - only_coding_cdr3 - only select options that actually code for something
# - same_scale - for different metadata comparisons
Rscript $script_path/VDJ_analysis.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--VDJ_annotation_path $main/'data/VDJ_OTUs' \
--columns_metadata 'infection,organ' \
--top_TCRs '10' \
--paired_only 'true' \
--only_coding_cdr3 'true' \
--same_scale 'true' \
--assay 'RNA' \
--output_path $main/'analysis/06_VDJ_analysis_paired' \
2>&1 | tee $main/log/'14_VDJ_analysis_paired_log.txt'

# run again, but unpaired
Rscript $script_path/VDJ_analysis.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--VDJ_annotation_path $main/'data/VDJ_OTUs' \
--columns_metadata 'infection,organ' \
--top_TCRs '10' \
--paired_only 'false' \
--only_coding_cdr3 'true' \
--same_scale 'true' \
--assay 'RNA' \
--output_path $main/'analysis/06_VDJ_analysis_unpaired' \
2>&1 | tee $main/log/'15_VDJ_analysis_unpaired_log.txt'



conda deactivate
