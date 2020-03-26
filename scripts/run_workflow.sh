#! /bin/bash -l

########################
### DEFINE VARIABLES ###
########################
var_to_plot='dataset,mouse_nr,day_post_infection,organ,infection'
var_to_regress='nFeature_RNA,percent_mito,S.Score,G2M.Score'
main='/home/jonrob/projects/d_angeletti_1910'
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

At this point, we choose a clustering result yielding many small clusters
to help remove all non-B-cells from the data


########################################
### RUN CLUSTER CORRELATION ANALYSIS ###
########################################
Rscript $script_path/05_cluster_correlation.R \
--Seurat_object_path $main/'analysis/02_cluster/seurat_object.rds' \
--clustering_use 'HC_12' \
--exclude_cluster 'NONE' \
--merge_cluster '0.95,0.9,0.85,0.8,0.75,0.7' \
--output_path $main/'analysis/02_cluster/cluster_correlations' \
2>&1 | tee $main/'log/02_clust_corr.txt'



###################################
### RUN DIFFERENTIAL EXPRESSION ###
###################################
Rscript $script_path/04_diff_gene_expr.R \
--Seurat_object_path $main/'analysis/02_cluster/seurat_object.rds' \
--clustering_use 'HC_12' \
--metadata_use 'organ,infection' \
--exclude_cluster 'NONE' \
--assay 'RNA' \
--output_path $main/'analysis/03_diff_expr' \
2>&1 | tee $main/'log/03_diff_expr_log.txt'



############################
### CELL TYPE PREDICTION ###
############################
# - "main_cell_types" contains some cells, but many more types can be added (the group may have some in mind)
# - Don't need to use DE to predict which cell types are which, we can use this annotation to define them ourselves
# - Prediction is currently based on correlation, but other packages may be implemented later
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $main/'analysis/02_cluster/seurat_object.rds' \
--marker_lists $script_path/../'support_files/cell_markers/main_cell_types.csv' \
--cluster_use 'HC_12' \
--assay 'RNA' \
--output_path $main/'analysis/02_cluster/cell_type_prediction' \
2>&1 | tee $main/'log/02_cell_type_prediction_log.txt'

# The resulting differentially expressed genes associated with each cluster were examined to
# determine which clusters are unlikely to be B-cells.
# Cluster  3: NK cells (Nkg7, Gzma, Klra4), T-cells (Ctsw, Ms4a4b, Ctla2a, Cd7),   
# Cluster 12: NK cells (Nkg7, Gzma, Klra4), T-cells (Ctsw, Ms4a4b, Cd3d, Ctla2a, Cd7), 
#
# Cluster 8 is strange. Specific expression of Siglech suggests pDCs, and also Grn which 
# seems to be specific to DCs based on data from the Blood Atlas. Klk1, Ctsl, and Bst2
# (and Ctsb to a lesser extent) are quite specific to this cluster, though there is not much
# information on these genes being specific for one cell type. Also, Cd7 is highly expressed
# in cluster 8, though it is usually associated with NK and T-cells.
# Cd79a, Cd79b, and Ebf1 are generally B-cell markers, but were lowly expressed in cluster 8
# providing further evidence that this cluster is unlikely to be B-cells. 
#
# Cells belonging to these clusters (3, 8, and 12) will be removed from the data, and the pipeline
# re-run starting from the DATA INTEGRATION step.



##############################################################
### RUN DATA INTEGRATION, NORMALIZE AND GET VARIABLE GENES ###
##############################################################

# Re-integrate and normalize data with suspected non-B-cell clusters 3, 8, and 12 removed
# (specified using the 'cluster_use' input)
Rscript $script_path/02_integrate.R \
--Seurat_object_path $main/'analysis/02_cluster/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--var_genes 'seurat' \
--integration_method 'mnn,dataset' \
--cluster_use 'HC_12,1,2,4,5,6,7,9,10,11' \
--assay 'RNA' \
--output_path $main/'analysis/04_cluster' \
2>&1 | tee $main/log/'04_integrate_log.txt'


###################################################
### RUN DIMENSIONALITY REDUCTION AND CLUSTERING ###
###################################################

# Do not need to specify the B-cell only clusters for this run because we are using
# the output from the previous step which has already removed the non-B-cells
Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'var,1' \
--var_genes 'seurat' \
--dim_reduct_use 'umap' \
--cluster_use 'none' \
--cluster_method 'HC,louvain' \
--assay 'mnn' \
--output_path $main/'analysis/04_cluster' \
2>&1 | tee $main/log/'04_dr_and_cluster_log.txt'

# HC cluster #15 appears to yield a reasonable grouping of cells


########################################
### RUN CLUSTER CORRELATION ANALYSIS ###
########################################

Rscript $script_path/05_cluster_correlation.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--clustering_use 'HC_15' \
--exclude_cluster 'NONE' \
--merge_cluster '0.95,0.9,0.85,0.8,0.75,0.7' \
--output_path $main/'analysis/04_cluster/cluster_correlations' \
2>&1 | tee $main/'log/04_clust_corr.txt'



###################################
### RUN DIFFERENTIAL EXPRESSION ###
###################################
Rscript $script_path/04_diff_gene_expr.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--clustering_use 'HC_15' \
--metadata_use 'organ,infection' \
--exclude_cluster 'NONE' \
--assay 'RNA' \
--output_path $main/'analysis/05_diff_expr' \
2>&1 | tee $main/'log/05_diff_expr_log.txt'



############################
### CELL TYPE PREDICTION ###
############################
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--marker_lists $script_path/../'support_files/cell_markers/main_cell_types.csv' \
--cluster_use 'HC_15' \
--assay 'RNA' \
--output_path $main/'analysis/04_cluster/cell_type_prediction' \
2>&1 | tee $main/'log/04_cell_type_prediction_log.txt'

# The only cluster that really stood out was #14, which was associated with elevated
# expression in Fcer1g, Tyrobp, Lgals3, Alox5ap, Ifitm3, Lyz2, Ctsb, and Ccl6,
# which are generally associated with macrophages, NK cells, and/or T-cells, rather
# than B-cells.
#
# The pipeline will be re-run after removing cluster #14.



##############################################################
### RUN DATA INTEGRATION, NORMALIZE AND GET VARIABLE GENES ###
##############################################################

# Re-integrate and normalize data with suspected non-B-cell cluster #14 removed
Rscript $script_path/02_integrate.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--var_genes 'seurat' \
--integration_method 'mnn,dataset' \
--cluster_use 'HC_15,1,2,3,4,5,6,7,8,9,10,11,12,13,15' \
--assay 'RNA' \
--output_path $main/'analysis/06_cluster' \
2>&1 | tee $main/log/'06_integrate_log.txt'


###################################################
### RUN DIMENSIONALITY REDUCTION AND CLUSTERING ###
###################################################
Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $main/'analysis/06_cluster/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'var,1' \
--var_genes 'seurat' \
--dim_reduct_use 'umap' \
--cluster_use 'none' \
--cluster_method 'HC,louvain' \
--assay 'mnn' \
--output_path $main/'analysis/06_cluster' \
2>&1 | tee $main/log/'06_dr_and_cluster_log.txt'

# HC cluster #11 appears to yield a reasonable grouping of cells


########################################
### RUN CLUSTER CORRELATION ANALYSIS ###
########################################

Rscript $script_path/05_cluster_correlation.R \
--Seurat_object_path $main/'analysis/06_cluster/seurat_object.rds' \
--clustering_use 'HC_11' \
--exclude_cluster 'NONE' \
--merge_cluster '0.95,0.9,0.85,0.8,0.75,0.7' \
--output_path $main/'analysis/06_cluster/cluster_correlations' \
2>&1 | tee $main/'log/06_clust_corr.txt'


###################################
### RUN DIFFERENTIAL EXPRESSION ###
###################################
Rscript $script_path/04_diff_gene_expr.R \
--Seurat_object_path $main/'analysis/06_cluster/seurat_object.rds' \
--clustering_use 'HC_11' \
--metadata_use 'organ,infection' \
--exclude_cluster 'NONE' \
--assay 'RNA' \
--output_path $main/'analysis/07_diff_expr' \
2>&1 | tee $main/'log/07_diff_expr_log.txt'


############################
### CELL TYPE PREDICTION ###
############################
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $main/'analysis/06_cluster/seurat_object.rds' \
--marker_lists $script_path/../'support_files/cell_markers/main_cell_types.csv' \
--cluster_use 'HC_11' \
--assay 'RNA' \
--output_path $main/'analysis/06_cluster/cell_type_prediction' \
2>&1 | tee $main/'log/06_cell_type_prediction_log.txt'

# It seems that all of the non-B-cells have been removed at this point


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
--Seurat_object_path $main/'analysis/06_cluster/seurat_object.rds' \
--VDJ_annotation_path $main/'data/VDJ_OTUs' \
--columns_metadata 'infection,organ' \
--top_TCRs '10' \
--paired_only 'true' \
--only_coding_cdr3 'true' \
--same_scale 'true' \
--assay 'RNA' \
--output_path $main/'analysis/08_VDJ_analysis_paired' \
2>&1 | tee $main/log/'08_VDJ_analysis_paired_log.txt'

# run again, but unpaired
Rscript $script_path/VDJ_analysis.R \
--Seurat_object_path $main/'analysis/06_cluster/seurat_object.rds' \
--VDJ_annotation_path $main/'data/VDJ_OTUs' \
--columns_metadata 'infection,organ' \
--top_TCRs '10' \
--paired_only 'false' \
--only_coding_cdr3 'true' \
--same_scale 'true' \
--assay 'RNA' \
--output_path $main/'analysis/08_VDJ_analysis_unpaired' \
2>&1 | tee $main/log/'08_VDJ_analysis_unpaired_log.txt'




conda deactivate
