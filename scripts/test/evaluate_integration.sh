#! /bin/bash -l

########################
### DEFINE VARIABLES ###
########################
var_to_plot='dataset,mouse_nr,day_post_infection,infection'
var_to_regress='nFeature_RNA,percent_mito,S.Score,G2M.Score'
# main='/home/jonrob/projects/d_angeletti_1910'
# script_path='/home/jonrob/repos/sauron/scripts'
main='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'
script_path='/Users/jonrob/Documents/NBIS/repos/sauron/scripts'
analysis_path=$main'/analysis/scRNAseq_pipeline'

# make and change directories
cd $main'/analysis'
mkdir scRNAseq_pipeline
cd scRNAseq_pipeline


##################################
### ACTIVATE CONDA ENVIRONMENT ###
##################################
conda activate Sauron.v1  # macOS (local)
# source activate Sauron.v1  # linux/unix (cluster)


#######################################
### LOAD AND PROCESS SPLEEN DATASET ###
#######################################
# load data
mkdir $analysis_path'/01_qc_spleen'
Rscript $script_path/00_load_data.R \
--input_path $main/'data/cellranger' \
--dataset_metadata_path $main/'data/metadata_spleen.csv' \
--species_use 'mmusculus' \
--assay 'RNA' \
--output_path $analysis_path'/01_qc_spleen' \
2>&1 | tee $analysis_path'/01_qc_spleen/00_load_data_log.txt'

# run QC and normalize data
Rscript $script_path/01_qc_filter.R \
--Seurat_object_path $analysis_path'/01_qc_spleen/raw_seurat_object.rds' \
--columns_metadata $var_to_plot \
--species_use 'mmusculus' \
--remove_non_coding 'True' \
--plot_gene_family 'RPS,RPL,mito,HB' \
--keep_genes "$(cat $main/data/gene_lists/igh_genes_to_keep.txt)" \
--remove_gene_family 'mito' \
--min_gene_count '5' \
--min_gene_per_cell '200' \
--normalization 'LogNormalize' \
--pct_mito_range '0,10' \
--pct_ribo_range '0,25' \
--assay 'RNA' \
--output_path $analysis_path'/01_qc_spleen' \
2>&1 | tee $analysis_path'/01_qc_spleen/01_qc_log.txt'

# skip integration, but calculate HVGs
mkdir $analysis_path'/02_cluster_spleen'
Rscript $script_path/02_integrate.R \
--Seurat_object_path $analysis_path/'01_qc_spleen/filt_seurat_object.rds' \
--var_genes 'seurat' \
--integration_method 'NONE' \
--cluster_use 'all' \
--assay 'RNA' \
--output_path $analysis_path/'02_cluster_spleen' \
2>&1 | tee $analysis_path/'02_cluster_spleen/02_integrate_log.txt'

# run dimensionality reduction and clustering
Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $analysis_path/'02_cluster_spleen/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'var,1' \
--var_genes 'seurat' \
--dim_reduct_use 'umap' \
--dim_reduct_params 'umap, n.neighbors=30, min.dist=0.01, spread=3, n.epochs=500, learning.rate=0.5; umap10, n.neighbors=30, min.dist=0.01, n.epochs=500' \
--pre_dim_reduct 'pca' \
--cluster_use 'all' \
--cluster_method 'louvain' \
--assay 'RNA' \
--output_path $analysis_path/'02_cluster_spleen' \
2>&1 | tee $analysis_path/'02_cluster_spleen/03_dr_and_cluster_log.txt'

# run cell type prediction
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $analysis_path/'02_cluster_spleen/seurat_object.rds' \
--marker_lists $main/'data/gene_lists/main_cell_types.csv' \
--assay 'RNA' \
--clustering_use 'louvain_1' \
--output_path $analysis_path'/02_cluster_spleen/cell_type_prediction' \
2>&1 | tee $analysis_path'/02_cluster_spleen/cell_type_prediction/06_cell_type_prediction_log.txt'


#####################################
### LOAD AND PROCESS LUNG DATASET ###
#####################################
# load data
mkdir $analysis_path'/01_qc_lung'
Rscript $script_path/00_load_data.R \
--input_path $main/'data/cellranger' \
--dataset_metadata_path $main/'data/metadata_lung.csv' \
--species_use 'mmusculus' \
--assay 'RNA' \
--output_path $analysis_path'/01_qc_lung' \
2>&1 | tee $analysis_path'/01_qc_lung/00_load_data_log.txt'

# run QC and normalize data
Rscript $script_path/01_qc_filter.R \
--Seurat_object_path $analysis_path'/01_qc_lung/raw_seurat_object.rds' \
--columns_metadata $var_to_plot \
--species_use 'mmusculus' \
--remove_non_coding 'True' \
--plot_gene_family 'RPS,RPL,mito,HB' \
--keep_genes "$(cat $main/data/gene_lists/igh_genes_to_keep.txt)" \
--remove_gene_family 'mito' \
--min_gene_count '5' \
--min_gene_per_cell '200' \
--normalization 'LogNormalize' \
--pct_mito_range '0,10' \
--pct_ribo_range '0,25' \
--assay 'RNA' \
--output_path $analysis_path'/01_qc_lung' \
2>&1 | tee $analysis_path'/01_qc_lung/01_qc_log.txt'

# skip integration, but calculate HVGs
mkdir $analysis_path'/02_cluster_lung'
Rscript $script_path/02_integrate.R \
--Seurat_object_path $analysis_path/'01_qc_lung/filt_seurat_object.rds' \
--var_genes 'seurat' \
--integration_method 'NONE' \
--cluster_use 'all' \
--assay 'RNA' \
--output_path $analysis_path/'02_cluster_lung' \
2>&1 | tee $analysis_path/'02_cluster_lung/02_integrate_log.txt'

# run dimensionality reduction and clustering
Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $analysis_path/'02_cluster_lung/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'var,1' \
--var_genes 'seurat' \
--dim_reduct_use 'umap' \
--dim_reduct_params 'umap, n.neighbors=30, min.dist=0.01, spread=3, n.epochs=500, learning.rate=0.5; umap10, n.neighbors=30, min.dist=0.01, n.epochs=500' \
--pre_dim_reduct 'pca' \
--cluster_use 'all' \
--cluster_method 'louvain' \
--assay 'RNA' \
--output_path $analysis_path/'02_cluster_lung' \
2>&1 | tee $analysis_path/'02_cluster_lung/03_dr_and_cluster_log.txt'

# run cell type prediction
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $analysis_path/'02_cluster_lung/seurat_object.rds' \
--marker_lists $main/'data/gene_lists/main_cell_types.csv' \
--assay 'RNA' \
--clustering_use 'louvain_1' \
--output_path $analysis_path'/02_cluster_lung/cell_type_prediction' \
2>&1 | tee $analysis_path'/02_cluster_lung/cell_type_prediction/06_cell_type_prediction_log.txt'


####################################
### LOAD AND PROCESS MLN DATASET ###
####################################
# load data
mkdir $analysis_path'/01_qc_mln'
Rscript $script_path/00_load_data.R \
--input_path $main/'data/cellranger' \
--dataset_metadata_path $main/'data/metadata_mln.csv' \
--species_use 'mmusculus' \
--assay 'RNA' \
--output_path $analysis_path'/01_qc_mln' \
2>&1 | tee $analysis_path'/01_qc_mln/00_load_data_log.txt'

# run QC and normalize data
Rscript $script_path/01_qc_filter.R \
--Seurat_object_path $analysis_path'/01_qc_mln/raw_seurat_object.rds' \
--columns_metadata $var_to_plot \
--species_use 'mmusculus' \
--remove_non_coding 'True' \
--plot_gene_family 'RPS,RPL,mito,HB' \
--keep_genes "$(cat $main/data/gene_lists/igh_genes_to_keep.txt)" \
--remove_gene_family 'mito' \
--min_gene_count '5' \
--min_gene_per_cell '200' \
--normalization 'LogNormalize' \
--pct_mito_range '0,10' \
--pct_ribo_range '0,25' \
--assay 'RNA' \
--output_path $analysis_path'/01_qc_mln' \
2>&1 | tee $analysis_path'/01_qc_mln/01_qc_log.txt'

# skip integration, but calculate HVGs
mkdir $analysis_path'/02_cluster_mln'
Rscript $script_path/02_integrate.R \
--Seurat_object_path $analysis_path/'01_qc_mln/filt_seurat_object.rds' \
--var_genes 'seurat' \
--integration_method 'NONE' \
--cluster_use 'all' \
--assay 'RNA' \
--output_path $analysis_path/'02_cluster_mln' \
2>&1 | tee $analysis_path/'02_cluster_mln/02_integrate_log.txt'

# run dimensionality reduction and clustering
Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $analysis_path/'02_cluster_mln/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'var,1' \
--var_genes 'seurat' \
--dim_reduct_use 'umap' \
--dim_reduct_params 'umap, n.neighbors=30, min.dist=0.01, spread=3, n.epochs=500, learning.rate=0.5; umap10, n.neighbors=30, min.dist=0.01, n.epochs=500' \
--pre_dim_reduct 'pca' \
--cluster_use 'all' \
--cluster_method 'louvain' \
--assay 'RNA' \
--output_path $analysis_path/'02_cluster_mln' \
2>&1 | tee $analysis_path/'02_cluster_mln/03_dr_and_cluster_log.txt'

# run cell type prediction
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $analysis_path/'02_cluster_mln/seurat_object.rds' \
--marker_lists $main/'data/gene_lists/main_cell_types.csv' \
--assay 'RNA' \
--clustering_use 'louvain_1' \
--output_path $analysis_path'/02_cluster_mln/cell_type_prediction' \
2>&1 | tee $analysis_path'/02_cluster_mln/cell_type_prediction/06_cell_type_prediction_log.txt'

# run cell type prediction for B-cell subtypes
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $analysis_path/'02_cluster_mln/seurat_object.rds' \
--marker_lists $main/'data/gene_lists/bcell_types.csv' \
--assay 'RNA' \
--clustering_use 'louvain_1' \
--output_path $analysis_path'/02_cluster_mln/cell_type_prediction' \
2>&1 | tee $analysis_path'/02_cluster_mln/cell_type_prediction/06_cell_type_prediction_log.txt'




