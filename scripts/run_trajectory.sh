#! /bin/bash -l

########################
### DEFINE VARIABLES ###
########################
main='/home/jonrob/projects/d_angeletti_1910'
# main='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'

cd $main
# mkdir analysis
# mkdir log



##################################
### ACTIVATE CONDA ENVIRONMENT ###
##################################
# conda activate trajectory-env  # macOS (local)
source activate trajectory-env  # linux/unix (cluster)



################################
### RUN TRAJECTORY INFERENCE ###
################################

run_name='trajectory_01'
Rscript $main/scripts/trajectory_inference.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--reduction_use 'dm' \
--reduction_visualize 'umap' \
--destiny_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.95,2,4,5,7,8,9,11,13,14,15' \
--start_cluster '13' \
--end_cluster '14' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/'$run_name \
2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


# run_name='trajectory_01b'
# Rscript $main/scripts/trajectory_inference.R \
# --Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
# --reduction_use 'dm' \
# --reduction_visualize 'umap' \
# --destiny_params 'k=30, n_eigs=20' \
# --cluster_use 'louvain_0.95,2,4,5,7,8,9,11,12,13,14,15' \
# --start_cluster '13' \
# --end_cluster '14' \
# --diff_testing 'true' \
# --assay 'RNA' \
# --output_path $main/'analysis/'$run_name \
# 2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


# # Note: we want to try and force this trajectory through GC somehow.
# run_name='trajectory_02'
# Rscript $main/scripts/trajectory_inference.R \
# --Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
# --reduction_use 'dm' \
# --reduction_visualize 'umap' \
# --destiny_params 'k=30, n_eigs=20' \
# --cluster_use 'louvain_0.95,2,3,4,5,7,8,9,11,13,15' \
# --start_cluster '13' \
# --end_cluster '3' \
# --diff_testing 'true' \
# --assay 'RNA' \
# --output_path $main/'analysis/'$run_name \
# 2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


# # Note: we want to try and force this trajectory through GC somehow.
# run_name='trajectory_02b'
# Rscript $main/scripts/trajectory_inference.R \
# --Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
# --reduction_use 'dm' \
# --reduction_visualize 'umap' \
# --destiny_params 'k=30, n_eigs=20' \
# --cluster_use 'louvain_0.95,2,3,4,5,7,8,9,11,12,13,15' \
# --start_cluster '13' \
# --end_cluster '3' \
# --diff_testing 'true' \
# --assay 'RNA' \
# --output_path $main/'analysis/'$run_name \
# 2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


run_name='trajectory_03'
Rscript $main/scripts/trajectory_inference.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--reduction_use 'dm' \
--reduction_visualize 'umap' \
--destiny_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.95,2,4,5,7,8,9,11,12,13,15' \
--start_cluster '13' \
--end_cluster '12' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/'$run_name \
2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


run_name='trajectory_04'
Rscript $main/scripts/trajectory_inference.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--reduction_use 'dm' \
--reduction_visualize 'umap' \
--destiny_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.95,2,4,5,7,8,9,11,13,15' \
--start_cluster '13' \
--end_cluster 'none' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/'$run_name \
2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


run_name='trajectory_04b'
Rscript $main/scripts/trajectory_inference.R \
--Seurat_object_path $main/'analysis/04_cluster/seurat_object.rds' \
--reduction_use 'dm' \
--reduction_visualize 'umap' \
--destiny_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.95,2,4,5,7,8,9,11,12,13,15' \
--start_cluster '13' \
--end_cluster 'none' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/'$run_name \
2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'



# deactivate conda environment
conda deactivate





