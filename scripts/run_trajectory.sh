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
# source activate trajectory-env  # linux/unix (cluster)



################################
### RUN TRAJECTORY INFERENCE ###
################################
# Rscript $main/../scripts/trajectory_inference.R \
Rscript $main/scripts/trajectory_inference.R \
--Seurat_object_path $main/'analysis/06_trajectory/seurat_object.rds' \
--reduction_use 'dm' \
--reduction_visualize 'umap' \
--cluster_use 'louvain_0.8' \
--start_cluster '1' \
--end_cluster 'none' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/06_trajectory' \
2>&1 | tee $main/log/'trajectory_log.txt'



# deactivate conda environment
# conda deactivate





