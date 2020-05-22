#! /bin/bash -l

########################
### DEFINE VARIABLES ###
########################
# main='/home/jonrob/projects/d_angeletti_1910'
main='/cephyr/users/jonrob/Hebbe/projects/d_angeletti_1910'
# main='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'

seurat_obj_path=$main/'analysis/06_cluster_scale/seurat_object.rds'

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
mkdir $main'/analysis/'$run_name
Rscript $main/scripts/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.9,0,1,4,5,6,8,10,11,12' \
--start_cluster '0' \
--end_cluster '12' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/'$run_name \
2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


run_name='trajectory_01b'
mkdir $main'/analysis/'$run_name
Rscript $main/scripts/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.9,0,1,4,5,6,8,9,10,11,12' \
--start_cluster '0' \
--end_cluster '12' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/'$run_name \
2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


# Note: we want to try and force this trajectory through GC somehow.
run_name='trajectory_02'
mkdir $main'/analysis/'$run_name
Rscript $main/scripts/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.9,0,1,3,4,5,6,8,10,11' \
--start_cluster '0' \
--end_cluster '3' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/'$run_name \
2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


# Note: we want to try and force this trajectory through GC somehow.
run_name='trajectory_02b'
mkdir $main'/analysis/'$run_name
Rscript $main/scripts/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.9,0,1,3,4,5,6,8,9,10,11' \
--start_cluster '0' \
--end_cluster '3' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/'$run_name \
2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


run_name='trajectory_03'
mkdir $main'/analysis/'$run_name
Rscript $main/scripts/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.9,0,1,4,5,6,8,9,10,11' \
--start_cluster '0' \
--end_cluster '9' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/'$run_name \
2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


run_name='trajectory_04'
mkdir $main'/analysis/'$run_name
Rscript $main/scripts/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.9,0,1,4,5,6,8,10,11' \
--start_cluster '0' \
--end_cluster 'auto' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/'$run_name \
2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


run_name='trajectory_04b'
mkdir $main'/analysis/'$run_name
Rscript $main/scripts/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.9,0,1,4,5,6,8,9,10,11' \
--start_cluster '0' \
--end_cluster 'auto' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/'$run_name \
2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


run_name='trajectory_05'
mkdir $main'/analysis/'$run_name
Rscript $main/scripts/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.9,0,1,3,4,5,6,8,10,11,12' \
--start_cluster '0' \
--end_cluster '3,12' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/'$run_name \
2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


run_name='trajectory_05b'
mkdir $main'/analysis/'$run_name
Rscript $main/scripts/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.9,0,1,3,4,5,6,8,9,10,11,12' \
--start_cluster '0' \
--end_cluster '3,12' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/'$run_name \
2>&1 | tee $main/'analysis/'$run_name'/trajectory_inference_log.txt'


# deactivate conda environment
conda deactivate





