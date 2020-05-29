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
mkdir -p analysis/trajectory_slingshot
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
mkdir $main'/analysis/trajectory_slingshot/'$run_name
Rscript $main/scripts/trajectory/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.95,2,3,4,5,6,7,10,11,12,13' \
--start_cluster '2' \
--end_cluster '13' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/trajectory_slingshot/'$run_name \
2>&1 | tee $main/'analysis/trajectory_slingshot/'$run_name'/trajectory_inference_log.txt'


# Note: we want to try and force this trajectory through GC somehow.
run_name='trajectory_02'
mkdir $main'/analysis/trajectory_slingshot/'$run_name
Rscript $main/scripts/trajectory/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.95,1,2,3,4,5,6,7,10,11,12' \
--start_cluster '2' \
--end_cluster '1' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/trajectory_slingshot/'$run_name \
2>&1 | tee $main/'analysis/trajectory_slingshot/'$run_name'/trajectory_inference_log.txt'


# Note: we want to try and force this trajectory through GC somehow.
run_name='trajectory_02b'
mkdir $main'/analysis/trajectory_slingshot/'$run_name
Rscript $main/scripts/trajectory/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.95,1,2,3,4,5,6,7,9,10,11,12' \
--start_cluster '2' \
--end_cluster '1' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/trajectory_slingshot/'$run_name \
2>&1 | tee $main/'analysis/trajectory_slingshot/'$run_name'/trajectory_inference_log.txt'


run_name='trajectory_04'
mkdir $main'/analysis/trajectory_slingshot/'$run_name
Rscript $main/scripts/trajectory/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.95,2,3,4,5,6,7,10,11,12' \
--start_cluster '2' \
--end_cluster 'auto' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/trajectory_slingshot/'$run_name \
2>&1 | tee $main/'analysis/trajectory_slingshot/'$run_name'/trajectory_inference_log.txt'


run_name='trajectory_05'
mkdir $main'/analysis/trajectory_slingshot/'$run_name
Rscript $main/scripts/trajectory/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.95,1,2,3,4,5,6,7,10,11,12,13' \
--start_cluster '2' \
--end_cluster '1,13' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/trajectory_slingshot/'$run_name \
2>&1 | tee $main/'analysis/trajectory_slingshot/'$run_name'/trajectory_inference_log.txt'


run_name='trajectory_05b'
mkdir $main'/analysis/trajectory_slingshot/'$run_name
Rscript $main/scripts/trajectory/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'louvain_0.95,1,2,3,4,5,6,7,9,10,11,12,13' \
--start_cluster '2' \
--end_cluster '1,13' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/trajectory_slingshot/'$run_name \
2>&1 | tee $main/'analysis/trajectory_slingshot/'$run_name'/trajectory_inference_log.txt'


# deactivate conda environment
conda deactivate





