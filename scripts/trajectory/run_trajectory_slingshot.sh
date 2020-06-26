#! /bin/bash -l

########################
### DEFINE VARIABLES ###
########################
main='/home/jonrob/projects/d_angeletti_1910'
# main='/cephyr/users/jonrob/Hebbe/projects/d_angeletti_1910'
# main='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'

seurat_obj_path=$main/'analysis/06_cluster/seurat_object.rds'

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
--cluster_use 'HC_16,3,5,7,8,9,10,11,12,13,14,15,16' \
--start_cluster '3' \
--end_cluster '5' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/trajectory_slingshot/'$run_name \
2>&1 | tee $main/'analysis/trajectory_slingshot/'$run_name'/trajectory_inference_log.txt'


run_name='trajectory_02'
mkdir $main'/analysis/trajectory_slingshot/'$run_name
Rscript $main/scripts/trajectory/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'HC_16,3,5,8,9,10,11,12,13,14,15,16' \
--start_cluster '3' \
--end_cluster '5' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/trajectory_slingshot/'$run_name \
2>&1 | tee $main/'analysis/trajectory_slingshot/'$run_name'/trajectory_inference_log.txt'


run_name='trajectory_03'
mkdir $main'/analysis/trajectory_slingshot/'$run_name
Rscript $main/scripts/trajectory/trajectory_slingshot_tradeSeq.R \
--Seurat_object_path $seurat_obj_path \
--pre_dim_reduct 'mnn' \
--dim_reduct_use 'dm' \
--dim_reduct_vis 'umap' \
--diffusion_params 'k=30, n_eigs=20' \
--cluster_use 'HC_16,3,4,5,7,8,9,10,11,12,13,14,15,16' \
--start_cluster '4' \
--end_cluster '5' \
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
--cluster_use 'HC_16,3,4,5,8,9,10,11,12,13,14,15,16' \
--start_cluster '4' \
--end_cluster '5' \
--diff_testing 'true' \
--assay 'RNA' \
--output_path $main/'analysis/trajectory_slingshot/'$run_name \
2>&1 | tee $main/'analysis/trajectory_slingshot/'$run_name'/trajectory_inference_log.txt'



# deactivate conda environment
conda deactivate





