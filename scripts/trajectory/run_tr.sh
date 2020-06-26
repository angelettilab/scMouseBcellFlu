#!/bin/bash -l
#SBATCH -A snic2019-8-255
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 0-24:00:00
#SBATCH --mail-user=jonrob@chalmers.se
#SBATCH --mail-type=END
#SBATCH -o run_tr.log

singularity run ~/projects/d_angeletti_1910/envs/trajectory/trajectory.sif bash ~/projects/d_angeletti_1910/scripts/run_trajectory.sh


