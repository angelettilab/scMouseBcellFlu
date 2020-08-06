#!/bin/bash -l
#SBATCH -A snic2019-8-255
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 0-20:00:00
#SBATCH --mail-user=jonrob@chalmers.se
#SBATCH --mail-type=END
#SBATCH -o run_velo.log

singularity run ~/projects/d_angeletti_1910/scMouseBcellFlu/envs/rna_velocity/velocyto.sif bash ~/projects/d_angeletti_1910/scMouseBcellFlu/scripts/velocyto_get_spliced_unspliced.sh


