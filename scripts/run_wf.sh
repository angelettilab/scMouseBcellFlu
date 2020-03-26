#!/bin/bash -l
#SBATCH -A snic2019-8-255
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 0-4:00:00
#SBATCH --mail-user=jonrob@chalmers.se
#SBATCH --mail-type=END
#SBATCH -o run_wf.log

singularity run ~/projects/d_angeletti_1910/scripts/sauron/sauron.sif bash ~/projects/d_angeletti_1910/scripts/run_workflow.sh


