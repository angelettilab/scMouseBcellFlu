#! /bin/bash -l

########################
### DEFINE VARIABLES ###
########################
main='/home/jonrob/projects/d_angeletti_1910'
data_path='/crex/proj/snic2019-8-255/private/jay/filezila'
# main='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'
# data_path=$main/'data/cellranger_full'
annot_path=$main/'data/rna_velocity/genome_annotation_files'


cd $main

# activate conda environment
source activate velocyto-env


####################
### RUN VELOCYTO ###
####################
velocyto run10x -m $annot_path/mm10_rmsk.gtf $data_path/G20-001_SampleID_9_24feb20 $annot_path/refdata-cellranger-mm10-3.0.0/genes/genes.gtf




