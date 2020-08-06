#! /bin/bash -l

########################
### DEFINE VARIABLES ###
########################
main='/home/jonrob/projects/d_angeletti_1910/scMouseBcellFlu'
data_path='/crex/proj/snic2019-8-255/private/jay/filezila/2019/10x_cellranger_26_apr'
annot_path=$main/'data/rna_velocity/genome_annotation_files'


cd $main

# activate conda environment
source activate velocyto-env


####################
### RUN VELOCYTO ###
####################

# get list of sample folders in data_path
sample_list=(`ls -d $data_path/*_SampleID_*/ | xargs -n 1 basename`)

# process each sample individually
for sample in ${sample_list[@]}
do
    velocyto run10x -m $annot_path/mm10_rmsk.gtf $data_path/$sample $annot_path/refdata-cellranger-mm10-3.0.0/genes/genes.gtf
    wait
done

