#!/bin/bash

# Immcantation VDJ analysis workflow

##########################
### DEFINE DIRECTORIES ###
##########################
# specify main project directory
main='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'


##################################
### ACTIVATE CONDA ENVIRONMENT ###
##################################
source activate immcant-env
PATH="$PATH:"$main"/scripts/immcantation"  # add immcantation scripts directory to the PATH


#############################################
### RETRIEVE IG BLAST REFERENCE DATABASES ###
#############################################
# Download reference databases
fetch_igblastdb.sh -o $main/data/immcantation/igblast
fetch_imgtdb.sh -o $main/data/immcantation/germlines/imgt

# Build IgBLAST database from IMGT reference sequences
# NOTE! This script uses the "realpath" command, which can be installed on MacOS with `brew install coreutils`
imgt2igblast.sh -i $main/data/immcantation/germlines/imgt -o $main/data/immcantation/igblast


#######################################
### RUN IG BLAST ON VDJ FASTA FILES ###
#######################################
AssignGenes.py igblast \
-s $main/data/VDJ_OTUs/*/filtered_contig_*.fasta \
-b $main/data/immcantation/igblast \
--organism mouse \
--loci ig \
--format blast


##############################################
### PARSE VDJ FILES INTO CHANGEO DB FORMAT ###
##############################################
# Get list of all samples available
sample_list=(`ls -d $main/data/VDJ_OTUs/*/ | xargs -n 1 basename`)

# Create filtered VDJ seq database files for each sample
for sample in ${sample_list[@]}
do
    # Create tab-delimited database file to store seq alignment info
    MakeDb.py igblast \
    -i $main'/data/VDJ_OTUs/'$sample'/filtered_contig_'$sample'_igblast.fmt7' \
    -s $main'/data/VDJ_OTUs/'$sample'/filtered_contig_'$sample'.fasta' \
    -r $main'/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IG'*'.fasta' \
    -o $main'/data/VDJ_OTUs/'$sample'/filtered_contig_'$sample'_igblast_db-pass.tab' \
    --10x $main'/data/VDJ_OTUs/'$sample'/filtered_contig_annotations_'$sample'.csv' \
    --format changeo \
    --extended

    # Filter database to only keep functional sequences
    ParseDb.py select -d $main'/data/VDJ_OTUs/'$sample'/filtered_contig_'$sample'_igblast_db-pass.tab' \
    -f FUNCTIONAL \
    -u T \
    --outname $sample'_functional'

    # Parse database output into light and heavy chain change-o files
    ParseDb.py select -d $main'/data/VDJ_OTUs/'$sample'/'$sample'_functional_parse-select.tab' \
    -f LOCUS \
    -u "IGH" \
    --logic all \
    --regex \
    --outname $sample'_heavy'

    ParseDb.py select -d $main'/data/VDJ_OTUs/'$sample'/'$sample'_functional_parse-select.tab' \
    -f LOCUS \
    -u "IG[LK]" \
    --logic all \
    --regex \
    --outname $sample'_light'
done


#########################################################
### INFER IG GENOTYPE AND ESTIMATE SEQ DIST THRESHOLD ###
#########################################################
# Run R-script to infer sequence genotype and estimate VDJ seq hamming distance threshold
# Results are exported as "IGHV-genotyped_M#.tab", where "M#" is "M1", "M2", etc. for each mouse,
# and a fasta file of the V-segment germline sequences: "IGHV_genotype_M#.fasta".
# A .csv file of estimated threshold values will also be written to the output directory.
Rscript $main/scripts/VDJ_analysis/01_VDJ_genotype_and_threshold.R \
--VDJ_data_path $main'/data/VDJ_OTUs' \
--metadata_file $main'/data/metadata.csv' \
--germline_path $main'/data/immcantation/germlines' \
--density_method 'density' \
--default_threshold '0.1' \
--output_path $main'/analysis/immcantation'


############################################
### DEFINE CLONES AND GERMLINE SEQUENCES ###
############################################
# extract mouse numbers and threshold values (columns 1 and 2, respectively) from predicted_thresholds.csv
mouse_nums=(`awk -F "\"*,\"*" 'FNR > 1 {print $1}' $main/'analysis/immcantation/threshold_estimation/predicted_thresholds.csv'`)
thresholds=(`awk -F "\"*,\"*" 'FNR > 1 {print $2}' $main/'analysis/immcantation/threshold_estimation/predicted_thresholds.csv'`)

# create sequence 0 to #mice
indx=($(seq 0 $(( ${#mouse_nums[@]} - 1 )) ))

# for mouse in ${mouse_nums[@]}
for i in ${indx[@]}
do
    # Define clones (dist = distance threshold)
    # Output file is named "IGHV-genotyped_M#_clone-pass.tab"
    DefineClones.py -d $main'/analysis/immcantation/genotyping/IGHV-genotyped_'${mouse_nums[$i]}'.tab' \
    --act set \
    --model ham \
    --norm len \
    --dist ${thresholds[$i]} \
    --format changeo \
    --outname 'IGHV-genotyped_'${mouse_nums[$i]} \
    --log $main'/analysis/immcantation/genotyping/IGHV-genotyped_'${mouse_nums[$i]}'_DefineClones.log'

    # Create germline sequences using genotyped sequences from TIgGER
    # Output file is named "IGHV-genotyped_M#_germ-pass.tab"
    CreateGermlines.py -d $main'/analysis/immcantation/genotyping/IGHV-genotyped_'${mouse_nums[$i]}'_clone-pass.tab' \
    -g dmask \
    --cloned \
    -r $main'/analysis/immcantation/genotyping/IGHV_genotype_'${mouse_nums[$i]}'.fasta' \
    $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHD.fasta \
    $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHJ.fasta \
    --vf V_CALL_GENOTYPED \
    --format changeo \
    --outname 'IGHV-genotyped_'${mouse_nums[$i]} \
    --log $main'/analysis/immcantation/genotyping/IGHV-genotyped_'${mouse_nums[$i]}'_CreateGermlines.log'
done



####################################
### QUANTIFY VDJ MUTATION BURDEN ###
####################################
# exports results as ChangeO database file "VDJseq_mutation_quant.tab"
Rscript $main/scripts/VDJ_analysis/02_VDJ_mutation_quant.R \
--genotyped_path $main'/analysis/immcantation/genotyping' \
--output_path $main'/analysis/immcantation/mutation'


# Change conda environment to Sauron.v1
conda deactivate
source activate Sauron.v1


#####################################
### ADD VDJ DATA TO SEURAT OBJECT ###
#####################################
Rscript $main/scripts/VDJ_analysis/03_VDJ_RNAseq_integration.R \
--Seurat_object_path $main'/analysis/06_cluster/seurat_object.rds' \
--changeo_db_path $main'/analysis/immcantation/mutation/VDJseq_mutation_quant.tab' \
--output_path $main'/analysis/immcantation/seurat_object_VDJannot.rds'


conda deactivate





