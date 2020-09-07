#!/bin/bash

# Immcantation VDJ analysis workflow

##########################
### DEFINE DIRECTORIES ###
##########################
# specify main project directory
main='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/scMouseBcellFlu'


###########################
### PREPARE ENVIRONMENT ###
###########################
source activate immcant-env

# !IMPORTANT!
# The Immcantation pipeline requires some accessory scripts that are not
# included in a nicely packaged conda or python package, but instead must be
# retrieved from their Bitbucket repository.
# All files in "https://bitbucket.org/kleinstein/immcantation/src/master/scripts/"
# should be downloaded to your local "$main/scripts/immcantation" directory.
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
# Run R-script to infer sequence genotype and estimate heavy chain hamming distance threshold.
# Results are exported into sample subfolders as "IGH_genotyped.tab" and "IGKL_genotyped.tab",
# and fasta files of the V-segment germline sequences: "IGHV_genotype.fasta" and "IGKLV_genotype.fasta".
# A .csv file of estimated threshold values is also written to the output directory.
Rscript $main/scripts/VDJ_analysis/01_VDJ_genotype_and_threshold.R \
--VDJ_data_path $main'/data/VDJ_OTUs' \
--germline_path $main'/data/immcantation/germlines' \
--density_method 'density' \
--default_threshold '0.1' \
--output_path $main'/analysis/immcantation/clone_assignment'


############################################
### DEFINE CLONES AND GERMLINE SEQUENCES ###
############################################
# If True, use V_CALL_GENOTYPED column generated from TIgGER. If False, the default V_CALL column is used.
use_genotyped=False

# extract sample names and threshold values (columns 1 and 2, respectively) from predicted_thresholds.csv
samples=(`awk -F "\"*,\"*" 'FNR > 1 {print $1}' $main/'analysis/immcantation/clone_assignment/predicted_thresholds.csv'`)
thresholds=(`awk -F "\"*,\"*" 'FNR > 1 {print $2}' $main/'analysis/immcantation/clone_assignment/predicted_thresholds.csv'`)

# create sequence 0 to #samples
indx=($(seq 0 $(( ${#samples[@]} - 1 )) ))

for i in ${indx[@]}
do
    # Define clones based on heavy chain (dist = distance threshold)
    # Output file is named "IGH_genotyped_clone-pass.tab"
    DefineClones.py -d $main'/analysis/immcantation/clone_assignment/'${samples[$i]}'/IGH_genotyped.tab' \
    --act set \
    --model ham \
    --norm len \
    --dist ${thresholds[$i]} \
    --format changeo \
    --outname 'IGH_genotyped' \
    --log $main'/analysis/immcantation/clone_assignment/'${samples[$i]}'/IGH_genotyped_DefineClones.log'

    # Incorporate light chain information to update the clone assignments
    # and remove cells assigned to more than one heavy chain.
    # Note that using the "light_parse-select.tab" or "IGKL_genotyped.tab" file yields identical results.
    light_cluster.py -d $main'/analysis/immcantation/clone_assignment/'${samples[$i]}'/IGH_genotyped_clone-pass.tab' \
    -e $main'/data/VDJ_OTUs/'${samples[$i]}'/'${samples[$i]}'_light_parse-select.tab' \
    --format changeo \
    -o $main'/analysis/immcantation/clone_assignment/'${samples[$i]}'/IGH_genotyped_clone_corrected-pass.tab'

    if $use_genotyped
    then
        # Create germline sequences for heavy chains using genotyped sequences (V_CALL_GENOTYPED) from TIgGER
        # Output file is named "IGH_genotyped_germ-pass.tab"
        CreateGermlines.py -d $main'/analysis/immcantation/clone_assignment/'${samples[$i]}'/IGH_genotyped_clone_corrected-pass.tab' \
        -g dmask \
        --cloned \
        -r $main'/analysis/immcantation/clone_assignment/'${samples[$i]}'/IGHV_genotype.fasta' \
        $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHD.fasta \
        $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHJ.fasta \
        --vf V_CALL_GENOTYPED \
        --format changeo \
        --outname 'IGH_genotyped' \
        --log $main'/analysis/immcantation/clone_assignment/'${samples[$i]}'/IGH_genotyped_CreateGermlines.log'

        # Create germline sequences for light chains using genotyped sequences (V_CALL_GENOTYPED) from TIgGER
        # Output file is named "IGKL_genotyped_germ-pass.tab"
        CreateGermlines.py -d $main'/analysis/immcantation/clone_assignment/'${samples[$i]}'/IGKL_genotyped.tab' \
        -g dmask \
        -r $main'/analysis/immcantation/clone_assignment/'${samples[$i]}'/IGKLV_genotype.fasta' \
        $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGKJ.fasta \
        $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGLJ.fasta \
        --vf V_CALL_GENOTYPED \
        --format changeo \
        --outname 'IGKL_genotyped' \
        --log $main'/analysis/immcantation/clone_assignment/'${samples[$i]}'/IGKL_genotyped_CreateGermlines.log'
    else
        # Create germline sequences for heavy chains using original V_CALL columns
        # Output file is named "IGH_genotyped_germ-pass.tab"
        CreateGermlines.py -d $main'/analysis/immcantation/clone_assignment/'${samples[$i]}'/IGH_genotyped_clone_corrected-pass.tab' \
        -g dmask \
        --cloned \
        -r $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHV.fasta \
        $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHD.fasta \
        $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHJ.fasta \
        --format changeo \
        --outname 'IGH_genotyped' \
        --log $main'/analysis/immcantation/clone_assignment/'${samples[$i]}'/IGH_genotyped_CreateGermlines.log'

        # Create germline sequences for light chains using original V_CALL columns
        # Output file is named "IGKL_genotyped_germ-pass.tab"
        CreateGermlines.py -d $main'/analysis/immcantation/clone_assignment/'${samples[$i]}'/IGKL_genotyped.tab' \
        -g dmask \
        -r $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGKV.fasta \
        $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGLV.fasta \
        $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGKJ.fasta \
        $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGLJ.fasta \
        --format changeo \
        --outname 'IGKL_genotyped' \
        --log $main'/analysis/immcantation/clone_assignment/'${samples[$i]}'/IGKL_genotyped_CreateGermlines.log'
    fi
done


####################################
### QUANTIFY VDJ MUTATION BURDEN ###
####################################
# heavy chain: exports results as ChangeO database file "IGH_mutation_quant.tab"
Rscript $main/scripts/VDJ_analysis/02_VDJ_mutation_quant.R \
--genotyped_path $main'/analysis/immcantation/clone_assignment' \
--chain 'heavy' \
--output_path $main'/analysis/immcantation/mutation'

# light chain: exports results as ChangeO database file "IGKL_mutation_quant.tab"
Rscript $main/scripts/VDJ_analysis/02_VDJ_mutation_quant.R \
--genotyped_path $main'/analysis/immcantation/clone_assignment' \
--chain 'light' \
--output_path $main'/analysis/immcantation/mutation'


# Change conda environment to Sauron.v1
conda deactivate
source activate Sauron.v1


#################################################################
### MERGE HEAVY AND LIGHT CHAIN DATA AND ADD TO SEURAT OBJECT ###
#################################################################
Rscript $main/scripts/VDJ_analysis/03_VDJ_RNAseq_integration.R \
--Seurat_object_path $main'/analysis/06_cluster/seurat_object.rds' \
--changeo_db_path $main'/analysis/immcantation/mutation' \
--output_path $main'/analysis/immcantation/seurat_object_VDJannot.rds'


conda deactivate





