#!/bin/bash

# Immcantation VDJ analysis workflow

##########################
### DEFINE DIRECTORIES ###
##########################

# specify main project directory
main='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'

# add immcantation and ChangeO scripts to system path
PATH="$PATH:/Users/jonrob/Documents/NBIS/repos/immcantation/scripts"
PATH="$PATH:/Users/jonrob/Library/Python/3.7/bin"


#############################################
### RETRIEVE IG BLAST REFERENCE DATABASES ###
#############################################

# Download reference databases
fetch_igblastdb.sh -o $main/data/immcantation/igblast
fetch_imgtdb.sh -o $main/data/immcantation/germlines/imgt

# Build IgBLAST database from IMGT reference sequences
imgt2igblast.sh -i $main/data/immcantation/germlines/imgt -o $main/data/immcantation/igblast

# Get list of all samples available
sample_list=(`ls -d $main/data/VDJ_OTUs/*/ | xargs -n 1 basename`)

# Run IgBLAST
AssignGenes.py igblast \
-s $main/data/VDJ_OTUs/*/filtered_contig_*.fasta \
-b $main/data/immcantation/igblast \
--organism mouse \
--loci ig \
--format blast


for sample in ${sample_list[@]}
do
    # Create tab-delimited database file to store seq alignment info
    MakeDb.py igblast \
    -i $main'/data/VDJ_OTUs/'$sample'/filtered_contig_'$sample'_igblast.fmt7' \
    -s $main'/data/VDJ_OTUs/'$sample'/filtered_contig_'$sample'.fasta' \
    -r $main'/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IG'*'.fasta' \
    -o $main'/data/VDJ_OTUs/'$sample'/filtered_contig_'$sample'_igblast_db-pass.tab' \
    --10x $main'/data/VDJ_OTUs/'$sample'/filtered_contig_annotations_'$sample'.csv' \
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


# RUN R-SCRIPT TO IDENTIFY THRESHOLD FOR GROUPING B-CELL CLONES
# Results are exported as "IGHV-genotyped_M#.tab", where "M#" is "M1", "M2", etc. for each mouse,
# and a fasta file of the V-segment germline sequences: "IGHV_genotype_M#.fasta"


# Specify distance threshold for trimming the hierarchical clustering into B cell clones
# thresholds=


# Get list of mouse numbers to process
mouse_nums=(`ls $main/analysis/immcantation/genotyping/IGHV-genotyped_M{?,??}.tab | awk -F '[_.]' '{print $(NF-1)}'`)

for mouse in ${mouse_nums[@]}
do
    # Define clones (dist = distance threshold)
    # Output file is named "IGHV-genotyped_M#_clone-pass.tab"
    DefineClones.py -d $main'/analysis/immcantation/genotyping/IGHV-genotyped_'$mouse'.tab' \
    --act set \
    --model ham \
    --norm len \
    --dist 0.1 \
    --outname 'IGHV-genotyped_'$mouse \
    --log 'IGHV-genotyped_'$mouse'_DefineClones.log'


    # Create germline sequences using genotyped sequences from TIgGER
    # Output file is named "IGHV-genotyped_M#_germ-pass.tab"
    CreateGermlines.py -d $main'/analysis/immcantation/genotyping/IGHV-genotyped_'$mouse'_clone-pass.tab' \
    -g dmask \
    --cloned \
    -r $main'/analysis/immcantation/genotyping/IGHV_genotype_'$mouse'.fasta' \
    $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHD.fasta \
    $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHJ.fasta \
    --vf V_CALL_GENOTYPED \
    --outname 'IGHV-genotyped_'$mouse \
    --log 'IGHV-genotyped_'$mouse'_CreateGermlines.log'
done





