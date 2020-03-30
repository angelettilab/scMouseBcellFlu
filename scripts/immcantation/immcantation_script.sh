# Immcantation workflow

# specify main project directory
main='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'

# add immcantation scripts to path
PATH="$PATH:/Users/jonrob/Documents/NBIS/repos/immcantation/scripts"

# Download reference databases
fetch_igblastdb.sh -o $main/data/immcantation/igblast
fetch_imgtdb.sh -o $main/data/immcantation/germlines/imgt

# Build IgBLAST database from IMGT reference sequences
imgt2igblast.sh -i $main/data/immcantation/germlines/imgt -o $main/data/immcantation/igblast

# Run IgBLAST
AssignGenes.py igblast \
-s $main/data/VDJ_OTUs/lung0_2/filtered_contig_lung0_2.fasta \
-b $main/data/immcantation/igblast \
--organism mouse \
--loci ig \
--format blast

# Create tab-delimited database file to store seq alignment info
MakeDb.py igblast \
-i $main/data/VDJ_OTUs/lung0_2/filtered_contig_lung0_2_igblast.fmt7 \
-s $main/data/VDJ_OTUs/lung0_2/filtered_contig_lung0_2.fasta \
-r $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IG*.fasta \
-o $main/data/VDJ_OTUs/lung0_2/filtered_contig_lung0_2_igblast_db-pass.tab \
--10x $main/data/VDJ_OTUs/lung0_2/filtered_contig_annotations_lung0_2.csv \
--extended

# Filter database to only keep functional sequences
ParseDb.py select -d $main/data/VDJ_OTUs/lung0_2/filtered_contig_lung0_2_igblast_db-pass.tab \
-f FUNCTIONAL \
-u T \
--outname lung0_2_functional

# Parse database output into light and heavy chain change-o files
ParseDb.py select -d $main/data/VDJ_OTUs/lung0_2/lung0_2_functional_parse-select.tab \
-f LOCUS \
-u "IGH" \
--logic all \
--regex \
--outname lung0_2_heavy

ParseDb.py select -d $main/data/VDJ_OTUs/lung0_2/lung0_2_functional_parse-select.tab \
-f LOCUS \
-u "IG[LK]" \
--logic all \
--regex \
--outname lung0_2_light


# RUN R-SCRIPT TO IDENTIFY THRESHOLD FOR GROUPING B-CELL CLONES
# RESULTS EXPORTED AS heavy_lung0_2_parse-select-genotyped.tab
# AND A FASTA FILE OF THE V-SEGMENT GERMLINE SEQUENCES: lung0_2_v_genotype.fasta


# Define clones
DefineClones.py -d $main/data/VDJ_OTUs/lung0_2/lung0_2_functional_parse-select.tab \
--act set \
--model ham \
--sym min \
--norm len \
--dist 0


# Create germline sequences
# Option1: Using genotyped sequences from TIgGER
CreateGermlines.py -d $main/data/VDJ_OTUs/lung0_2/lung0_2_heavy_parse-select-genotyped.tab \
-r $main/data/VDJ_OTUs/lung0_2/lung0_2_v_genotype.fasta \
$main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGH[DJ].fasta \
-g dmask \
--vf V_CALL_GENOTYPED

# Option2: Using the clonal assignments from above
CreateGermlines.py -d $main/data/VDJ_OTUs/lung0_2/lung0_2_functional_parse-select_clone-pass.tab \
 -r $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGH[VDJ].fasta \
-g dmask \
--cloned




