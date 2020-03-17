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
--outname functional_lung0_2

# Parse database output into light and heavy chain change-o files
ParseDb.py select -d $main/data/VDJ_OTUs/lung0_2/functional_lung0_2_parse-select.tab \
-f LOCUS \
-u "IGH" \
--logic all \
--regex \
--outname heavy_lung0_2

ParseDb.py select -d $main/data/VDJ_OTUs/lung0_2/functional_lung0_2_parse-select.tab \
-f LOCUS \
-u "IG[LK]" \
--logic all \
--regex \
--outname light_lung0_2


