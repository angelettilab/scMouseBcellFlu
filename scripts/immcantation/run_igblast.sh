#!/usr/bin/env bash
# Simple IgBLAST wrapper
#
# Author:  Jason Anthony Vander Heiden
# Date:    2017.08.10
#
# Arguments:
#   -s = FASTA sequence file.
#   -o = Output directory.
#   -g = Species name. One of human or mouse. Defaults to human.
#   -t = Receptor type. One of ig or tr. Defaults to ig.
#   -b = IGDATA directory, which contains the IgBLAST database, optional_file
#        and auxillary_data directories. Defaults to /usr/local/share/igblast.
#   -i = Input type. One of nt (DNA) or aa (protein). Defaults to nt."
#   -n = Number of IgBLAST threads. Defaults to 1.
#   -h = Display help.

# Default argument values
IGDATA="/usr/local/share/igblast"
OUTDIR="."
INPUT_TYPE="nt"
NPROC=1


# Print usage
usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -s  FASTA sequence file."
    echo -e "  -o  Output directory."
    echo -e "  -g  Species name. One of human or mouse. Defaults to human."
    echo -e "  -t  Receptor type. One of ig or tr. Defaults to ig."
    echo -e "  -b  IGDATA directory, which contains the IgBLAST database,\n" \
            "     optional_file and auxillary_data directories.\n" \
            "     Defaults to /usr/local/share/igblast."
    echo -e "  -i  Input type. One of nt (DNA) or aa (protein). Defaults to nt."
    echo -e "  -n  Number of IgBLAST threads. Defaults to 1."
    echo -e "  -h  This message."
}

# Validation variables
READFILE_SET=false
OUTDIR_SET=false
SPECIES_SET=false
RECEPTOR_SET=false

# Get commandline arguments
while getopts "s:o:g:t:b:i:n:h" OPT; do
    case "$OPT" in
    s)  READFILE=$(realpath $OPTARG)
        READFILE_SET=true
        ;;
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
        ;;
    g)  SPECIES=$OPTARG
        SPECIES_SET=true
        ;;
    t)  RECEPTOR=$OPTARG
        RECEPTOR_SET=true
        ;;
    b)  IGDATA=$OPTARG
        IGDATA_SET=true
        ;;
    i)  INPUT_TYPE=$OPTARG
        ;;
    n)  NPROC=$OPTARG
        ;;
    h)  usage
        exit
        ;;
    \?) echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    :)  echo "Option -$OPTARG requires an argument" >&2
        exit 1
        ;;
    esac
done

# Exit if no input file provided
if ! $READFILE_SET; then
    echo "You must specify a FASTA file using the -s option" >&2
    exit 1
fi

# Set and check species
if ! ${SPECIES_SET}; then
    SPECIES="human"
elif [ ${SPECIES} != "human" ] && [ ${SPECIES} != "mouse" ]; then
    echo "Species (-g) must be one of human or mouse" >&2
    exit 1
fi

# Set and check receptor type
if ! ${RECEPTOR_SET}; then
    RECEPTOR="ig"
elif [ ${RECEPTOR} != "ig" ] && [ ${RECEPTOR} != "tr" ]; then
    echo "Receptor type (-t) must be one of ig or tr" >&2
    exit 1
fi

# Make output directory if it does not exist
if $OUTDIR_SET && [ ! -d "${OUTDIR}" ]; then
    mkdir -p $OUTDIR
fi

# Define IgBLAST directories and base command
export IGDATA
AUXILIARY="${SPECIES}_gl.aux"
IGBLAST_DB="${IGDATA}/database"

# Define seqtype argument. No associative assays in bash v3
case "$RECEPTOR" in
    "ig")
        SEQTYPE="Ig"
        ;;
    "tr")
        SEQTYPE="TCR"
        ;;
esac

if [ ${INPUT_TYPE} == "nt" ]; then
    IGBLAST_CMD="igblastn \
        -germline_db_V ${IGBLAST_DB}/imgt_${SPECIES}_${RECEPTOR}_v \
        -germline_db_D ${IGBLAST_DB}/imgt_${SPECIES}_${RECEPTOR}_d \
        -germline_db_J ${IGBLAST_DB}/imgt_${SPECIES}_${RECEPTOR}_j \
        -auxiliary_data ${IGDATA}/optional_file/${AUXILIARY} \
        -ig_seqtype ${SEQTYPE} -organism ${SPECIES} \
        -domain_system imgt -outfmt '7 std qseq sseq btop'"
elif [ ${INPUT_TYPE} == "aa" ]; then
    IGBLAST_CMD="igblastp \
        -germline_db_V ${IGBLAST_DB}/imgt_aa_${SPECIES}_${RECEPTOR}_v \
        -ig_seqtype ${SEQTYPE} -organism ${SPECIES} \
        -domain_system imgt -outfmt '7 std qseq sseq btop'"
fi

# Set run commmand
OUTFILE=$(basename ${READFILE})
OUTFILE="${OUTDIR}/${OUTFILE%.fasta}.fmt7"
IGBLAST_VER=$(${IGBLAST_CMD} -version | grep 'Package' |sed s/'Package: '//)
IGBLAST_RUN="${IGBLAST_CMD} -query ${READFILE} -out ${OUTFILE} -num_threads ${NPROC}"

echo $IGBLAST_RUN
# Align V(D)J segments using IgBLAST
echo -e "   START> igblast"
echo -e " VERSION> ${IGBLAST_VER}"
echo -e "  IGDATA> ${IGDATA}"
echo -e "  GERMDB> ${SPECIES}_${RECEPTOR}"
echo -e "    FILE> $(basename ${READFILE})\n"
echo -e "PROGRESS> [Running]"
eval $IGBLAST_RUN
echo -e "PROGRESS> [Done   ]\n"
echo -e "  OUTPUT> $(basename ${OUTFILE})"
echo -e "     END> igblast\n"