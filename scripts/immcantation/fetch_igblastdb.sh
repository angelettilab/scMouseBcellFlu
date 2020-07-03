#!/usr/bin/env bash
# Download IgBLAST database files
#
# Author:  Jason Vander Heiden
# Date:    2016.11.21
#
# Arguments:
#   -o = Output directory for downloaded files. Defaults to current directory.
#   -h = Display help.

# Default argument values
OUTDIR="."

# Print usage
usage () {
    echo "Usage: `basename $0` [OPTIONS]"
    echo "  -o  Output directory for downloaded files. Defaults to current directory."
    echo "  -h  This message."
}

# Get commandline arguments
while getopts "o:h" OPT; do
    case "$OPT" in
    o)  OUTDIR=$OPTARG
        OUTDIR_SET=true
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

# Make output directory if it does not exist
if $OUTDIR_SET && [ ! -d "${OUTDIR}" ]; then
    mkdir -p $OUTDIR
fi

# Fetch internal_data
wget -q -r -nH --cut-dirs=5 --no-parent \
    ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data \
    -P ${OUTDIR}/internal_data

# Fetch database
wget -q -r -nH --cut-dirs=5 --no-parent \
    ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/database \
    -P ${OUTDIR}/database

# Fetch optional_file
wget -q -r -nH --cut-dirs=5 --no-parent \
    ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/optional_file  \
    -P ${OUTDIR}/optional_file