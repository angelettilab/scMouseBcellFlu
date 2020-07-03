#!/usr/bin/env Rscript

# Writes a TSV containing column used by VDJTools.
# Takes as input a list of Change-O db files as positional arguments.
# Conversion is approximate, as source alignment information may differ
# between MiGMAP and vanilla IgBLAST.
# Does not perform collapse of records into unique junctions.
# Note, VDJTools uses the annotation CDR3 for the junction region,
# so cdr3nt/cdr3aa contain the conserved residues flanking the CDR3.

# Imports
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(alakazam))

# Arguments
args <- commandArgs(trailingOnly=TRUE)

# Loop through files and write output
for (f in args) {
    # Load Change-O db
    db <- readChangeoDb(f)

    # Define VDJTools table
    vdj <- db %>%
        transmute(count=DUPCOUNT,
                  freq=DUPCOUNT/sum(DUPCOUNT, na.rm=TRUE),
                  cdr3nt=JUNCTION,
                  cdr3aa=translateDNA(JUNCTION),
                  v=getGene(V_CALL),
                  d=getGene(D_CALL),
                  j=getGene(J_CALL),
                  VEnd=(V_GERM_LENGTH_IMGT - 312 + 3),
                  DStart=VEnd + NP1_LENGTH,
                  DEnd=DStart + D_GERM_LENGTH,
                  JStart=DEnd + NP2_LENGTH)

    # Write output file
    n <- file_path_sans_ext(f)
    output_file <- sub(paste0("^", n), paste0(n, "_vdjtools"), f)
    cat("WRITING> ", output_file, "\n", sep="")
    write.table(vdj, file=output_file, quote=FALSE, sep="\t", row.names=FALSE)
}

