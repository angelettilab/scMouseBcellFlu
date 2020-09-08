#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(optparse)))

##################################
### DEFINE PATH TO LOCAL FILES ###
##################################
cat("\nRunning VDJ MUTATION QUANTIFICATION with the following parameters ...\n")
option_list = list(
  make_option(c("-d", "--changeo_db_path"), type = "character",   metavar="character",   default='none',    help="Path of the directory containing the mouse-specific folders which themselves contain Change-O seq db files."),
  make_option(c("-c", "--chain"),           type = "character",   metavar="character",   default='none',    help="BCR IG chain type, valid options are 'heavy' (IGH) or 'light' (IGK and IGL)."),
  make_option(c("-o", "--output_path"),     type = "character",   metavar="character",   default='none',    help="Output directory.")
)
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if (!dir.exists(opt$output_path)) { dir.create(opt$output_path, recursive=T) }
setwd(opt$output_path)
#---------


##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading libraries ...\n\n")
suppressMessages(suppressWarnings(library(alakazam)))
suppressMessages(suppressWarnings(library(shazam)))
suppressMessages(suppressWarnings(library(tigger)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
#---------


###############################
### LOAD GERMLINE SEQ FILES ###
###############################
# get list of available chain seq db files and iterate through each
seq_files <- file.path(list.dirs(opt$changeo_db_path, full.names=T, recursive=F), paste0('seqdb_', opt$chain, '_germ-pass.tab'))
db_full <- NULL
for (s_file in seq_files) {

  cat('Processing mouse:', basename(dirname(s_file)), '...\n')
  
  # load the Change-O database file with germline sequence information (*_germ-pass.tab file)
  db <- readChangeoDb(s_file)
  
  # calculate mutation counts and frequencies
  # mutaion counts per individual regions
  db <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT", germlineColumn="GERMLINE_IMGT_D_MASK",
                          frequency=F, regionDefinition=IMGT_V_BY_REGIONS, combine=F)
  # total mutation counts
  db <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT", germlineColumn="GERMLINE_IMGT_D_MASK",
                          frequency=F, combine=T)
  colnames(db)[colnames(db) == 'MU_COUNT'] <- 'MU_COUNT_TOT'
  
  # mutation frequency per individual regions
  db <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT", germlineColumn="GERMLINE_IMGT_D_MASK",
                          frequency=T, regionDefinition=IMGT_V_BY_REGIONS, combine=F)
  # total mutation frequency
  db <- observedMutations(db, sequenceColumn="SEQUENCE_IMGT", germlineColumn="GERMLINE_IMGT_D_MASK",
                          frequency=T, combine=T)
  colnames(db)[colnames(db) == 'MU_FREQ'] <- 'MU_FREQ_TOT'
  
  # extract some metadata from the SEQUENCE_ID column
  organ_day <- unlist(lapply(db$SEQUENCE_ID, function(x) tail(unlist(strsplit(x, '-|_')), 2)[1]))
  db$ORGAN <- gsub('\\d+', '', organ_day)
  db$DAY_POST_INFECTION <- gsub('[a-z]+', '', organ_day)  # each mouse only has one day

  # combine data with full database
  if (is.null(db_full)) {
    db_full <- db
  } else {
    db_full <- rbind(db_full, db)
  }
  
}

# write the merged ChangeO database containing mutation data to a file
out_file <- paste0('mutation_quant_', opt$chain, 'chain.tab')
cat('\nWriting merged ChangeO database file:', out_file, '...\n')
if (!dir.exists(opt$output_path)) { dir.create(opt$output_path, recursive=T) }
writeChangeoDb(db_full, file.path(opt$output_path, out_file))
cat('Done!\n\n')
#---------


##########################
### PRINT SESSION INFO ###
##########################
cat('R SESSION INFO:\n')
Sys.info()
cat('\n\n\n\n')
sessionInfo()
#---------





