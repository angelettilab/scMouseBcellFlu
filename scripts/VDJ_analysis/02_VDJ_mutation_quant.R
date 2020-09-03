#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(optparse)))

##################################
### DEFINE PATH TO LOCAL FILES ###
##################################
cat("\nRunning VDJ MUTATION QUANTIFICATION with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--genotyped_path"),      type = "character",   metavar="character",   default='none',    help="Path of the directory containing the output genotype (*_germ-pass.tab) files generated using the ChangeO CreateGermlines function."),
  make_option(c("-o", "--output_path"),         type = "character",   metavar="character",   default='none',    help="Output directory.")
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

# get list of available germline files and iterate through each
geno_files <- file.path(list.dirs(opt$genotyped_path, full.names=T, recursive=F), 'IGHLK-genotyped_germ-pass.tab')
db_full <- NULL
for (g_file in geno_files) {

  cat('Processing sample:', basename(dirname(g_file)), '...\n')
  
  # load the Change-O database file with germline sequence information (*_germ-pass.tab file)
  db <- readChangeoDb(g_file)
  
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
  
  # plot differences in mutation frequency between organs
  # ggplot(db, aes(x=ORGAN, y=MU_FREQ, fill=ORGAN)) +
  #   theme_bw() + ggtitle("Mutation quantification") +
  #   xlab("Organ") + ylab("Mutation frequency") +
  #   geom_boxplot()

  # combine data with full database
  if (is.null(db_full)) {
    db_full <- db
  } else {
    db_full <- rbind(db_full, db)
  }
  
}

# write the merged ChangeO database containing mutation data to a file
out_file <- 'VDJseq_mutation_quant.tab'
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





