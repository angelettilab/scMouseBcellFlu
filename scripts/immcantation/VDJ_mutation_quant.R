# Script for B-cell VDJ mutation quantification using Immcantation package

library(alakazam)
library(shazam)
library(tigger)
library(dplyr)
library(ggplot2)


# specify directory containing genotyped VDJ files
proj_dir <- '/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/'
geno_dir <- paste0(proj_dir, 'analysis/immcantation/genotyping/')
mut_dir <- paste0(proj_dir, 'analysis/immcantation/mutation/')


###############################
### LOAD GERMLINE SEQ FILES ###
###############################

# get list of available files and iterate through each file
geno_files <- dir(geno_dir, 'germ-pass[.]tab', full.names=T)
db_full <- NULL
for (g_file in geno_files) {

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
  
  
  # view first few rows of new mutation columns
  # db %>% select(SEQUENCE_ID, starts_with("MU_")) %>% head(n=4)
  
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

# write the merged database to a file
if (!dir.exists(mut_dir)) { dir.create(mut_dir) }
writeChangeoDb(db_full, paste0(mut_dir, 'VDJseq_mutation_quant.tab'))









