#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(optparse)))

##################################
### DEFINE PATH TO LOCAL FILES ###
##################################
cat("\nRunning NEAREST NEIGHBOR ANALYSIS with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--VDJ_data_path"),       type = "character",   metavar="character",   default='none',    help="Path of the directory containing VDJ fasta files."),
  make_option(c("-m", "--metadata_file"),       type = "character",   metavar="character",   default='none',    help="Filename (including path) of sample metadata.csv file."),
  make_option(c("-d", "--density_method"),      type = "character",   metavar="character",   default='density', help="Method for finding the VDJ sequence hamming distance threshold defining distinct clones. Options are 'density', 'gmm', or 'none'. For 'gmm', one can also specify the model after a comma; e.g., 'gmm,gamma-gamma' (default), 'gmm,gamma-norm', 'gmm,norm-norm'."),
  make_option(c("-t", "--default_threshold"),   type = "character",   metavar="character",   default='0.1',     help="Default nearest neighbor hamming distance threshold to use if the threshold estimation fails or is not performed."),
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
cat("\nLoading libraries ...\n")
suppressMessages(suppressWarnings(library(alakazam)))
suppressMessages(suppressWarnings(library(shazam)))
suppressMessages(suppressWarnings(library(tigger)))
#---------


################################################
### PREPARE DIRECTORIES AND COLLECT METADATA ###
################################################
# get list of all available samples
samples <- list.dirs(opt$VDJ_data_path, full.names=F, recursive=F)

# get sample metadata
sample.meta <- read.csv(opt$metadata_file, stringsAsFactors=F)
sample.meta <- sample.meta[match(samples, sample.meta$dataset), ]  # align to sample list
rownames(sample.meta) <- seq(nrow(sample.meta))

# create directory for results output
invisible( if (!dir.exists(opt$output_path)) {dir.create(opt$output_path, recursive=T)} )
#---------


##########################################
### LOAD AND PREPARE SEQUENCE DB FILES ###
##########################################
# import ChangeO-formatted sequence database files
seqdb_heavy <- seqdb_light <- NULL
for (i in seq(length(samples))){
  # The sample name is appended to the SEQUENCE_ID and CELL identifiers to ensure
  # that there is no accidental overlap between samples.
  
  # heavy chain (IGH)
  tmp_heavy <- readChangeoDb(file.path(opt$VDJ_data_path, samples[i], paste0(samples[i], '_heavy_parse-select.tab')))
  tmp_heavy$SEQUENCE_ID <- paste0(tmp_heavy$SEQUENCE_ID, '-', samples[i])
  tmp_heavy$CELL <- paste(unlist(lapply(tmp_heavy$CELL, function(x) unlist(strsplit(x, '-'))[-2])), samples[i], sep='_')
  tmp_heavy$MOUSE_NR <- sample.meta$mouse_nr[i]
  
  # remove cells from DB if they contain multiple heavy chains
  dup_cells <- tmp_heavy$CELL[duplicated(tmp_heavy$CELL)]
  tmp_heavy <- tmp_heavy[!(tmp_heavy$CELL %in% dup_cells), ]
  seqdb_heavy <- rbind(seqdb_heavy, tmp_heavy)
  
  # light chain (IGK and IGL)
  tmp_light <- readChangeoDb(file.path(opt$VDJ_data_path, samples[i], paste0(samples[i], '_light_parse-select.tab')))
  tmp_light$SEQUENCE_ID <- paste0(tmp_light$SEQUENCE_ID, '-', samples[i])
  tmp_light$CELL <- paste(unlist(lapply(tmp_light$CELL, function(x) unlist(strsplit(x, '-'))[-2])), samples[i], sep='_')
  tmp_light$MOUSE_NR <- sample.meta$mouse_nr[i]
  seqdb_light <- rbind(seqdb_light, tmp_light)
}
#---------


###################################
### CALCULATE NEAREST NEIGHBORS ###
###################################
mice <- unique(sample.meta$mouse_nr)
mice <- mice[order(as.numeric(gsub('M','',mice)))]
predicted_thresholds <- data.frame(mouse=mice, threshold=rep(as.numeric(opt$default_threshold), length(mice)))
for (m in mice) {
  
  cat('\nPROCESSING MOUSE:', m, '\n')
  mouse_path <- file.path(opt$output_path, m)
  invisible( if (!dir.exists(mouse_path)) {dir.create(mouse_path, recursive=T)} )
  seqdb_heavy_mouse <- seqdb_heavy[seqdb_heavy$MOUSE_NR %in% m, ]
  seqdb_light_mouse <- seqdb_light[seqdb_light$MOUSE_NR %in% m, ]
  
  # export combined sequence DB for mouse
  writeChangeoDb(seqdb_heavy_mouse, file.path(mouse_path, 'seqdb_heavy.tab'))
  writeChangeoDb(seqdb_light_mouse, file.path(mouse_path, 'seqdb_light.tab'))

  # calculate distances to nearest neighbors
  seqdb_mouse <- rbind(seqdb_heavy_mouse, seqdb_light_mouse)
  dist_ham <- distToNearest(seqdb_mouse, vCallColumn='V_CALL', model='ham', normalize='len', nproc=1,
                            cellIdColumn='CELL', locusColumn='LOCUS', groupUsingOnlyIGH=F)

  # plot distance distribution
  png(filename=file.path(mouse_path, 'distToNearestNeighbor.png'), units='mm', height=150, width=180, res=300)
  print(ggplot(subset(dist_ham, !is.na(DIST_NEAREST)), aes(x=DIST_NEAREST)) +
          theme_bw() +
          xlab("Hamming distance") +
          ylab("Count") +
          scale_x_continuous(breaks=seq(0, 1, 0.1)) +
          geom_histogram(color="white", binwidth=0.02))
  invisible(dev.off())
  
  # perform automatic threshold estimation
  dist_args <- trimws(unlist(strsplit(casefold(opt$density_method), ',')))
  if (dist_args[1] == 'none') {
    
    cat('\tSkipping automatic threshold estimation. Threshold value set to default:', opt$default_threshold, '\n')
    
  } else if (dist_args[1] %in% c('density','gmm')) {
    
    # Find threshold using the density or gmm (mixture model) methods
    if (dist_args[1] == 'density') {
      output <- findThreshold(dist_ham$DIST_NEAREST, method=dist_args[1])
    } else {
      if (is.na(dist_args[2])) { dist_args[2] <- 'gamma-gamma' }
      output <- findThreshold(dist_ham$DIST_NEAREST, method=dist_args[1], model=dist_args[2])  
    }
    
    if (is.null(output) || is.na(output@threshold)) {  
      cat('\tThreshold estimation failed. Reverting to default threshold:', opt$default_threshold, '\n')
    } else {
      cat('\tEstimated hamming distance threshold:', round(output@threshold, 5), '\n')
      predicted_thresholds$threshold[predicted_thresholds$mouse == m] <- output@threshold
      png(file.path(mouse_path, paste0('distToNearestNeighbor_', paste(dist_args, collapse='_'), '_fit.png')), units='mm', height=120, width=150, res=300)
      plot(output, binwidth=0.02, title=paste0('Threshold Prediction [threshold = ', round(output@threshold, 3), ']'))
      invisible(dev.off())
    }
    
  } else {
    stop(paste0('Invalid density_method: "', dist_args[1], '". Valid options are "density", "gmm", or "none".'))
  }

}

# export table of predicted thresholds
write.csv(predicted_thresholds, file=file.path(opt$output_path, 'predicted_thresholds.csv'), quote=F, row.names=F)
#---------


##########################
### PRINT SESSION INFO ###
##########################
cat('R SESSION INFO:\n')
Sys.info()
cat('\n\n\n\n')
sessionInfo()
#---------





