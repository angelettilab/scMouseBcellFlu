#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(optparse)))

##################################
### DEFINE PATH TO LOCAL FILES ###
##################################
cat("\nRunning GENOTYPE INFERENCE and NEAREST NEIGHBOR ANALYSIS with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--VDJ_data_path"),       type = "character",   metavar="character",   default='none',    help="Path of the directory containing VDJ fasta files."),
  make_option(c("-m", "--metadata_file"),       type = "character",   metavar="character",   default='none',    help="Filename (including path) of sample metadata.csv file."),
  make_option(c("-g", "--germline_path"),       type = "character",   metavar="character",   default='none',    help="Path of the directory containing IMGT germline data (should contain the 'imgt' folder)"),
  make_option(c("-d", "--density_method"),      type = "character",   metavar="character",   default='density', help="Method for finding the VDJ sequence hamming distance threshold defining distinct clones. Options are 'density', 'gmm', or 'none'. For 'gmm', one can also specify the model after a comma; e.g., 'gmm,gamma-gamma' (default), 'gmm,gamma-norm', 'gmm,norm-norm'."),
  make_option(c("-t", "--default_threshold"),   type = "character",   metavar="character",   default='0.1',     help="Default nearest neighbor hamming distance threshold to use if the threshold estimation fails or is not performed."),
  make_option(c("-o", "--output_path"),         type = "character",   metavar="character",   default='none',    help="Output directory.")
)
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

# TODO: if some path(s) are not provided (i.e., option = 'none'), try to construct a default using relative paths.

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

# create directories for results output
genotype_dir <- paste0(opt$output_path,'/genotyping')
threshold_dir <- paste0(opt$output_path,'/threshold_estimation')
invisible( if (!dir.exists(genotype_dir)) {dir.create(genotype_dir, recursive=T)} )
invisible( if (!dir.exists(threshold_dir)) {dir.create(threshold_dir, recursive=T)} )
#---------


#####################################################
### FIND NOVEL VDJ SEQ ALLELES AND INFER GENOTYPE ###
#####################################################

# load V-segment germline sequences
ighv <- readIgFasta(paste0(opt$germline_path, '/imgt/mouse/vdj/imgt_mouse_IGHV.fasta'))

# import ChangeO-formatted sequence database files (heavy chain)
seqdb <- NULL
for (i in seq(length(samples))){
  tmp <- readChangeoDb(file.path(opt$VDJ_data_path, samples[i], paste0(samples[i], '_heavy_parse-select.tab')))
  tmp$SEQUENCE_ID <- paste0(tmp$SEQUENCE_ID, '-', samples[i])
  tmp$MOUSE_NR <- sample.meta$mouse_nr[i]
  seqdb <- rbind(seqdb, tmp)
}

# infer genotype (performed separately for each mouse)
mice <- unique(seqdb$MOUSE_NR)
mice <- mice[order(as.numeric(gsub('M','',mice)))]
predicted_thresholds <- data.frame(mouse=mice, threshold=rep(as.numeric(opt$default_threshold), length(mice)))
for (m in mice) {
  
  cat('\n\nPROCESSING MOUSE', m, '\n\n')
  seqdb_mouse <- seqdb[seqdb$MOUSE_NR %in% m, ]
  
  # find novel alleles (if any)
  cat('\tSearching for novel alleles ...\n')
  novel_rows <- NULL
  try (nv <- findNovelAlleles(seqdb_mouse, ighv))
  try (novel_rows <- selectNovel(nv))
  
  # Extract and view the rows that contain successful novel allele calls
  if (!is.null(novel_rows) && (nrow(novel_rows) > 0)) {
    png(filename=paste0(genotype_dir, '/novel_alleles_', m, '.png'), units='mm', height=250, width=180, res=300)
    plotNovel(seqdb_mouse, novel_rows[1, ])  # only plot first novel allele
    invisible(dev.off())
  }
  
  # infer mouse genotype
  gt_mouse <- inferGenotype(seqdb_mouse, germline_db=ighv, novel=nv)
  png(filename=paste0(genotype_dir, '/IGHV_genotype_plot_', m, '.png'), units='mm', height=250, width=100, res=300)
  plotGenotype(gt_mouse, gene_sort="position", text_size=8) #, facet_by='ALLELES')
  invisible(dev.off())
  
  # convert genotype table to vector of nucleotide sequences
  gtseq_mouse <- genotypeFasta(gt_mouse, germline_db=ighv, novel=nv)
  writeFasta(gtseq_mouse, paste0(genotype_dir, '/IGHV_genotype_', m, '.fasta'))
  
  # correct allele calls based on the personalized genotype
  seqdb_mouse <- reassignAlleles(seqdb_mouse, gtseq_mouse)
  #---------
  
  
  ###################################
  ### CALCULATE NEAREST NEIGHBORS ###
  ###################################
  
  # calculate distances to nearest neighbors
  dist_ham <- distToNearest(seqdb_mouse, vCallColumn="V_CALL_GENOTYPED", model="ham", normalize="len", nproc=1)

  # plot distance distribution
  png(filename=paste0(threshold_dir, '/distToNearestNeighbor_', m, '.png'), units='mm', height=150, width=180, res=300)
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
    
    cat('\n\tSkipping automatic threshold estimation. Threshold value set to default:', opt$default_threshold, '\n')
    
  } else if (dist_args[1] %in% c('density','gmm')) {
    
    # Find threshold using the density or gmm (mixture model) methods
    if (dist_args[1] == 'density') {
      output <- findThreshold(dist_ham$DIST_NEAREST, method=dist_args[1])
    } else {
      if (is.na(dist_args[2])) { dist_args[2] <- 'gamma-gamma' }
      output <- findThreshold(dist_ham$DIST_NEAREST, method=dist_args[1], model=dist_args[2])  
    }
    
    if (is.null(output) || is.na(output@threshold)) {  
      cat('\n\tThreshold estimation failed. Reverting to default threshold:', opt$default_threshold, '\n')
    } else {
      cat('\n\tEstimated hamming distance threshold:', round(output@threshold, 5), '\n')
      predicted_thresholds$threshold[predicted_thresholds$mouse == m] <- output@threshold
      png(paste0(threshold_dir, '/distToNearestNeighbor_', paste(dist_args, collapse='_'), '_fit_', m, '.png'), units='mm', height=120, width=150, res=300)
      plot(output, binwidth=0.02, title=paste0('Threshold Prediction [threshold = ', round(output@threshold, 3), ']'))
      invisible(dev.off())
    }
    
  } else {
    stop(paste0('Invalid density_method: "', dist_args[1], '". Valid options are "density", "gmm", or "none".'))
  }
  #---------
  
  
  ###########################
  ### EXPORT DATA TO FILE ###
  ###########################
  writeChangeoDb(seqdb_mouse, paste0(genotype_dir, '/IGHV-genotyped_', m, '.tab'))
  write.csv(predicted_thresholds, file=paste0(threshold_dir, '/predicted_thresholds.csv'), quote=F, row.names=F)
  
}
#---------


##########################
### PRINT SESSION INFO ###
##########################
cat('R SESSION INFO:\n')
Sys.info()
cat('\n\n\n\n')
sessionInfo()
#---------






