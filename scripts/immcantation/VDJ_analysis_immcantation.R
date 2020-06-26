# Script for B-cell VDJ analysis using Immcantation package

library(alakazam)
library(shazam)
library(tigger)


#################################
### SPECIFY SCRIPT PARAMETERS ###
#################################

proj_dir <- '/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'

# get list of all available samples
samples <- list.dirs(paste0(proj_dir, '/data/VDJ_OTUs'), full.names=F, recursive=F)

# get sample metadata
sample.meta <- read.csv(paste0(proj_dir, '/data/metadata.csv'), stringsAsFactors=F)
sample.meta <- sample.meta[match(samples, sample.meta$dataset), ]  # align to sample list
rownames(sample.meta) <- seq(nrow(sample.meta))

# create directories for results output
invisible(ifelse(!dir.exists(file.path(proj_dir, 'analysis')), dir.create(file.path(proj_dir, 'analysis')), FALSE))
invisible(ifelse(!dir.exists(file.path(paste0(proj_dir,'/analysis'), 'immcantation')), dir.create(file.path(paste0(proj_dir,'/analysis'), 'immcantation')), FALSE))
invisible(ifelse(!dir.exists(file.path(paste0(proj_dir,'/analysis/immcantation'), 'genotyping')), dir.create(file.path(paste0(proj_dir,'/analysis/immcantation'), 'genotyping')), FALSE))
invisible(ifelse(!dir.exists(file.path(paste0(proj_dir,'/analysis/immcantation'), 'threshold_estimation')), dir.create(file.path(paste0(proj_dir,'/analysis/immcantation'), 'threshold_estimation')), FALSE))


#####################################################
### FIND NOVEL VDJ SEQ ALLELES AND INFER GENOTYPE ###
#####################################################

# load V-segment germline sequences
ighv <-readIgFasta(paste0(proj_dir, '/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHV.fasta'))

# import ChangeO-formatted sequence database files (heavy chain)
seqdb <- NULL
for (i in seq(length(samples))){
  x <- readChangeoDb(paste0(proj_dir, '/data/VDJ_OTUs/', samples[i], '/', samples[i], '_heavy_parse-select.tab'))
  x$SEQUENCE_ID <- paste0(x$SEQUENCE_ID, '-', samples[i])
  x$MOUSE_NR <- sample.meta$mouse_nr[i]
  seqdb <- rbind(seqdb, x)
}

# infer genotype (should be done separately for each mouse)
mice <- unique(seqdb$MOUSE_NR)
for (m in mice) {
  
  seqdb_mouse <- seqdb[seqdb$MOUSE_NR %in% m, ]
  
  # find novel alleles (if any)
  novel_rows <- NULL
  try (nv <- findNovelAlleles(seqdb_mouse, ighv))
  try (novel_rows <- selectNovel(nv))
  
  # Extract and view the rows that contain successful novel allele calls
  
  if (!is.null(novel_rows) && (nrow(novel_rows) > 0)) {
    png(filename=paste0(proj_dir, '/analysis/immcantation/genotyping/novel_alleles_', m, '.png'), height=2500, width=1800, res=300)
    plotNovel(seqdb_mouse, novel_rows[1, ])  # only plot first novel allele
    invisible(dev.off())
  }
  
  # infer mouse genotype
  gt_mouse <- inferGenotype(seqdb_mouse, germline_db=ighv, novel=nv)
  png(filename=paste0(proj_dir, '/analysis/immcantation/genotyping/IGHV_genotype_plot_', m, '.png'), units='mm', height=250, width=100, res=300)
  plotGenotype(gt_mouse, gene_sort="position", text_size=8) #, facet_by='ALLELES')
  invisible(dev.off())
  
  # convert genotype table to vector of nucleotide sequences
  gtseq_mouse <- genotypeFasta(gt_mouse, germline_db=ighv, novel=nv)
  writeFasta(gtseq_mouse, paste0(proj_dir, '/analysis/immcantation/genotyping/IGHV_genotype_', m, '.fasta'))
  
  # correct allele calls based on the personalized genotype
  seqdb_mouse <- reassignAlleles(seqdb_mouse, gtseq_mouse)
  
  
  ###################################
  ### CALCULATE NEAREST NEIGHBORS ###
  ###################################
  
  # calculate distances to nearest neighbors
  dist_ham <- distToNearest(seqdb_mouse, vCallColumn="V_CALL_GENOTYPED", model="ham", normalize="len", nproc=1)

  # plot distance distribution
  png(filename=paste0(proj_dir, '/analysis/immcantation/threshold_estimation/distToNearestNeighbor_', m, '.png'), units='mm', height=150, width=180, res=300)
  print(ggplot(subset(dist_ham, !is.na(DIST_NEAREST)), aes(x=DIST_NEAREST)) +
          theme_bw() +
          xlab("Hamming distance") +
          ylab("Count") +
          scale_x_continuous(breaks=seq(0, 1, 0.1)) +
          geom_histogram(color="white", binwidth=0.02))
  invisible(dev.off())
  
  
  # Find threshold using the density and gmm (mixture model) methods
  output_density <- findThreshold(dist_ham$DIST_NEAREST, method="density")
  threshold_density <- output_density@threshold

  png(filename=paste0(proj_dir, '/analysis/immcantation/threshold_estimation/distToNearestNeighbor_densityFit_', m, '.png'), units='mm', height=120, width=150, res=300)
  plot(output_density, binwidth=0.02, title=paste0('Density Method [threshold = ', round(threshold_density,3), ']'))
  invisible(dev.off())
  
  # output_gmm <- findThreshold(dist_ham$DIST_NEAREST, method="gmm", model="gamma-gamma")
  # if (!is.null(output_gmm)) {
  #   threshold_gmm <- output_gmm@threshold
  #   png(filename=paste0(proj_dir, '/analysis/immcantation/threshold_estimation/distToNearestNeighbor_gmmFit_', m, '.png'), units='mm', height=120, width=150, res=300)
  #   plot(output_gmm, binwidth=0.02, title=paste0('GMM Method: gamma-gamma [threshold = ', round(threshold_gmm,3), ']'))
  #   invisible(dev.off())
  # } else {
  #   threshold_gmm <- NULL
  # }
  # output_gmm_norm <- findThreshold(dist_ham$DIST_NEAREST, method="gmm", model="gamma-norm")
  # if (!is.null(output_gmm_norm)) {
  #   threshold_gmm_norm <- output_gmm_norm@threshold
  #   png(filename=paste0(proj_dir, '/analysis/immcantation/threshold_estimation/distToNearestNeighbor_gmmNormFit_', m, '.png'), units='mm', height=120, width=150, res=300)
  #   plot(output_gmm_norm, binwidth=0.02, title=paste0('GMM Method: gamma-normal [threshold = ', round(threshold_gmm_norm,3), ']'))
  #   invisible(dev.off())
  # } else {
  #   threshold_gmm_norm <- NULL
  # }

  
  ###########################
  ### EXPORT DATA TO FILE ###
  ###########################
  writeChangeoDb(seqdb_mouse, paste0(proj_dir, '/analysis/immcantation/genotyping/IGHV-genotyped_', m, '.tab'))

}










