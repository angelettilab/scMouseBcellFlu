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

# create directory for results output
invisible(ifelse(!dir.exists(file.path(proj_dir, 'analysis')), dir.create(file.path(proj_dir, 'analysis')), FALSE))
invisible(ifelse(!dir.exists(file.path(paste0(proj_dir,'/analysis'), 'immcantation')), dir.create(file.path(paste0(proj_dir,'/analysis'), 'immcantation')), FALSE))

#####################################################
### FIND NOVEL VDJ SEQ ALLELES AND INFER GENOTYPE ###
#####################################################

# load V-segment germline sequences
ighv <-readIgFasta(paste0(proj_dir, '/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHV.fasta'))

# import ChangeO-formatted sequence database files (heavy chain)
seqdb <- NULL
for (s in samples){
  x <- readChangeoDb(paste0(proj_dir, '/data/VDJ_OTUs/', s, '/', s, '_heavy_parse-select.tab'))
  x$SEQUENCE_ID <- paste0(x$SEQUENCE_ID, '-', s)
  seqdb <- rbind(seqdb, x)
}

# find novel alleles
nv <- findNovelAlleles(seqdb, ighv, germline_min=100)

# Extract and view the rows that contain successful novel allele calls
novel_rows <- selectNovel(nv)
png(filename=paste0(proj_dir, '/analysis/immcantation/novel_seq_plot.png'), height=2500, width=1800, res=300)
plotNovel(seqdb, novel_rows[1, ])
invisible(dev.off())

# infer genotype
gt <- inferGenotype(seqdb, germline_db=ighv, novel=nv, fraction_to_explain = 0.875,
                    gene_cutoff = 1e-04, find_unmutated = TRUE)
png(filename=paste0(proj_dir, '/analysis/immcantation/IGHV_genotype_plot.png'), height=2500, width=1000, res=300)
plotGenotype(gt, gene_sort="position", text_size=8) #, facet_by='ALLELES')
invisible(dev.off())

# convert genotype table to vector of nucleotide sequences
gtseq <- genotypeFasta(gt, germline_db=ighv, novel=nv)
writeFasta(gtseq, paste0(proj_dir, '/analysis/immcantation/IGHV_genotype.fasta'))

# correct allele calls based on a personalized genotype
seqdb <- reassignAlleles(seqdb, gtseq)


#####################################
### CALCULATING NEAREST NEIGHBORS ###
#####################################

# calculate distances to nearest neighbors
dist_ham <- distToNearest(seqdb, vCallColumn="V_CALL_GENOTYPED", model="ham", normalize="len", nproc=1)

png(filename=paste0(proj_dir, '/analysis/immcantation/distToNearestNeighbor.png'), height=1500, width=1800, res=300)
ggplot(subset(dist_ham, !is.na(DIST_NEAREST)), aes(x=DIST_NEAREST)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02)
  # geom_vline(xintercept=0.1, color="firebrick", linetype=2)
invisible(dev.off())

# Find threshold using the density and gmm (mixture model) methods
output_density <- findThreshold(dist_ham$DIST_NEAREST, method="density")
threshold_density <- output_density@threshold
output_gmm <- findThreshold(dist_ham$DIST_NEAREST, method="gmm", model="gamma-gamma")
threshold_gmm <- output_gmm@threshold
output_gmm_norm <- findThreshold(dist_ham$DIST_NEAREST, method="gmm", model="gamma-norm")
threshold_gmm_norm <- output_gmm_norm@threshold

# Plot distance histogram, fitted curves, and optimum threshold
png(filename=paste0(proj_dir, '/analysis/immcantation/distToNearestNeighbor_densityFit.png'), height=1200, width=1500, res=300)
plot(output_density, binwidth=0.02, title=paste0('Density Method [threshold = ', round(threshold_density,3), ']'))
invisible(dev.off())
png(filename=paste0(proj_dir, '/analysis/immcantation/distToNearestNeighbor_gmmFit.png'), height=1200, width=1500, res=300)
plot(output_gmm, binwidth=0.02, title=paste0('GMM Method: gamma-gamma [threshold = ', round(threshold_gmm,3), ']'))
invisible(dev.off())
png(filename=paste0(proj_dir, '/analysis/immcantation/distToNearestNeighbor_gmmNormFit.png'), height=1200, width=1500, res=300)
plot(output_gmm_norm, binwidth=0.02, title=paste0('GMM Method: gamma-normal [threshold = ', round(threshold_gmm_norm,3), ']'))
invisible(dev.off())

###########################
### EXPORT DATA TO FILE ###
###########################
writeChangeoDb(lung0_2, paste0(proj_dir, '/analysis/immcantation/IGHV-genotyped.tab'))



