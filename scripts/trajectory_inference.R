#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nRunning TRAJECTORY ANALYSIS with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object FILE."),
  make_option(c("-m", "--metadata_use"),          type = "character",   metavar="character",   default='none',  help="Column names of the metadata table to plot in the trajectory map."),
  make_option(c("-d", "--reduction_use"),         type = "character",   metavar="character",   default='dm',    help="Dimensionality reduction method to base the trajectories on. It could be a pre-computed slot within your Seurat Object (in case: pca, umap, umap10, tsne) or it will compute additional others (ica, dm)"),
  make_option(c("-r", "--reduction_visualize"),   type = "character",   metavar="character",   default='umap',  help="Dimensionality reduction method to visualize trajectories. It could be a pre-computed slot within your Seurat Object (in case: pca, umap, umap10, tsne) or it will compute additional others (ica, dm)"),
  make_option(c("-n", "--cluster_use"),           type = "character",   metavar="character",   default='none',  help="The cluster of cells to be used for analysis. Should be defined as the clustering name followed by the cluster names to be used, comma-separated. E.g.: 'louvain_0.2,1,2,3,5,6'. A clustering method MUST be specified, though the cluster names can be omitted to include all clusters."),
  make_option(c("-s", "--start_cluster"),         type = "character",   metavar="character",   default='none',  help="Cluster from which the trajectories will start."),
  make_option(c("-x", "--end_cluster"),           type = "character",   metavar="character",   default='none',  help="Cluster at which the trajectories will end"),
  make_option(c("-z", "--diff_testing"),          type = "character",   metavar="character",   default='none',  help="Whether to test for diffential expression across branches"),
  make_option(c("-a", "--assay"),                 type = "character",   metavar="character",   default='RNA',   help="Assay to use for trajectory differential expression. Default is 'RNA'"),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output DIRECTORY.")
)

opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

# opt <- list(
#   Seurat_object_path  ="/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/TEST/analysis/02_cluster/seurat_object.rds",
#   metadata_use        ="dataset,organ"                                                                                        ,
#   reduction_use       ="dm"                                                                                               ,
#   reduction_visualize ="umap"                                                                                                 ,
#   method_use          ="slingshot"                                                                                              ,
#   no_traj_components  ="3"                                                                                                    ,
#   no_of_paths         ="none"                                                                                                 ,
#   cluster_use         ="louvain_0.45"                                                                                                  ,
#   start_cluster       ="none"                                                                                                 ,
#   diff_testing        ="true"                                                                                                 ,
#   assay               ="RNA"                                                                                                  ,
#   output_path         ="/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/TEST/analysis/X_trajectory"                ,
#   help                ="FALSE")

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)


##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
suppressMessages(suppressWarnings({
  library(Seurat)
  library(slingshot)
  library(viridisLite)
  library(tradeSeq)
  library(destiny)
  library(dplyr)
  library(SingleCellExperiment)
}))

# function for generating color hues (mimics default ggplot palette)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#---------

if (opt$cluster_use == 'none') {
  stop('A clustering method MUST be specified (--cluster_use)')
}


#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\n### LOADING Seurat.v3 OBJECT ###\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------



###################################
### SELECT CELLS FROM A CLUSTER ###
###################################
cat("\n### SELECTING CELLS FROM A CLUSTER ###")

clustering_use <- as.character(unlist(strsplit(opt$cluster_use,",")))[1]
if (clustering_use == 'none') {
  stop('A cluster method MUST be specified with the "cluster_use" input!')
} else if (!(clustering_use %in% colnames(DATA@meta.data))) {
  stop(paste0('The specified clustering method "', clustering_use, '" was not found in the Seurat object.'))
}

if (length(unlist(strsplit(opt$cluster_use,","))) >= 2 ){
  clusters_to_select <- as.character(unlist(strsplit(opt$cluster_use,",")))[-1]
  if(sum(clusters_to_select %in% unique(DATA@meta.data[,clustering_use])) == 0){
    cat("\nAll clusters will be used ...\n")
    cells_use <- rownames(DATA@meta.data)
  } else {
    cells_use <- rownames(DATA@meta.data)[factor(DATA@meta.data[,clustering_use]) %in% clusters_to_select]   #Filter out cells with no assigned clusters
    DATA <- SubsetData(DATA,assay = opt$assay,cells = cells_use)}
} else {
  cat("\nAll clusters will be used ...\n")
  cells_use <- rownames(DATA@meta.data) }
#---------



#########################
### DEFINE PARAMETERS ###
#########################
# NOTE: CURRENTLY CAN ONLY HANDLE ONE OPTION FOR EACH REDUCTION_USE AND VISUALIZE
reduction_use <- unlist(strsplit( opt$reduction_use , split = ","))
reduction_vis <- unlist(strsplit( opt$reduction_visualize , split = ","))
#---------



#########################
### RUN DIFFUSION MAP ###
#########################
cat("\n### RUNNING DIFFUSION MAP ###")
if ('dm' %in% reduction_use) {
  if ('dm' %in% names(DATA@reductions)) {
    cat('\nWARNING: An existing diffusion map (DM) reduction was found in the Seurat object, and will be overwritten.\n')
  }
  # dm <- DiffusionMap(DATA@reductions[["pca"]]@cell.embeddings[ , 1:8], k=4, n_eigs=5)
  dm <- DiffusionMap(DATA@reductions[["mnn"]]@cell.embeddings[, 1:8], k=4, n_eigs=5)
  rownames(dm@eigenvectors) <- colnames(DATA)
  DATA@reductions[["dm"]] <- CreateDimReducObject(embeddings=dm@eigenvectors, key="DC_", assay="RNA")
}
#---------



##########################################
### RUN SLINGSHOT TRAJECTORY INFERENCE ###
##########################################
cat("\n### RUNNING SLINGSHOT TRAJECTORY INFERENCE ###")

# extract fields
clustering <- DATA@meta.data[, clustering_use]
dimred_use <- DATA@reductions[[reduction_use]]@cell.embeddings
dimred_vis <- DATA@reductions[[reduction_vis]]@cell.embeddings
counts <- as.matrix( DATA@assays[[opt$assay]]@counts )
if (casefold(opt$start_cluster) != 'none') {
  start_clust <- trimws(unlist(strsplit(opt$start_cluster, ',')))
} else {
  start_clust <- NULL
}
if (casefold(opt$end_cluster) != 'none') {
  end_clust <- trimws(unlist(strsplit(opt$end_cluster, ',')))
} else {
  end_clust <- NULL
}

# define cell lineages with Singshot
set.seed(1)
lineages <- getLineages(data=dimred_use, clusterLabels=clustering, start.clus=start_clust, end.clus=end_clust)
# lineages@reducedDim <- DATA@reductions$dm@cell.embeddings[,c(1,3)]
lineages@reducedDim <- DATA@reductions$umap@cell.embeddings
lineages

# plot the lineages
pal <- gg_color_hue(nlevels(clustering))
png(filename=paste0(opt$output_path,'/trajectory_lineages.png'), res=300, units='mm', width=220, height=120)
par(mfrow = c(1, 2))
plot(dimred_vis[, 1:2], col=pal[clustering], cex=0.2, pch=16, las=1)
for (i in levels(clustering)) {
  text(mean(dimred_vis[clustering == i, 1]), mean(dimred_vis[clustering == i, 2]), labels=i, font=2)
}
plot(dimred_vis[, 1:2], col=pal[clustering], cex=0.2, pch=16, las=1)
lines(lineages, lwd=1.5, col="black", cex=1)
dev.off()

# define the principal curves for each lineage
curves <- getCurves(lineages, thresh=0.001, stretch=0.01, allow.breaks=F)
curves

# plot principal curves
png(filename=paste0(opt$output_path,'/trajectory_curves.png'), res=300, units='mm', width=120, height=120)
par(mfrow=c(1,1))
plot(dimred_vis, col='grey75', cex=0.2, pch=16)
lines(curves, lwd = 2, col = gg_color_hue(length(curves@curves)))
dev.off()


###########################################
### FIND DIFFERENTIALLY EXPRESSED GENES ###
###########################################

# Fit an additive model (GAM) to the trajectory curves
# ====================================================

# limit to highly variable genes to avoid excessive calculation times
filt_counts <- counts[rownames(counts) %in% DATA@assays[[opt$assay]]@var.features[1:300], ]

# find differentially expressed genes (this takes a while to run)
cat("\n### FITTING GAM ###\n")
sce <- fitGAM(counts=filt_counts, sds=curves, nknots=10)
saveRDS(sce, file=paste0(opt$output_path,'/sce_object.rds'))


# 
# # plot curves
# plotGeneCount(curves, filt_counts, clusters=clustering, models=sce)
# 
# # define plotting function
# plot_differential_expression <- function(feature_id) {
#   cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene=feature_id[1], clusters=clustering, models=sce) +
#                        ggplot2::theme(legend.position = "none"), 
#                      plotSmoothers(sce, as.matrix(counts), gene=feature_id[1]))
# }
# 
# 
# # Identify genes associated with pseudotime
# # =========================================
# 
# # calculate pseudotime association significance for each gene
# rowData(sce)$tradeSeq$beta <- list(rowData(sce)$tradeSeq$beta)  # necessary to avoid weird bug
# pseudotime_association <- associationTest(sce)
# pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method="fdr")
# pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
# pseudotime_association$feature_id <- rownames(pseudotime_association)
# 
# # plot gene with most significant association to pseudotime
# feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, waldStat) %>% pull(feature_id)
# plot_differential_expression(feature_id)
# 
# 
# 
# # Compare specific pseudotime values within a lineage
# # ===================================================
# 
# pseudo_start_end <- startVsEndTest(sce, pseudotimeValues = c(0, 1), lineages=T)
# feature_id <- rownames(pseudo_start_end)[which.max(pseudo_start_end$waldStat)]
# plot_differential_expression(feature_id)
# 
# 
# # find genes associated with early differences in lineages
# earlyDERes <- earlyDETest(sce, knots = c(2, 3))
# oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
# head(rownames(earlyDERes)[oEarly])
# 
# plot_differential_expression(rownames(earlyDERes)[oEarly][1])
# 
# 
# # cluster genes according to their expression pattern
# library(clusterExperiment)
# nPointsClus <- 20
# clusPat <- clusterExpressionPatterns(sce, nPoints=nPointsClus, genes=rownames(filt_counts))
# clusterLabels <- primaryCluster(clusPat$rsec)
# 
# # to fix missing rownames
# rownames(clusPat$yhatScaled) <- rownames(filt_counts)
# 
# # plot clustered genes (show four clusters)
# library(cowplot)
# library(ggplot2)
# 
# cUniq <- unique(clusterLabels)
# cUniq <- cUniq[!cUniq == -1] # remove unclustered genes
# plots <- list()
# for (xx in cUniq[1:4]) {
#   cId <- which(clusterLabels == xx)
#   p <- ggplot(data = data.frame(x = 1:nPointsClus,
#                                 y = rep(range(clusPat$yhatScaled[cId, ]),
#                                         nPointsClus / 2)),
#               aes(x = x, y = y)) +
#     geom_point(alpha = 0) +
#     labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") +
#     theme_classic()
#   for (ii in 1:length(cId)) {
#     geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
#     p <- p +
#       geom_line(data = data.frame(x = rep(1:nPointsClus, 3),
#                                   y = clusPat$yhatScaled[geneId, ],
#                                   lineage = rep(0:2, each = nPointsClus)),
#                 aes(col = as.character(lineage), group = lineage), lwd = 0.5)
#   }
#   p <- p + guides(color = FALSE) +
#     scale_color_manual(values = c("orange", "darkseagreen3", "purple"),
#                        breaks = c("0", "1", "2"))  
#   plots[[as.character(xx)]] <- p
# }
# plots$ncol <- 2
# do.call(plot_grid, plots)
# 
# 


