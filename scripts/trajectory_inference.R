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
  make_option(c("-p", "--destiny_params"),        type = "character",   metavar="character",   default='default',  help="Optional adjustment of parameters used for calculating the diffusion map components with destiny. Should be provided as 'param=value' separated by commas, where available params are 'k' and 'n_eigs'. For example, the default settings are equivalent to 'k=30, n_eigs=20'."),
  make_option(c("-n", "--cluster_use"),           type = "character",   metavar="character",   default='none',  help="The cluster of cells to be used for analysis. Should be defined as the clustering name followed by the cluster names to be used, comma-separated. E.g.: 'louvain_0.2,1,2,3,5,6'. A clustering method MUST be specified, though the cluster names can be omitted to include all clusters."),
  make_option(c("-s", "--start_cluster"),         type = "character",   metavar="character",   default='none',  help="Cluster from which the trajectories will start."),
  make_option(c("-x", "--end_cluster"),           type = "character",   metavar="character",   default='none',  help="Cluster(s) at which the trajectories will end. By default, additional trajectories ending at other clusters may also be produced; include 'only' to exclude these additional trajectories. For example: '3, only' "),
  make_option(c("-z", "--diff_testing"),          type = "character",   metavar="character",   default='true',  help="Whether to test for diffential expression within lineages"),
  make_option(c("-a", "--assay"),                 type = "character",   metavar="character",   default='RNA',   help="Assay to use for trajectory differential expression. Default is 'RNA'"),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output DIRECTORY.")
)

opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

# opt <- list(Seurat_object_path  ="/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/analysis/04_cluster/seurat_object.rds",
#             metadata_use        ="none"                                                                        ,
#             reduction_use       ="dm"                                                                          ,
#             reduction_visualize ="umap"                                                                        ,
#             destiny_params      ="k=30, n_eigs=20"                                                             ,
#             cluster_use         ="louvain_0.95,2,4,5,7,8,9,11,12,13,14,15"                                     ,
#             start_cluster       ="13"                                                                          ,
#             end_cluster         ="14"                                                                          ,
#             diff_testing        ="true"                                                                        ,
#             assay               ="RNA"                                                                         ,
#             output_path         ="/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/analysis/trajectory_01b",
#             help                ="FALSE")

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
  library(ggplot2)
  library(cowplot)
}))

# function for generating color hues (mimics default ggplot palette)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# more basic gray to blue color scale
blue_pal <- c("grey85","navy")
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
cat("\n### SELECTING CELLS FROM A CLUSTER ###\n")

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
    DATA <- subset(DATA, cells = cells_use)}
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
cat("\n### RUNNING DIFFUSION MAP ###\n")
if ('dm' %in% reduction_use) {
  dm_filename <- paste0(opt$output_path,'/diffusion_map_coords.csv')
  if (file.exists(dm_filename)) {
    cat('Existing diffusion map reduction found (', dm_filename, ') - skipping calculation.\n\n', sep='')
    dm_coords <- as.matrix(read.csv(dm_filename, row.names=1))
    dm_coords <- dm_coords[match(colnames(DATA), rownames(dm_coords)), ]  # ensure same cell ordering
    DATA@reductions[["dm"]] <- CreateDimReducObject(embeddings=dm_coords, key="DC_", assay=opt$assay)
    
  } else {
    
    # get parameters for diffusion map calculation
    dim_reduct_params <- list(k=30, n_eigs=20)  # set defaults
    if (!(opt$destiny_params %in% c('default','defaults'))) {
      new_params <- eval(parse(text=paste0("list(", opt$destiny_params, ")")))
      if (!all(names(new_params) %in% names(dim_reduct_params))) {
        stop(paste('Invalid dim_reduct_param parameter(s):', paste(names(new_params)[!(names(new_params) %in% names(dim_reduct_params))], collapse=', ')))
      }
      dim_reduct_params <- modifyList(dim_reduct_params, new_params)
    }
    
    # dm <- DiffusionMap(DATA@reductions[["pca"]]@cell.embeddings[ , 1:8], k=4, n_eigs=5)
    dm <- DiffusionMap(DATA@reductions[["mnn"]]@cell.embeddings, k=dim_reduct_params$k, n_eigs=dim_reduct_params$n_eigs, verbose=T)
    rownames(dm@eigenvectors) <- colnames(DATA)
    DATA@reductions[["dm"]] <- CreateDimReducObject(embeddings=dm@eigenvectors, key="DC_", assay=opt$assay)
    
    cat('Saving diffusion map reduction coordinates... ')
    write.csv(DATA@reductions$dm@cell.embeddings, dm_filename)
    cat('Done.\n\n')
  }
}
#---------



##########################################
### RUN SLINGSHOT TRAJECTORY INFERENCE ###
##########################################
cat("\n### RUNNING SLINGSHOT TRAJECTORY INFERENCE ###\n")

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
allow_more_lineages <- TRUE  # default
if (casefold(opt$end_cluster) != 'none') {
  end_clust <- trimws(unlist(strsplit(opt$end_cluster, ',')))
  if ('only' %in% end_clust) {
    # WARNING: REMOVING LINEAGES IN THIS WAY DOES NOT SEEM TO WORK!
    # It causes problems with fitGAMs later on
    allow_more_lineages <- FALSE  # restrict to the requested end points
    end_clust <- end_clust[!(end_clust %in% 'only')]
  }
} else {
  end_clust <- NULL
}

# define cell lineages with Singshot
if (file.exists(paste0(opt$output_path,'/lineages_object.rds'))) {
  cat('Found existing lineages object - skipping lineages calculation.\n')
  lineages <- readRDS(paste0(opt$output_path,'/lineages_object.rds'))
} else {
  set.seed(1)
  cat('Defining cell lineages with Slingshot... ')
  lineages <- getLineages(data=dimred_use, clusterLabels=clustering, start.clus=start_clust, end.clus=end_clust)
  # lineages@reducedDim <- DATA@reductions$dm@cell.embeddings[,c(1,3)]
  
  # remove additional lineages if requested
  if (!allow_more_lineages) {
    keep_lineages <- lineages@slingParams$end.given
    lineages@lineages[!keep_lineages] <- NULL
  }
  
  lineages@reducedDim <- dimred_vis
  saveRDS(lineages, file=paste0(opt$output_path,'/lineages_object.rds'))
  cat('Done.\n')
}

# plot the lineages
pal <- gg_color_hue(nlevels(clustering))
png(filename=paste0(opt$output_path,'/trajectory_lineages.png'), res=300, units='mm', width=120, height=120)
plot(dimred_vis[, 1:2], col=pal[clustering], cex=0.2, pch=16, las=1)
lines(lineages, lwd=2, col="black", cex=1.7)
invisible(lapply(levels(clustering), function(x) text(mean(dimred_vis[clustering == x, 1]), mean(dimred_vis[clustering == x, 2]), labels=x, font=1, cex=0.6, col='white')))
invisible(dev.off())

# define the principal curves for each lineage
if (file.exists(paste0(opt$output_path,'/curves_object.rds'))) {
  cat('Found existing curves object - skipping curves calculation.\n')
  curves <- readRDS(paste0(opt$output_path,'/curves_object.rds'))
} else {
  cat('Defining the principal curves for each lineage... ')
  curves <- getCurves(lineages, thresh=0.01, stretch=0.0001, allow.breaks=T, extend='n')
  saveRDS(curves, file=paste0(opt$output_path,'/curves_object.rds'))
  cat('Done.\n\n')
}
curves  # show the curve data


# get curve points and endpoints (used for plotting)
curve_points <- NULL
curve_ends <- NULL
for (i in 1:length(curves@curves)) {
  temp <- as.data.frame(curves@curves[[i]]$s[curves@curves[[i]]$ord, ])
  curve_ends <- rbind(curve_ends, tail(temp, n=1))
  temp$curve_num <- i
  curve_points <- rbind(curve_points, temp)
}
colnames(curve_points)[1:2] <- c('x','y')
colnames(curve_ends) <- c('x','y')
rownames(curve_ends) <- 1:nrow(curve_ends)

# set up curve color palette
gray_pal <- gray.colors(length(curves@curves), start=0, end=0.4)

# plot principal curves
png(filename=paste0(opt$output_path,'/trajectory_curves.png'), res=300, units='mm', width=120, height=120)
plot(rbind(dimred_vis, as.matrix(curve_ends)), col=pal[clustering], cex=0.2, pch=16)
lines(curves, lwd=2, col=gray_pal)
points(curve_ends, pch=16, cex=1.7)
text(curve_ends, labels=rownames(curve_ends), cex=0.6, col='white')
invisible(dev.off())

# plot principal curves, colored by pseudotime
pseudo_time <- slingPseudotime(curves)
curve_labels <- colnames(pseudo_time)
n_plot_cols <- min(3, length(curve_labels))
n_plot_rows <- ceiling(length(curve_labels)/n_plot_cols)
pal <- viridis(100, end=0.95)
png(filename=paste0(opt$output_path,'/trajectory_curves_pseudotime.png'), res=300, units='mm', width=120*n_plot_cols, height=120*n_plot_rows)
par(mfrow=c(n_plot_rows, n_plot_cols), mgp=c(0.5, 0, 0), mar=c(2,2,2,2), cex=1.3)
for (i in curve_labels) {
  colors <- pal[cut(pseudo_time[,i], breaks=100)]
  # use "reducedDims(curves)" in case its ordering differs from dimred_vis
  plot(reducedDims(curves), xlim=range(dimred_vis[,1]), ylim=range(dimred_vis[,2]), col='grey85', pch=16, cex=0.3, yaxt='n', xaxt='n')
  par(new=T)
  plot(reducedDims(curves), xlim=range(dimred_vis[,1]), ylim=range(dimred_vis[,2]), col=colors, pch=16, cex=0.3, main=i, yaxt='n', xaxt='n')
  tmp <- curves; tmp@curves <- tmp@curves[i]
  lines(tmp, lwd=3, col='black', type='curves')
  box(lwd=1.5)
}
invisible(dev.off())
write.csv(pseudo_time, file=paste0(opt$output_path,'/trajectory_curves_pseudotime.csv'))


###########################################
### FIND DIFFERENTIALLY EXPRESSED GENES ###
###########################################

if (!(casefold(opt$diff_testing) %in% c('true','yes'))) {
  cat('\nSkipping lineage differential expression analysis as requested!\n\n')
  return()
}

# Fit an additive model (GAM) to the trajectory curves
# ====================================================

if (file.exists(paste0(opt$output_path,'/sce_object.rds'))) {
  
  # load existing SCE object
  cat('Found existing SCE object - skipping GAM fitting.\n')
  sce <- readRDS(paste0(opt$output_path,'/sce_object.rds'))
  filt_counts <- sce@assays@data@listData$counts
  
} else {
  # limit to highly variable genes to avoid excessive calculation times
  filt_counts <- counts[rownames(counts) %in% DATA@assays[[opt$assay]]@var.features[1:500], ]
  
  # find differentially expressed genes (this takes a while to run)
  cat("\n### FITTING GAM ###\n")
  sce <- fitGAM(counts=filt_counts, sds=curves, nknots=10)
  saveRDS(sce, file=paste0(opt$output_path,'/sce_object.rds'))
}

# plot curves
png(filename=paste0(opt$output_path,'/trajectory_curves_tradeseq.png'), res=300, units='mm', width=120, height=120)
plotGeneCount(curves, clusters=clustering, models=sce) +
  theme(legend.position='none')
invisible(dev.off())

# define plotting function
plot_feat_curves <- function(feature_id, curve_id=1:length(curves@curves)) {
  # feature_id - name(s) of feature(s) (genes) whose expression should be plotted
  # curve_id - index of curve(s) to show (default = show all curves)
  
  # select curves
  sel_curves <- curves
  sel_curves@curves <- sel_curves@curves[curve_id]
  sel_curve_ends <- curve_ends[curve_id, ]
  
  # plot cells with mapped gene expression
  color_pal <- c("grey90","grey80","royalblue","navy","navy")
  n_plot_col = ceiling(sqrt(length(feature_id)))
  n_plot_row = ceiling(length(feature_id)/n_plot_col)
  par(mgp=c(0.5, 0, 0), mfrow=c(n_plot_row, n_plot_col), mar=c(2,2,2,2))
  
  for (i in 1:length(feature_id)) {
    feat <- filt_counts[feature_id[i], ]
    feat <- feat / ( sort(feat, T, na.last=T)[ min(10, sum(feat!=0, na.rm=T)) ] )
    feat[feat > 1] <- 1
    o <- order(feat, na.last=T)
    feat_col <- c( color_pal[1],colorRampPalette(color_pal[-1])(10))[round(feat*9)+1][o]
    
    plot(rbind(dimred_vis[o, ], as.matrix(sel_curve_ends)), col=feat_col, cex=0.2, pch=16, yaxt='n', xaxt='n')
    title(feature_id[i], adj=0, line=0.5, cex=0.8)
    lines(sel_curves, lwd=1.5, col='black')
    points(sel_curve_ends, pch=16, cex=1.7)
    text(sel_curve_ends, labels=rownames(sel_curve_ends), cex=0.6, col='white')
  }
}


# Identify genes associated with pseudotime
# =========================================

# calculate overall pseudotime association significance for each gene
rowData(sce)$tradeSeq$beta <- list(rowData(sce)$tradeSeq$beta)  # necessary to avoid weird bug
pseudotime_association <- associationTest(sce)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method="fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$waldStat, decreasing=T), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)

# plot gene with most significant association to pseudotime overall (all lineages)
feature_id <- pseudotime_association %>% top_n(1, waldStat) %>% pull(feature_id)
png(filename=paste0(opt$output_path,'/top_pseudotime_assoc_gene_umap.png'), res=300, units='mm', width=120, height=120)
plot_feat_curves(feature_id)
invisible(dev.off())

# due to a bug in plotSmoothers, cannot plot more than 9 lineages
if (length(curves@curves) < 10) {
  png(filename=paste0(opt$output_path,'/top_pseudotime_assoc_gene_traject_expr.png'), res=300, units='mm', width=120, height=120)
  plotSmoothers(sce, as.matrix(filt_counts), gene=feature_id) + ggtitle(feature_id)
  invisible(dev.off())
}


# Calculate pseudotime association for each lineage separately
lineage_association <- associationTest(sce, lineages=T)
lineage_association$fdr <- p.adjust(lineage_association$pvalue, method="fdr")
lineage_association <- lineage_association[order(lineage_association$waldStat, decreasing=T), ]
lineage_association$feature_id <- rownames(lineage_association)

# export lineage association data to csv
write.csv(lineage_association, file=paste0(opt$output_path,'/lineage_association.csv'))

cat('Plotting lineage expression profiles... ')
for (L in 1:length(curves@curves)) {
  o <- order(lineage_association[, paste0('waldStat_', L)], decreasing=T)
  top_feat <- rownames(lineage_association[o, ])[1:6]
  
  png(filename=paste0(opt$output_path,'/curve',L,'_assoc_gene_umap.png'), res=300, units='mm', width=200, height=120)
  plot_feat_curves(top_feat, L)
  invisible(dev.off())
  
  if (length(curves@curves) < 10) {
    plot_list <- lapply(top_feat, function(x) {plotSmoothers(sce, as.matrix(filt_counts), gene=x, lwd=1, size=0.1) + ggtitle(x)})
    p <- cowplot::plot_grid(plotlist=plot_list, ncol=3)
    ggplot2::ggsave(p, filename=paste0('/curve',L,'_assoc_gene_traject_expr.png'), path=opt$output_path,
                    dpi=300, units='mm', width=300, height=150, limitsize=F)
  }
}
cat('Done.\n\n')


# Compare specific pseudotime values within a lineage
# ===================================================

pseudo_start_end <- startVsEndTest(sce, pseudotimeValues = c(0, 1), lineages=T)
feature_id <- rownames(pseudo_start_end)[which.max(pseudo_start_end$waldStat)]

png(filename=paste0(opt$output_path,'/top_start_end_assoc_gene_umap.png'), res=300, units='mm', width=120, height=120)
plot_feat_curves(feature_id)
invisible(dev.off())

if (length(curves@curves) < 10) {
  png(filename=paste0(opt$output_path,'/top_start_end_assoc_gene_traject_expr.png'), res=300, units='mm', width=120, height=120)
  plotSmoothers(sce, as.matrix(filt_counts), gene=feature_id) + ggtitle(feature_id)
  invisible(dev.off())
}

# find genes associated with early differences in lineages
earlyDERes <- earlyDETest(sce, knots = c(3, 5))
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
# head(rownames(earlyDERes)[oEarly])
feature_id <- rownames(earlyDERes)[oEarly][1]

png(filename=paste0(opt$output_path,'/top_early_diff_assoc_gene_umap.png'), res=300, units='mm', width=120, height=120)
plot_feat_curves(feature_id)
invisible(dev.off())

if (length(curves@curves) < 10) {
  png(filename=paste0(opt$output_path,'/top_early_diff_assoc_gene_traject_expr.png'), res=300, units='mm', width=120, height=120)
  plotSmoothers(sce, as.matrix(filt_counts), gene=feature_id) + ggtitle(feature_id)
  invisible(dev.off())
}


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
# nLin <- length(curves@lineages)
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
#       geom_line(data = data.frame(x = rep(1:nPointsClus, nLin),
#                                   y = clusPat$yhatScaled[geneId, ],
#                                   lineage = rep(0:(nLin-1), each = nPointsClus)),
#                 aes(col = as.character(lineage), group = lineage), lwd = 0.5)
#   }
#   p <- p + guides(color = FALSE) +
#     scale_color_manual(values = gg_color_hue(nLin),
#                        breaks = as.character(0:(nLin-1)))
#   plots[[as.character(xx)]] <- p
# }
# plots$ncol <- 2
# png(filename=paste0(opt$output_path,'/clustered_trajectory_genes.png'), res=300, units='mm', width=240, height=240)
# do.call(plot_grid, plots)
# invisible(dev.off())




