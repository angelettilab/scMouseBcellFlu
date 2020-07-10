trajectory_pseudotime_heatmap <- function(sce, clustering, lineage_num=1, n_top_genes=50, out_file=NULL) {

library(Seurat)
library(slingshot)
library(viridisLite)
library(tradeSeq)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
# library(heatmap3)
library(SingleCellExperiment)

# sce <- readRDS('/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/analysis/trajectory_slingshot/trajectory_05/sce_object.rds')
# DATA <- readRDS('/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/analysis/06_cluster/seurat_object.rds')
# DATA <- subset(DATA, cells=sce$slingshot@rownames)
# clustering <- as.factor(DATA@meta.data[, 'HC_16']) %>% setNames(colnames(DATA))

  
# function for generating color hues (mimics default ggplot palette)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# pre-process clustering data
# if (is.factor(clustering)) {
#   clustering <- as.numeric(as.character(clustering)) %>% setNames(names(clustering))
# }
clustering <- data.frame(cluster=clustering) %>% rownames_to_column(var='cell')
if (!any(clustering$cell %in% sce$slingshot@rownames)) {
  stop('Names of clustering input should be cell IDs found in the SCE input.')
}

# calculate lineage pseudotime association significance for each gene
if (is.data.frame(rowData(sce)$tradeSeq$beta)) {
  rowData(sce)$tradeSeq$beta <- list(rowData(sce)$tradeSeq$beta)  # necessary bug fix
}
lineage_association <- associationTest(sce, lineages=T)

# get top N genes with most significant association to lineage X
top_genes <- rownames(lineage_association)[order(lineage_association[paste0('waldStat_',lineage_num)], decreasing=T)][1:n_top_genes]

# extract lineage pseudotime (PT) and initialize expr_data
PT <- sce$slingshot[, paste0('pseudotime.curve', lineage_num)] %>% setNames(sce$slingshot@rownames)
PT <- PT[sce$slingshot[, paste0('cellWeights.curve', lineage_num)] > 0]
expr_data <- data.frame(pseudotime=PT) %>% rownames_to_column(var='cell') %>% arrange(pseudotime)

# add corresponding count data
filt_counts <- sce@assays@data@listData$counts
count_data <- as.data.frame(t(filt_counts)[, top_genes]) %>% rownames_to_column(var='cell') %>% inner_join(clustering, by='cell')
expr_data <- inner_join(expr_data, count_data, by='cell') %>% select(-cell)

# calculate smoothed (Loess-fit) curves for each gene
expr_data <- gather(expr_data, variable, value, -pseudotime, -cluster) %>%
  group_by(variable) %>%
  mutate(Loess = loess(value ~ pseudotime, span=0.1)$fitted)

# convert back to wide format and keep only smoothed curves
loess_curves <- expr_data %>%
  select(-value) %>% 
  distinct() %>% 
  pivot_wider(names_from=variable, values_from=Loess)

# reduce number of points along the curves (via interpolation).
# for each pseudotime point set the cluster identity to the majority cluster for that time point
timepoints <- seq(0, max(loess_curves$pseudotime), length.out=500)
timestep <- timepoints[2] - timepoints[1]
Mode <- function(x) {
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}
clusters_interp <- data.frame(cluster=factor(rep(NA, length(timepoints)), levels=levels(clustering$cluster)), row.names=seq(length(timepoints)))
for (i in seq(length(timepoints))) {
  indx <- between(loess_curves$pseudotime, timepoints[i]-timestep/2, timepoints[i]+timestep/2)
  if (any(indx)) {
    clusters_interp$cluster[i] <- Mode(loess_curves$cluster[indx])
  }
}
loess_curves <- select(loess_curves, -cluster) %>% distinct()
curves_interp <- apply(loess_curves[,-1], 2, function(x) approx(loess_curves$pseudotime, x, xout=timepoints)$y) %>% as.matrix() %>% t()
colnames(curves_interp) <- seq(ncol(curves_interp))



# # generate heatmap of gene expression as a function of pseudotime
# heatmap(curves_interp, Colv=NA, distfun=function(x) as.dist(1-cor(t(x))),
#         labCol=NA, labRow=rownames(curves_interp), col=magma(100,direction=-1),
#         keep.dendro=F)


# terribly inefficient way to create vector to label x-axis
x_ticks <- seq(0, 1, 0.1)
keep <- unlist(lapply(x_ticks , function(i) which.min(abs(seq(min(x_ticks), max(x_ticks), length.out=ncol(curves_interp)) - i))))
xvals <- rep("", ncol(curves_interp))
xvals[keep] <- x_ticks

# annotation colors for clusters
annColors <- list(cluster = gg_color_hue(nlevels(clusters_interp$cluster)) %>% setNames(levels(clusters_interp$cluster)))

# curves_scaled <- t(scale(t(curves_interp), center=F, scale=T))

# png('test_heatmap.png', units='mm', width=200, height=200, res=150)
# svg('test_heatmap.svg', width=7, height=7)
p <- pheatmap(curves_interp,
              cluster_cols=F,
              scale='row',
              color=magma(200,direction=-1),
              clustering_distance_rows='correlation',
              breaks=seq(-1, 4, length.out=200),
              cutree_rows=3,
              labels_col=xvals,
              angle_col=0,
              border_color=NA,
              annotation_col=clusters_interp,
              annotation_colors=annColors,
              annotation_names_col=F)

if (!is.null(out_file)) {
  if (endsWith(casefold(out_file), '.svg')) {
    svg(out_file, width=7, height=7)
    print(p)
    invisible(dev.off())
  } else {
    png(out_file, units='mm', width=200, height=200, res=150)
    print(p)
    invisible(dev.off())
  }
}

return(p)


}





# 
# heatmap3(curves_interp,
#          Colv=NA,
#          distfun=function(x) as.dist(1-cor(t(x))),
#          labCol=xvals,
#          labRow=rownames(curves_interp),
#          col=magma(100, direction=-1),
#          xlab='Pseudotime',
#          lasCol=1,
#          ColSideColors=as.character(clusters_interp)
#          )



