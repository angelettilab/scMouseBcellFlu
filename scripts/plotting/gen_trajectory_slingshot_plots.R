# This script generates some additional plots from the slingshot/tradeSeq trajectory inference analysis results.
# Use conda environment: trajectory-env

# load packages
suppressMessages(suppressWarnings({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(slingshot)
}))

# specify relevant directories
pdir <- '/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/'  # project directiory
adir <- paste0(pdir,'analysis/')  # analysis subdirectory


############################
### TRAJECTORY INFERENCE ###
############################

# specify directories and clustering type(s)
cdir <- paste0(adir, '06_cluster/')  # clustering analysis direcotry
parent_tdir <- paste0(adir, 'trajectory_slingshot/')  # slingshot trajectory analysis directory
tdirs <- dir(parent_tdir, '^trajectory_\\d')  # individual trajectory run sub-directories
clustering_use <- 'HC_16'  # cluster grouping to show on plots

# initialize summary table
table_cols <- c('ID', 'Clustering', 'Clusters', 'Start', 'End')
traj_summary_tab <- data.frame(matrix(nrow=length(tdirs), ncol=length(table_cols))) %>% setNames(table_cols)
traj_summary_tab$ID <- tdirs
traj_summary_tab$Clustering <- clustering_use  # same clustering for all

# load Seurat object
DATA <- readRDS(paste0(cdir, 'seurat_object.rds'))

# loop through each trajectory result subdirectory
for (i in seq(length(tdirs))) {
  tdir <- paste0(adir, 'trajectory_slingshot/', tdirs[i], '/')
  
  # load lineages object and extract run information
  lineages <- readRDS(paste0(tdir, 'lineages_object.rds'))
  clusters <- colnames(lineages@clusterLabels)
  traj_summary_tab$Clusters[i] <- clusters[order(as.numeric(clusters))] %>% paste0(collapse=',')
  if (any(lineages@slingParams$start.given)) {
    start_clusters <- lineages@slingParams$start.clus[lineages@slingParams$start.given]
    traj_summary_tab$Start[i] <- start_clusters[order(as.numeric(start_clusters))] %>% paste0(collapse=',')
  } else {
    traj_summary_tab$Start[i] <- 'auto'
  }
  if (any(lineages@slingParams$end.given)) {
    end_clusters <- lineages@slingParams$end.clus[lineages@slingParams$end.given]
    traj_summary_tab$End[i] <- end_clusters[order(as.numeric(end_clusters))] %>% paste0(collapse=',')
  } else {
    traj_summary_tab$End[i] <- 'auto'
  }
  
  # load diffusion map coordinates
  dm_filename <- paste0(tdir, 'diffusion_map_coords.csv')
  dm_coords <- as.matrix(read.csv(dm_filename, row.names=1))
  
  # subset Seurat object and add diffusion map reduction
  DATAsub <- subset(DATA, cells=rownames(dm_coords))
  dm_coords <- dm_coords[match(colnames(DATAsub), rownames(dm_coords)), ]  # ensure same cell ordering
  DATAsub@reductions[["dm"]] <- CreateDimReducObject(embeddings=dm_coords, key="DC_", assay='RNA')
  
  # generate plot showing UMAP and first several diffusion components
  umap_plot <- list(DimPlot(DATAsub, dims=1:2, reduction='umap', group.by=clustering_use, pt.size=0.1, label=T) +
                      theme(legend.position='none'))
  plot_nums <- as.list(seq(min(5, ncol(DATA@reductions$dm)-1)))
  plot_list <- lapply(plot_nums, function(x) { 
    DimPlot(DATAsub, dims=x:(x+1), reduction='dm', group.by=clustering_use, pt.size=0.1, ncol=3, label=T) +
      ggplot2::theme(legend.position='none')})
  p <- cowplot::plot_grid(plotlist=c(umap_plot, plot_list), ncol=3)
  ggplot2::ggsave(p, filename='diffusion_map_components.png', path=tdir, dpi=300,
                  units='mm', width=120*3, height=100*2, limitsize=F)
  
  # for trajectory runs with specified start and end points
  if (!(traj_summary_tab$Start[i] == 'auto') & !(traj_summary_tab$End[i] == 'auto')) {
    
    # generate UMAP with pseudotime color overlay for specified trajectory curve
    curves <- readRDS(paste0(tdir, 'curves_object.rds'))
    pseudo_time <- read.csv(paste0(tdir, 'trajectory_curves_pseudotime.csv'), row.names=1)
    keep_curves <- which(lineages@slingParams$end.given)
    
    for (kcurve in keep_curves) {
      pal <- viridis(100, end=0.95)
      dimred <- reducedDims(curves)
      colors <- pal[cut(pseudo_time[,kcurve], breaks=100)]
      tmp <- curves; tmp@curves <- tmp@curves[kcurve]
      
      png(filename=paste0(tdir, 'trajectory_curve', kcurve, '_pseudotime.png'), res=300, units='mm', width=120, height=120)
      par(mgp=c(0.5, 0, 0), mar=c(2,2,2,2), cex=1.3)
      plot(dimred, xlim=range(dimred[,1]), ylim=range(dimred[,2]), col='grey85', pch=16, cex=0.3, yaxt='n', xaxt='n')
      par(new=T)
      plot(dimred, xlim=range(dimred[,1]), ylim=range(dimred[,2]), col=colors, pch=16, cex=0.3, main='Trajectory Pseudotime', yaxt='n', xaxt='n')
      lines(tmp, lwd=3, col='black', type='curves')
      box(lwd=1.5)
      invisible(dev.off())
    }
  }
}

# write trajectory summary table to csv
write.csv(traj_summary_tab, file=paste0(parent_tdir, 'trajectory_runs_summary.csv'))
