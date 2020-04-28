
# load packages
suppressMessages(suppressWarnings({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(viridis)
}))

# specify relevant directories
pdir <- '/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/'  # parent directiory
adir <- paste0(pdir,'analysis/')  # analysis subdirectory
ddir <- paste0(pdir,'data/')  # data subdirectory


###############################################
### QC, DIM REDUCTION, CLUSTERING, DIFF EXP ###
###############################################

# directories and clustering types to iterate through
cdir <- paste0(adir, c('02_cluster/','04_cluster/'))
edir <- paste0(adir, c('03_diff_expr/','05_diff_expr/'))
clust_type <- c('louvain_0.7','louvain_0.8')

# iterate through each pass
for (i in 1:length(cdir)){
  
  # load Seurat object
  DATA <- readRDS(paste0(cdir[i],'seurat_object.rds'))
  
  # UMAP plot with metadata columns mapped as colors
  metadata_fields <- c('organ','mouse_nr','day_post_infection')
  umap_data <- as.data.frame(DATA@reductions$umap@cell.embeddings)
  for (mf in metadata_fields){
    png(filename=paste0(cdir[i],'umap_plots/umap_metadata_',mf,'.png'), res=300, units='mm', width=190, height=130)
    print(ggplot(umap_data, aes(x=UMAP_1, y=UMAP_2, color=DATA@meta.data[[mf]])) +
            geom_point(size=0.4, stroke=0) +
            theme_classic() +
            theme(text = element_text(size=14, color='black'),
                  axis.text = element_text(size=12, color='black'),
                  legend.text = element_text(size=12)) +
            guides(color=guide_legend('',override.aes=list(size=3))))
    invisible(dev.off())
  }

  # UMAP plot with QC feature (continuous data) mapped as color intensity
  col_scale <- c("grey85","navy")
  metadata_fields <- c('nCount_RNA','perc_mito')
  p <- FeaturePlot(object=DATA, features=metadata_fields, cols=col_scale, pt.size=0.1, reduction='umap', dims=1:2, order=T, ncol=2)
  ggsave(p, filename=paste0('umap_metadata_',paste0(metadata_fields,collapse='_'),'.png'), path=paste0(cdir[i],'umap_plots/'),
         dpi=300, units="mm", width=330, height=130, limitsize=F)
  
  # UMAP plots with comparing cell cycle stage, mapped as color intensity
  metadata_fields <- c('G1.Score','S.Score','G2M.Score')
  p <- FeaturePlot(object=DATA, features=metadata_fields, cols=col_scale, pt.size=0.1, reduction='umap', dims=1:2, order=T, ncol=3)
  ggsave(p, filename='umap_metadata_cellcycle.png', path=paste0(cdir[i],'umap_plots/'), dpi=300, units='mm', width=450, height=130, limitsize=F)
  

  # UMAP plot with clusters mapped as colors
  p <- DimPlot(DATA, dims=1:2, reduction='umap', group.by=clust_type[i], pt.size=0.1, label=T)
  ggsave(p, filename=paste0('clustering_',clust_type[i],'_umap.png'), path=paste0(cdir[i],'clustering/'), dpi=300, units="mm",
         width=190, height=150, limitsize=F)
  
  
  # determine top 3 DE genes per cluster
  DATA <- SetIdent(DATA, value=as.character(DATA@meta.data[, clust_type[i]]))
  DATA_markers <- read.csv2(paste0(edir[i],'Cluster_marker_genes.csv'))
  DATA_markers %>% group_by(cluster) %>% top_n(3, avg_logFC) -> top3
  # heatmap showing top DE genes per cluster
  png(filename=paste0(edir[i],'heatmap_genes_per_cluster_top3.png'), width=2000, height=1000, res=150)
  print(DoHeatmap(object=DATA, features=as.character(unique(top3$gene)), assay='RNA', label=F))
  invisible(dev.off())
  # violin plot showing top DE genes per cluster
  png(filename=paste0(edir[i],'violinPlot_genes_per_cluster_top3.png'), width=200*10, height=200*1.5*length(as.character(unique(top3$gene)))/3, res=150)
  print(VlnPlot(object=DATA, features=as.character(unique(top3$gene)), pt.size=0, ncol=3, assay='RNA'))
  invisible(dev.off())
}


##################################
### MAPPING OF B-CELL SUBTYPES ###
##################################

# specify directory of cell subtype prediction results
sdir <- paste0(cdir[length(cdir)], 'cell_type_prediction/')

# load latest Seurat object (if not already loaded)
if (!exists('DATA')) {
  DATA <- readRDS(paste0(cdir[length(cdir)],'seurat_object.rds'))
}


# B-Cell markers
# ==============

celltype_dir <- 'bcell_types'

# load b-cell gene marker list
cell_markers <- as.list(read.csv(paste0(ddir, 'gene_lists/', celltype_dir, '.csv')))
cell_markers <- lapply(cell_markers, function(x) casefold( as.character(x[x!=""]) ) )

# match markers to those in Seurat object (i.e., correct capitalization)
cell_markers <- lapply(cell_markers, function(x) rownames(DATA)[casefold(rownames(DATA)) %in% casefold(x)])
cell_markers[lapply(cell_markers, length) == 0] <- NULL  # remove empty list entries

# subset DATA to speed up calculation times
temp <- subset(DATA, features=unique(unlist(cell_markers)))

# calculate the sum of counts for all marker genes associated with each cell type
counts <- as.matrix(temp@assays$RNA@counts)
marker_counts <- NULL
for (i in 1:length(cell_markers)) {
  if (length(cell_markers[[i]]) > 1) {
    mcs <- colSums(counts[rownames(temp) %in% cell_markers[[i]], ])
  } else {
    mcs <- counts[rownames(temp) %in% cell_markers[[i]], ]
  }
  marker_counts <- cbind(marker_counts, mcs)
}
colnames(marker_counts) <- names(cell_markers)

# add cell marker counts to Seurat object
temp <- AddMetaData(temp, as.data.frame(marker_counts))

# generate UMAP plots with mapped count values for each cell type
col_scale <- c("grey85","navy")
p <- FeaturePlot(object=temp, features=names(cell_markers), cols=col_scale, pt.size=0.05, reduction='umap',
                 dims=1:2, order=T, ncol=3)
ggsave(p, filename='UMAP_cell_type_counts.png', path=paste0(sdir, celltype_dir), dpi=300, units='mm', width=400, height=200, limitsize=F)
rm(temp)


# import b-cell correlation results
if (!file.exists(paste0(sdir, celltype_dir, '/cell_pred_correlation_', celltype_dir, '.rds'))) {
  # csv import/export is extremely slow, so save as RDS to avoid this in the future
  cors <- read.csv2(paste0(sdir, celltype_dir, '/cell_pred_correlation_', celltype_dir, '.csv'), row.names=1)
  saveRDS(cors, file=paste0(sdir, celltype_dir, '/cell_pred_correlation_', celltype_dir, '.rds'))
} else {
  cors <- readRDS(paste0(sdir, celltype_dir, '/cell_pred_correlation_', celltype_dir, '.rds'))
}
celltypes <- rownames(cors)

# add correlation data to Seurat object
temp <- AddMetaData(DATA, as.data.frame(t(cors)))

# generate UMAP plots with mapped correlation values for each cell type
col_scale <- c("grey85","navy")
p <- FeaturePlot(object=temp, features=celltypes, cols=col_scale, pt.size=0.05, reduction='umap',
                 dims=1:2, order=T, ncol=3, min.cutoff=0, max.cutoff=1)
ggsave(p, filename='UMAP_cell_type_correlation.png', path=paste0(sdir, celltype_dir), dpi=300, units='mm', width=400, height=200, limitsize=F)
rm(temp)


# B-Cell GC subtype markers
# =========================

celltype_dir <- 'bcell_types_germsub_zonesub'

# subset DATA to only contain clusters associated with the germinal center (GC)
gc_clusters <- c('louvain_0.8','3','4','5','6','7','9','14')
temp <- subset(DATA, cells = rownames(DATA@meta.data)[ (DATA@meta.data[[gc_clusters[1]]] %in% gc_clusters[2:length(gc_clusters)]) ])

# load gene marker list
cell_markers <- as.list(read.csv(paste0(ddir, 'gene_lists/', celltype_dir, '.csv')))
cell_markers <- cell_markers[startsWith(names(cell_markers), 'GC_')]
cell_markers <- lapply(cell_markers, function(x) casefold( as.character(x[x!=""]) ) )

# match markers to those in Seurat object (i.e., correct capitalization)
cell_markers <- lapply(cell_markers, function(x) rownames(temp)[casefold(rownames(DATA)) %in% casefold(x)])
cell_markers[lapply(cell_markers, length) == 0] <- NULL  # remove empty list entries

# subset data to speed up calculation times
temp <- subset(temp, features=unique(unlist(cell_markers)))

# calculate the sum of counts for all marker genes associated with each cell type
counts <- as.matrix(temp@assays$RNA@counts)
marker_counts <- NULL
for (i in 1:length(cell_markers)) {
  if (length(cell_markers[[i]]) > 1) {
    mcs <- colSums(counts[rownames(temp) %in% cell_markers[[i]], ])
  } else {
    mcs <- counts[rownames(temp) %in% cell_markers[[i]], ]
  }
  marker_counts <- cbind(marker_counts, mcs)
}
colnames(marker_counts) <- names(cell_markers)

# add cell marker counts to Seurat object
temp <- AddMetaData(temp, as.data.frame(marker_counts))

# generate UMAP plots with mapped count values for each cell type
col_scale <- c("grey85","navy")
p <- FeaturePlot(object=temp, features=names(cell_markers), cols=col_scale, pt.size=0.05, reduction='umap',
                 dims=1:2, order=T, ncol=2)
ggsave(p, filename='UMAP_cell_type_counts.png', path=paste0(sdir, celltype_dir), dpi=300, units='mm', width=400, height=200, limitsize=F)
rm(temp)


# import b-cell subtype correlation results, with germinal center subtyping
temp <- subset(DATA, cells = rownames(DATA@meta.data)[ (DATA@meta.data[[gc_clusters[1]]] %in% gc_clusters[2:length(gc_clusters)]) ])
if (!file.exists(paste0(sdir, celltype_dir, '/cell_pred_correlation_', celltype_dir, '.rds'))) {
  # csv import/export is extremely slow, so save as RDS to avoid this in the future
  cors <- read.csv2(paste0(sdir, celltype_dir, '/cell_pred_correlation_', celltype_dir, '.csv'), row.names=1)
  saveRDS(cors, file=paste0(sdir, celltype_dir, '/cell_pred_correlation_', celltype_dir, '.rds'))
} else {
  cors <- readRDS(paste0(sdir, celltype_dir, '/cell_pred_correlation_', celltype_dir, '.rds'))
}
celltypes <- rownames(cors)
celltypes <- celltypes[startsWith(celltypes, 'GC_')]
temp <- AddMetaData(temp, as.data.frame(t(cors)))

# generate UMAP plots with mapped correlation values for each cell type
col_scale <- c("grey85","navy")
p <- FeaturePlot(object=temp, features=celltypes, cols=col_scale, pt.size=0.05, reduction='umap',
                 dims=1:2, order=T, ncol=2, min.cutoff=0, max.cutoff=1)
ggsave(p, filename='UMAP_GC_cell_type_correlation.png', path=paste0(sdir, celltype_dir), dpi=300, units='mm', width=267, height=200, limitsize=F)
rm(temp)



############################
### TRAJECTORY INFERENCE ###
############################

# specify directories and clustering type(s)
cdir <- paste0(adir, '04_cluster/')  # clustering analysis direcotry
tdirs <- dir(adir, '^trajectory_\\d')  # trajectory analysis directories
clustering_use <- 'louvain_0.95'  # cluster grouping to show on plots

# create trajectory_summary directory if it does not exist
summary_dir <- paste0(adir, 'trajectory_summary/')
if (!dir.exists(summary_dir)) {dir.create(summary_dir)}

# initialize summary table
table_cols <- c('ID', 'Clustering', 'Clusters', 'Start', 'End')
traj_summary_tab <- data.frame(matrix(nrow=length(tdirs), ncol=length(table_cols))) %>% setNames(table_cols)
traj_summary_tab$ID <- tdirs
traj_summary_tab$Clustering <- clustering_use  # same clustering for all
                             
# load Seurat object
DATA <- readRDS(paste0(cdir, 'seurat_object.rds'))
       
# loop through each trajectory result subdirectory
for (i in seq(length(tdirs))) {
  tdir <- paste0(adir, tdirs[i], '/')
  
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
    keep_curve <- which(lineages@slingParams$end.given)
    
    pal <- viridis(100, end=0.95)
    dimred <- reducedDims(curves)
    colors <- pal[cut(pseudo_time[,keep_curve], breaks=100)]
    tmp <- curves; tmp@curves <- tmp@curves[keep_curve]
    
    png(filename=paste0(tdir, 'trajectory_curve', keep_curve, '_pseudotime.png'), res=300, units='mm', width=120, height=120)
    par(mgp=c(0.5, 0, 0), mar=c(2,2,2,2), cex=1.3)
    plot(dimred, xlim=range(dimred[,1]), ylim=range(dimred[,2]), col='grey85', pch=16, cex=0.3, yaxt='n', xaxt='n')
    par(new=T)
    plot(dimred, xlim=range(dimred[,1]), ylim=range(dimred[,2]), col=colors, pch=16, cex=0.3, main='Trajectory Pseudotime', yaxt='n', xaxt='n')
    lines(tmp, lwd=3, col='black', type='curves')
    box(lwd=1.5)
    invisible(dev.off())
  }
}

# write trajectory summary table to csv
write.csv(traj_summary_tab, file=paste0(summary_dir, 'trajectory_runs_summary.csv'))
          


