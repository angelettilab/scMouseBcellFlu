
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
clust_type <- c('louvain_0.7','louvain_0.5')
filtered_clusters <- c(0, 11, 12, 13, 15, 17)

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
  p <- DimPlot(DATA, dims=1:2, reduction='umap', group.by=clust_type[i], pt.size=0.1)
  ggsave(p, filename=paste0('clustering_',clust_type[i],'_umap.png'), path=paste0(cdir[i],'clustering/'), dpi=300, units="mm",
         width=190, height=150, limitsize=F)
  
  # UMAP plot showing filtered clusters
  if (i == 1){
    DATA@meta.data$filter_status <- factor(DATA@meta.data[[clust_type[i]]] %in% filtered_clusters, labels=c('Keep','Remove'))
    p <- DimPlot(object=DATA, group.by='filter_status', cols=c('lightgray','firebrick'), pt.size=0.1, reduction='umap')
    ggsave(p, filename='removed_clusters_umap.png', path=paste0(cdir[i],'cell_type_prediction/main_cell_types/'), dpi=300, units="mm",
           width=190, height=150, limitsize=F)
  }

  
  # determine top 3 DE genes per cluster
  DATA <- SetIdent(DATA, value=as.character(DATA@meta.data[, clust_type[i]]))
  DATA_markers <- read.csv2(paste0(edir[i],'/Cluster_marker_genes.csv'))
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

# import b-cell correlation results
celltype_dir <- 'bcell_types'
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


# import b-cell subtype correlation results, with germinal center subtyping
celltype_dir <- 'bcell_types_germsub_zonesub'
if (!file.exists(paste0(sdir, celltype_dir, '/cell_pred_correlation_', celltype_dir, '.rds'))) {
  # csv import/export is extremely slow, so save as RDS to avoid this in the future
  cors <- read.csv2(paste0(sdir, celltype_dir, '/cell_pred_correlation_', celltype_dir, '.csv'), row.names=1)
  saveRDS(cors, file=paste0(sdir, celltype_dir, '/cell_pred_correlation_', celltype_dir, '.rds'))
} else {
  cors <- readRDS(paste0(sdir, celltype_dir, '/cell_pred_correlation_', celltype_dir, '.rds'))
}
celltypes <- rownames(cors)
celltypes <- celltypes[startsWith(celltypes, 'GC_')]

# subset DATA to only contain clusters associated with the germinal center (GC) and add correlation values
gc_clusters <- c('louvain_0.5','1','2','4','5','10','11')
# temp <- subset(DATA, cells = rownames(DATA@meta.data)[ (DATA@meta.data[[gc_clusters[1]]] %in% gc_clusters[2:length(gc_clusters)]) ]
temp <- subset(DATA, cells = rownames(DATA@meta.data)[ (DATA@reductions$umap@cell.embeddings[,2] < 10) & 
                                                         (DATA@reductions$umap@cell.embeddings[,2] > -4.6667*DATA@reductions$umap@cell.embeddings[,1] - 19)])
temp <- AddMetaData(temp, as.data.frame(t(cors)))

# generate UMAP plots with mapped correlation values for each cell type
col_scale <- c("grey85","navy")
p <- FeaturePlot(object=temp, features=celltypes, cols=col_scale, pt.size=0.05, reduction='umap',
                 dims=1:2, order=T, ncol=2, min.cutoff=0, max.cutoff=1)
ggsave(p, filename='UMAP_GC_cell_type_correlation.png', path=paste0(sdir, celltype_dir), dpi=300, units='mm', width=267, height=200, limitsize=F)
rm(temp)




