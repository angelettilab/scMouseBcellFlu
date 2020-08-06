# This script generates some additional plots from the scRNA-seq analysis pipeline results.
# Use conda environment: Sauron.v1

# load packages
suppressMessages(suppressWarnings({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(viridis)
}))

# specify relevant directories
pdir <- '/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/scMouseBcellFlu/'  # project directiory
adir <- paste0(pdir,'analysis/')  # analysis subdirectory
ddir <- paste0(pdir,'data/')  # data subdirectory


###############################################
### QC, DIM REDUCTION, CLUSTERING, DIFF EXP ###
###############################################

# directories and clustering types to iterate through
cdir <- paste0(adir, c('02_cluster/','06_cluster/'))
edir <- paste0(adir, c('03_diff_expr/','07_diff_expr/'))
clust_type <- c('louvain_0.65','HC_16')
gc_clusters <- c(tail(clust_type,n=1),'3','8','9','10','11','12','13','14','15','16')

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
sdir <- paste0(tail(cdir,1), 'cell_type_prediction/')

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
invisible(gc())


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
invisible(gc())



# B-Cell GC subtype markers
# =========================

celltype_dir <- 'bcell_types_germsub_zonesub'

# subset DATA to only contain clusters associated with the germinal center (GC)
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
invisible(gc())






