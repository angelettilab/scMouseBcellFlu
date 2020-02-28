
library(Seurat)
library(ggplot2)
library(dplyr)
# source('/Users/jonrob/Documents/NBIS/repos/niceRplots/R/plotting_functions.R')

# # color palette function
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }

# specify relevant directories
pdir <- '/Users/jonrob/Documents/NBIS/LTS_projects/2020_DAngeletti/'  # parent directiory
adir <- paste0(pdir,'analysis/')  # analysis subdirectory
ddir <- paste0(pdir,'data/')  # data subdirectory

# directories and clustering types to iterate through
cdir <- paste0(adir, c('02_cluster/','04_cluster/','06_cluster/'))
edir <- paste0(adir, c('03_diff_expr/','05_diff_expr/','07_diff_expr/'))
clust_type <- c('HC_12','HC_15','HC_11')

# iterate through each pass
for (i in 1:3){
  # load Seurat object
  DATA <- readRDS(paste0(cdir[i],'seurat_object.rds'))
  
  # UMAP plot with metadata columns mapped as colors
  metadata_fields <- c('organ','infection','day_post_infection')
  umap_data <- as.data.frame(DATA@reductions$umap@cell.embeddings)
  for (mf in metadata_fields){
    png(filename=paste0(cdir[i],'umap_plots/umap_metadata_',mf,'.png'), res=300, units='mm', width=170, height=130)
    print(ggplot(umap_data, aes(x=UMAP_1, y=UMAP_2, color=DATA@meta.data[[mf]])) +
            geom_point(size=0.3, stroke=0) +
            theme_classic() +
            theme(text = element_text(size=14, color='black'),
                  axis.text = element_text(size=12, color='black'),
                  legend.text = element_text(size=12)) +
            guides(color=guide_legend('',override.aes=list(size=3))))
    invisible(dev.off())
  }

  # UMAP plot with QC feature (continuous data) mapped as color intensity
  metadata_fields <- c('nCount_RNA','nFeature_RNA','S.Score')
  col_scale <- c("grey85","navy")
  umap_data <- as.data.frame(DATA@reductions$umap@cell.embeddings)
  for (mf in metadata_fields){
    p <- FeaturePlot(object=DATA, features=mf, cols=col_scale, pt.size=0.1, reduction='umap', dims=1:2, order=T)
    ggsave(p, filename=paste0('umap_metadata_',gsub('[.]','_',mf),'.png'), path=paste0(cdir[i],'umap_plots/'), dpi=300, units="mm", width=170, height=130, limitsize=F)
  }

  # UMAP plot with HC clusters mapped as colors
  p <- DimPlot(DATA, dims=1:2, reduction='umap', group.by=clust_type[i], pt.size=0.1)
  ggsave(p, filename=paste0('clustering_',clust_type[i],'_umap.png'), path=paste0(cdir[i],'clustering/'), dpi=300, units="mm",
         width=170, height=150, limitsize=F)

  
  # determine top 3 DE genes per cluster
  DATA <- SetIdent(DATA, value=as.character(DATA@meta.data[, clust_type[i]]))
  DATA_markers <- read.csv2(paste0(edir[i],'/Cluster_marker_genes.csv'))
  DATA_markers %>% group_by(cluster) %>% top_n(3, avg_logFC) -> top3
  # heatmap showing top DE genes per cluster
  png(filename=paste0(edir[i],'heatmap_genes_per_cluster_top3.png'), width=2000, height=800, res=150)
  print(DoHeatmap(object=DATA, features=as.character(unique(top3$gene)), assay='RNA', slot='data'))
  invisible(dev.off())
  # violin plot showing top DE genes per cluster
  png(filename=paste0(edir[i],'violinPlot_genes_per_cluster_top3.png'), width=200*10, height=200*1.5*length(as.character(unique(top3$gene)))/3, res=150)
  print(VlnPlot(object=DATA, features=as.character(unique(top3$gene)), pt.size=0, ncol=3, assay='RNA'))
  invisible(dev.off())
}







