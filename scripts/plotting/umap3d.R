
# 3D UMAP
library(Seurat)
library(plotly)
library(htmlwidgets)

# function for generating color hues (mimics default ggplot palette)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# load Seurat object
DATA <- readRDS('/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/scMouseBcellFlu/analysis/06_cluster/seurat_object.rds')

# run UMAP
DATA <- RunUMAP(object = DATA, assay="RNA", reduction="mnn", dims=1:50, n.components=3, n.neighbors=50, spread=5,
                repulsion.strength=0.5, min.dist=0.1, verbose=F, num_threads=0, n.epochs=500, metric = "euclidean",
                seed.use=42, learning.rate=0.5, negative.sample.rate=7, reduction.name="umap3", reduction.key="umap3_")

# extract 3D UMAP embedding into dataframe
df <- data.frame(DATA@reductions$umap3@cell.embeddings)
# add_data <- DATA$PlasmaCells
# add_data[add_data < 0] <- 0
df <- data.frame(df, cluster_use=as.factor(DATA@meta.data$HC_16))

pal <- gg_color_hue(length(unique(df$cluster_use)))
# pal <- c(scales::hue_pal()(8), RColorBrewer::brewer.pal(9,"Set1"), RColorBrewer::brewer.pal(8,"Set2") )
# pal <- viridisLite::magma(100)
p_State <- plot_ly(df, x=~umap3_1, y=~umap3_2, z=~umap3_3,color=~cluster_use,
                   colors=pal, alpha=1, alpha_stroke=0, size=0.1) %>%  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'UMAP_1'), yaxis = list(title = 'UMAP_2'),zaxis = list(title = 'UMAP_3')))
ggplotly(p_State, height=800, width=800)
saveWidget(p_State, "/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/scMouseBcellFlu/analysis/06_cluster/umap_plots/3d_umap.html")

