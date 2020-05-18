
# 3D UMAP
library(Seurat)
library(plotly)

DATA <- RunUMAP(object = DATA, assay="RNA", reduction="pca", dims=1:50, n.components=3, n.neighbors=50, spread=3,
                repulsion.strength=1, min.dist=0.01, verbose=T, num_threads=0, n.epochs=500, metric = "correlation",
                seed.use=42, learning.rate=0.5, negative.sample.rate=5, reduction.name="umap_3d", reduction.key="umap_")

df <- data.frame(DATA@reductions$umap_3d@cell.embeddings)
df <- data.frame(df, cluster_use=DATA$cell_pred_correlation_bcell_types_germsub)

pal <- c(scales::hue_pal()(8), RColorBrewer::brewer.pal(9,"Set1"), RColorBrewer::brewer.pal(8,"Set2") )
p_State <- plot_ly(df, x=~umap_1, y=~umap_2, z=~umap_3,color=~cluster_use,
                   colors=pal, alpha=1, alpha_stroke=0, size=0.1) %>%  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'UMAP_1'), yaxis = list(title = 'UMAP_2'),zaxis = list(title = 'UMAP_3')))
ggplotly(p_State, height=800, width=800)
htmlwidgets::saveWidget(p_State, "/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/analysis/05_cluster_scale/umap_plots/3d_umap.html")


