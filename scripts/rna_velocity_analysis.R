
# specify relevant directories
projdir <- '/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910'
loomdir <- paste0(projdir, '/data/rna_velocity/loom_files')
velodir <- paste0(projdir, '/analysis/rna_velocity')
if (!dir.exists(velodir)) { dir.create(velodir, recursive=T) }


# load required libraries
library(velocyto.R)
library(Matrix.utils)
library(Seurat)

# function for generating color hues (mimics default ggplot palette)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# define matrix merge function (from a newer version of Seurat)
mergeSparseMatrices <- function(mat1, mat2){
  mat1.names <- rownames(mat1)
  mat2.names <- rownames(mat2)
  if (length(mat1.names) == length(mat2.names) && all(mat1.names == mat2.names)) {
    new.mat <- cbind(mat1, mat2)
  } else {
    mat1 <- as(mat1, Class='sparseMatrix')
    mat2 <- as(mat2, Class='sparseMatrix')
    all.names <- union(x = mat1.names, y = mat2.names)
    new.mat1 <- as(matrix(0, length(all.names), ncol(mat1), dimnames=list(all.names, colnames(mat1))), Class='sparseMatrix')
    new.mat2 <- as(matrix(0, length(all.names), ncol(mat2), dimnames=list(all.names, colnames(mat2))), Class='sparseMatrix')
    new.mat1[match(mat1.names, rownames(new.mat1)), ] <- mat1
    new.mat2[match(mat2.names, rownames(new.mat2)), ] <- mat2
    new.mat <- cbind(new.mat1, new.mat2)
  }
  colnames(new.mat) <- c(colnames(mat1), colnames(mat2))
  return(new.mat)
}


# iterate over sample folders (ignore folder named "unused")
samples <- setdiff(list.dirs(loomdir, full.names=F, recursive=F), 'unused')
merged_data <- NULL
for (s in samples) {
  
  # load .loom file
  loomfile <- dir(paste0(loomdir, '/', s), '*[.]loom', full.names=T)
  if (length(loomfile) == 0) {
    cat('No loom file found in', s, 'directory - skipping.\n')
    next
  }
  ldat <- read.loom.matrices(loomfile)
  
  # strip all except the sequence from barcode
  barcodes <- unlist(lapply(colnames(ldat$spliced), function(b) tail(unlist(strsplit(b,':|x')), n=1)))
  
  # append barcode with sample name to match formatting in Seurat object and assign as column name
  barcodes <- paste0(barcodes, '_', s)
  for (i in seq(length(ldat))) { colnames(ldat[[i]]) <- barcodes }
  
  # combine with merged data
  if (is.null(merged_data)) {
    merged_data <- ldat
  } else {
    for (i in seq(length(merged_data))) {
      merged_data[[i]] <- mergeSparseMatrices(merged_data[[i]], ldat[[i]])
    }
  }
  cat('Loaded loom file for sample', s, '\n')
}

# extract spliced and unspliced counts
emat <- merged_data$spliced
nmat <- merged_data$unspliced

# load Seurat object
DATA <- readRDS('/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/analysis/04_cluster/seurat_object.rds')

# subset Seurat object to only include cells in loom files
DATA <- subset(DATA, cells=colnames(DATA)[colnames(DATA) %in% colnames(emat)])

# remove data for cells not present in Seurat object
emat <- emat[, colnames(DATA)[colnames(DATA) %in% colnames(emat)]]
nmat <- nmat[, colnames(DATA)[colnames(DATA) %in% colnames(nmat)]]

# specify cluster labels and remove clusters with less than 2 cells
cluster.label <- factor(DATA@meta.data$louvain_0.8)
clust_tab <- table(cluster.label)
if (any(clust_tab < 2)) {
  rem_clusters <- names(clust_tab)[clust_tab < 2]
  keep_cells <- rownames(DATA@meta.data)[!(cluster.label %in% rem_clusters)]
  DATA <- subset(DATA, cells=keep_cells)
  emat <- emat[, keep_cells]
  nmat <- nmat[, keep_cells]
  cluster.label <- factor(DATA@meta.data$louvain_0.8)
}
names(cluster.label) <- rownames(DATA@meta.data)

# specify color palette and cell colors
color.pal <- gg_color_hue(nlevels(cluster.label))
cell.colors <- color.pal[cluster.label]
names(cell.colors) <- names(cluster.label)

# get cell-cell distances and embeddings
ndims <- 30  # TODO: make this an input parameter
# cell.dist <- as.dist(1-cor(t(DATA@reductions$pca@cell.embeddings[, 1:ndims])))
cell.dist <- as.dist(1-cor(t(DATA@reductions$mnn@cell.embeddings[, 1:ndims])))
emb <- DATA@reductions$umap@cell.embeddings

# filter genes based on the minimum average expresion magnitude (in at least one of the clusters),
# and output total number of resulting valid genes
emat <- filter.genes.by.cluster.expression(emat, cluster.label, min.max.cluster.average=0.1)
nmat <- filter.genes.by.cluster.expression(nmat, cluster.label, min.max.cluster.average=0.01)
length(intersect(rownames(emat),rownames(nmat)))

# estimate RNA velocity 
if (file.exists(paste0(velodir, '/velocity_estimates.rds'))) {
  rvel.cd <- readRDS(paste0(velodir, '/velocity_estimates.rds'))
} else {
  # using gene-relative model with k=20 cell kNN pooling and using top/bottom 2% quantiles for gamma fit
  fit.quantile <- 0.02
  rvel.cd <- gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells=20, cell.dist=cell.dist, fit.quantile=fit.quantile)
  if (any(colSums(rvel.cd$current) == 0) | any(as.logical(duplicated(as.matrix(rvel.cd$current), MARGIN=2)))) {
    # cells with all zero counts or identical count vectors will result in error, so remove
    # these cells and re-run the velocity estimates function
    keep_cells <- !(colSums(rvel.cd$current) == 0) & !as.logical(duplicated(as.matrix(rvel.cd$current), MARGIN=2))
    emb <- emb[keep_cells, ]
    emat <- emat[, keep_cells]
    nmat <- nmat[, keep_cells]
    cell.dist <- as.dist(as.matrix(cell.dist)[keep_cells, keep_cells])
    cell.colors <- cell.colors[keep_cells]
    cluster.label <- cluster.label[keep_cells]
    
    rvel.cd <- gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells=20, cell.dist=cell.dist, fit.quantile=fit.quantile)
  }
  # export velocity estimates structure to file
  saveRDS(rvel.cd, paste0(velodir, '/velocity_estimates.rds'))
}

# visualize velocity on the embedding, using correlation-based transition matrix
if (file.exists(paste0(velodir, '/velocity_correlations.rds'))) {
  vel.vis <- readRDS(paste0(velodir, '/velocity_correlations.rds'))
} else {
  vel.vis <- show.velocity.on.embedding.cor(emb, rvel.cd, n=300, scale='sqrt', cell.colors=ac(cell.colors,alpha=0.5),
                                            cex=0.8, arrow.scale=10, show.grid.flow=TRUE, min.grid.cell.mass=0.5,
                                            grid.n=40, arrow.lwd=1, do.par=F, cell.border.alpha=0.1)
  saveRDS(vel.vis, paste0(velodir, '/velocity_correlations.rds'))
}

# to replot with different parameters (faster, uses results from prev plot to speed up)
png(filename=paste0(velodir, '/umap_velocity.png'), res=300, units='mm', width=150, height=150)
show.velocity.on.embedding.cor(emb, rvel.cd, cc=vel.vis$cc, n=300, scale='sqrt', cell.colors=ac(cell.colors,alpha=0.5),
                               cex=0.8, arrow.scale=10, show.grid.flow=TRUE, min.grid.cell.mass=0.5,
                               grid.n=70, arrow.lwd=1, do.par=F, cell.border.alpha=0.1)
invisible(dev.off())


# # visualize velocity on the embedding, using euclidean-based transition matrix
# vel.vis <- show.velocity.on.embedding.eu(emb, rvel.cd, n=300, scale='sqrt', cell.colors=ac(cell.colors,alpha=0.5),
#                                           cex=0.8, arrow.scale=10, show.grid.flow=TRUE, arrow.lwd=1, do.par=F)


