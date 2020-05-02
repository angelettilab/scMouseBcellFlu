

library(Seurat)
library(viridisLite)
library(ggplot2)
library(dplyr)


##########################
### LOAD SEURAT OBJECT ###
##########################
DATA <- readRDS('/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/analysis/04_cluster/seurat_object.rds')


########################################
### LOAD AND INTEGRATE MUTATION DATA ###
########################################
db <- read.delim('/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/analysis/immcantation/mutation/VDJseq_mutation_freq.tab',
                 stringsAsFactors=F)

# get cell IDs in the same format as used in Seurat object
cell_ids <- unlist(lapply(db$SEQUENCE_ID, function(x) paste(unlist(strsplit(x, '-'))[-2], collapse='_')))

# need to deal with duplicated cell IDs (TODO: improve this part!)
dup_ids <- cell_ids[duplicated(cell_ids)]
# for (id in dup_ids) {
#   indx <- which(cell_ids %in% id)
#   db$MU_FREQ[indx[1]] <- mean(db$MU_FREQ[indx])  # average mutation frequencies
# }
db <- db[!duplicated(cell_ids), ]
cell_ids <- cell_ids[!duplicated(cell_ids)]

# specify which columns to add to metadata
db_add <- db
rownames(db_add) <- cell_ids

# add mutation frequency data to Seurat object as metadata
DATA <- AddMetaData(DATA, db_add, col.name=colnames(db_add))


# define custom plotting function
plotFeat <- function(SeurObj, featName, featMax=Inf, combineMethod='sum', colorPalette=viridis(100)){
  # SeurObj: Seurat Object
  # featName: Column name of metadata to color cells by (NA values will be light gray)
  # featMax: Max value above which feature values will be trimmed
  # combineMethod: If multiple feature names are provided in featName,
  #                'sum': sum the feature values
  #                'sep': plot features in separate plots
  # colorPalette: color palette used to color cells
  
  umap_coords <- SeurObj@reductions$umap@cell.embeddings
  featData <- SeurObj@meta.data[featName]
  nPlots <- 1
  if (length(featName) > 1) {
    if (combineMethod == 'sum') {
      featName <- paste(featName, collapse=' + ')
      featData <- as.data.frame(rowSums(featData)) %>% setNames('Value')
    } else if (combineMethod == 'sep') {
      nPlots <- length(featName)
      if (length(featMax) == 1) {
        featMax <- rep(featMax, nPlots)
      } else if (length(featMax) != nPlots) {
        stop('Number of elements in featMax must be one, or equal to the number of elements in featName.')
      }
    } else {stop('Invalid combineMethod!')}
  }
  
  n_plot_cols <- min(3, nPlots)
  n_plot_rows <- ceiling(nPlots/n_plot_cols)
  par(mfrow=c(n_plot_rows, n_plot_cols), mgp=c(0.5, 0, 0), mar=c(2,2,2,2))
  
  for (i in seq(nPlots)) {
    plot_data <- as.data.frame(cbind(umap_coords, featData[, i]))
    plot_data <- plot_data[order(plot_data[, 3], na.last=NA), ]
    plot_data[plot_data[,3] > featMax[i], 3] <- featMax[i]
    
    colors <- colorPalette[cut(plot_data[,3], breaks=length(colorPalette))]
    plot(umap_coords, xlim=range(umap_coords[,1]), ylim=range(umap_coords[,2]), col='grey85', pch=16, cex=0.2, yaxt='n', xaxt='n')
    par(new=T)
    plot(plot_data[,c(1:2)], xlim=range(umap_coords[,1]), ylim=range(umap_coords[,2]), col=colors, pch=16, cex=0.2, main=featName[i], yaxt='n', xaxt='n')
  }
}


# plot mutation data on UMAP
# col_scale <- c("grey85","navy")
col_scale <- viridis(100, direction=-1)
col_scale <- magma(100, direction=-1)
# FeaturePlot(DATA, reduction='umap', features='MU_COUNT_CDR1_R', cols=col_scale, pt.size=0.05, order=T)
plotFeat(DATA, featName=featName, featMax=20, combineMethod='sum', colorPalette=col_scale)
# plotFeat(DATA, featName=featName, combineMethod='sep', colorPalette=col_scale)





















