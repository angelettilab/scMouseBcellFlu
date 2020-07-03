

library(Seurat)
library(viridisLite)
library(ggplot2)
library(dplyr)

analysis_dir <- '/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/analysis/immcantation/mutation'


##########################
### LOAD SEURAT OBJECT ###
##########################
DATA <- readRDS('/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/analysis/06_cluster_scale/seurat_object.rds')


########################################
### LOAD AND INTEGRATE MUTATION DATA ###
########################################
db <- read.delim('/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/analysis/immcantation/mutation/VDJseq_mutation_quant.tab',
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
cols <- colnames(db)
# keep_cols <- c('V_CALL','D_CALL','J_CALL','C_CALL','V_CALL_10X','D_CALL_10X','J_CALL_10X',
#                'V_CALL_GENOTYPED','CLONE','GERMLINE_V_CALL','GERMLINE_D_CALL','GERMLINE_J_CALL',
#                cols[startsWith(cols, 'MU_')])
keep_cols <- setdiff(cols, c('MOUSE_NR', 'CELL', 'CONSCOUNT', 'UMICOUNT', 'ORGAN', 'DAY_POST_INFECTION'))
db_add <- db[, keep_cols]
rownames(db_add) <- cell_ids

# add mutation frequency data to Seurat object as metadata
DATA <- AddMetaData(DATA, db_add, col.name=colnames(db_add))

# export Seurat object with added metadata
saveRDS(DATA, file='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/analysis/immcantation/seurat_object_VDJannot.rds')

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
      if (length(featName) > 2) {
        featName <- 'Total mutations'
      } else {
        featName <- paste(featName, collapse=' + ')
      }
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
    plot(umap_coords, xlim=range(umap_coords[,1]), ylim=range(umap_coords[,2]), col='grey85', pch=16, cex=0.3, yaxt='n', xaxt='n')
    par(new=T)
    plot(plot_data[,c(1:2)], xlim=range(umap_coords[,1]), ylim=range(umap_coords[,2]), col=colors, pch=16, cex=0.3, main=featName[i], yaxt='n', xaxt='n')
  }
}


# plot mutation data on UMAP
# col_scale <- c("grey85","navy")
col_scale <- viridis(100, direction=-1)
# col_scale <- magma(100, direction=-1)

featName <- colnames(DATA@meta.data)[grepl('MU_', colnames(DATA@meta.data))]

# FeaturePlot(DATA, reduction='umap', features='MU_COUNT_CDR1_R', cols=col_scale, pt.size=0.05, order=T)
png(paste0(analysis_dir, '/umap_VDJmut_all.png'), res=300, units='mm', width=120, height=100)
plotFeat(DATA, featName=featName, featMax=20, combineMethod='sum', colorPalette=col_scale)
invisible(dev.off())
# plotFeat(DATA, featName=featName, combineMethod='sep', colorPalette=col_scale)
















