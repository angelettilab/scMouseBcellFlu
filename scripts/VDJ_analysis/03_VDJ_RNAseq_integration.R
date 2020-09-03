#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(optparse)))

##################################
### DEFINE PATH TO LOCAL FILES ###
##################################
cat("\nRunning VDJ INTEGRATION WITH RNA-SEQ DATA OBJECT with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),  type = "character",   metavar="character",   default='none',  help="Path to the Seurat object."),
  make_option(c("-c", "--changeo_db_path"),     type = "character",   metavar="character",   default='none',  help="Path to the ChangeO database file (.tab) containing VDJ clonotype and mutation data."),
  make_option(c("-o", "--output_path"),         type = "character",   metavar="character",   default='none',  help="Output directory with optional filename (ending in '.rds'). If a filename is not provided, the output Seurat object will be named the same as the input object, but appended with _VDJannot.rds. If the path is not specified, the output directory will be the same as the directory containing the input Seurat object.")
)
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))
#---------


##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading libraries ...\n")
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(viridisLite)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(dplyr)))
#---------


##########################
### LOAD SEURAT OBJECT ###
##########################
cat('\n\nLoading Seurat object ...\n')
DATA <- readRDS(opt$Seurat_object_path)
#---------


########################################
### LOAD AND INTEGRATE MUTATION DATA ###
########################################
# load ChangeO database file
cat('Loading ChangeO database file ...\n')
db <- read.delim(opt$changeo_db_path, stringsAsFactors=F)

# get cell IDs in the same format as used in Seurat object
cell_ids <- unlist(lapply(db$SEQUENCE_ID, function(x) paste(unlist(strsplit(x, '-'))[-2], collapse='_')))

# remove duplicated cell IDs in ChangeO database (should not be any duplicates)
dup_ids <- cell_ids[duplicated(cell_ids)]
db <- db[!duplicated(cell_ids), ]
cell_ids <- cell_ids[!duplicated(cell_ids)]

# specify which columns to add to metadata
cols <- colnames(db)
keep_cols <- setdiff(cols, c('MOUSE_NR', 'CELL', 'CONSCOUNT', 'UMICOUNT', 'ORGAN', 'DAY_POST_INFECTION'))
db_add <- db[, keep_cols]
rownames(db_add) <- cell_ids

# add mutation frequency data to Seurat object as metadata
cat('Adding VDJ data to Seurat object metadata ...\n')
DATA <- AddMetaData(DATA, db_add, col.name=colnames(db_add))

# parse output file and path specifications
if (endsWith(casefold(opt$output_path), '.rds')) {
  out_name <- basename(opt$output_path)
  out_dir <- dirname(opt$output_path)
} else {
  out_name <- sub('.rds', '_VDJannot.rds', casefold(basename(opt$Seurat_object_path)))
  if (opt$output_path == 'none') {
    out_dir <- dirname(opt$Seurat_object_path)
  } else {
    out_dir <- opt$output_path
  }
}
out_full <- file.path(out_dir, out_name)

if (!dir.exists(out_dir)) { dir.create(out_dir, recursive=T) }
setwd(out_dir)
#---------


###################################
### EXPORT MERGED SEURAT OBJECT ###
###################################
cat('Writing Seurat object containing VDJ data:', out_name, '...\n')
saveRDS(DATA, file=out_full)
#---------


####################################
### PLOT UMAP WITH MUTATION FREQ ###
####################################

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
col_scale <- viridis(100, direction=-1)
png(file.path(dirname(opt$changeo_db_path), 'umap_VDJmutfreq_all.png'), res=300, units='mm', width=120, height=100)
plotFeat(DATA, featName='MU_FREQ_TOT', featMax=0.03, colorPalette=col_scale)
invisible(dev.off())
#---------


##########################
### PRINT SESSION INFO ###
##########################
cat('R SESSION INFO:\n')
Sys.info()
cat('\n\n\n\n')
sessionInfo()
#---------











