#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(optparse)))

##################################
### DEFINE PATH TO LOCAL FILES ###
##################################
cat("\nRunning VDJ INTEGRATION WITH RNA-SEQ DATA OBJECT with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),  type = "character",   metavar="character",   default='none',  help="Path to the Seurat object."),
  make_option(c("-c", "--changeo_db_path"),  type = "character",   metavar="character",   default='none',  help="Path to the ChangeO database file(s) (.tab) containing VDJ clonotype and mutation data."),
  make_option(c("-o", "--output_path"),         type = "character",   metavar="character",   default='none',  help="Output directory with optional filename (ending in '.rds'). If a filename is not provided, the output Seurat object will be named the same as the input object, but appended with _VDJannot.rds. If the path is not specified, the output directory will be the same as the directory containing the input Seurat object.")
)
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))
#---------

# opt <- list(Seurat_object_path='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/scMouseBcellFlu/analysis/06_cluster/seurat_object.rds',
#             changeo_db_path='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/scMouseBcellFlu/analysis/immcantation/mutation2',
#             output_path='/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/scMouseBcellFlu/analysis/immcantation/seurat_object_VDJannot.rds')


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


#############################################
### SPECIFY MUTATION DATA FILENAMES/PATHS ###
#############################################
mut_files <- list.files(opt$changeo_db_path, pattern='mutation_quant.tab')
heavy_file <- mut_files[startsWith(casefold(mut_files), 'igh')]
light_file <- mut_files[startsWith(casefold(mut_files), 'igkl')]

if (length(heavy_file) == 0) {
  stop(paste('No heavy chain Change-O DB file found in changeo_db_path.\n',
             'The path must contain a seq db file named "IGH_mutation_quant.tab".'))
}
#---------


################################
### PROCESS HEAVY CHAIN DATA ###
################################
cat('Loading heavy chain Change-O database file ...\n')
db <- read.delim(file.path(opt$changeo_db_path, heavy_file), stringsAsFactors=F)

# get cell IDs in the same format as used in Seurat object
cell_ids_heavy <- unlist(lapply(db$SEQUENCE_ID, function(x) paste(unlist(strsplit(x, '-'))[-2], collapse='_')))

# remove duplicated cell IDs in ChangeO database (should not be any duplicated cell IDs for heavy chain data)
dup_ids <- cell_ids_heavy[duplicated(cell_ids_heavy)]
if (length(dup_ids) > 0) {
  warning('Duplicate cell IDs found in heavy chain data. Duplicate entries will be removed.')
  db <- db[!duplicated(cell_ids_heavy), ]
  cell_ids_heavy <- cell_ids_heavy[!duplicated(cell_ids_heavy)]
}

# specify which columns to add to metadata
cols <- colnames(db)
keep_cols <- setdiff(cols, c('MOUSE_NR', 'CELL', 'CONSCOUNT', 'UMICOUNT', 'ORGAN', 'DAY_POST_INFECTION'))
db_heavy <- db[, keep_cols]
rownames(db_heavy) <- cell_ids_heavy

# prefix all heavy chain data columns with HC (Heavy Chain)
colnames(db_heavy) <- paste0('HC_', colnames(db_heavy))
#---------



################################
### PROCESS LIGHT CHAIN DATA ###
################################
# specify which columns of the light-chain data to keep
light_cols <- c('V_CALL', 'J_CALL', 'CDR3_IMGT', 'MU_COUNT_TOT', 'MU_FREQ_TOT')

if (length(light_file) == 0) {
  cat('No light chain Change-O database file found. Only heavy chain data will be integrated with Seurat object.\n')
  db_paired <- db_heavy
} else {
  cat('Loading light chain Change-O database file ...\n')
  db <- read.delim(file.path(opt$changeo_db_path, light_file), stringsAsFactors=F)
  
  # get cell IDs in the same format as used in Seurat object
  cell_ids_light <- unlist(lapply(db$SEQUENCE_ID, function(x) paste(unlist(strsplit(x, '-'))[-2], collapse='_')))
  
  # filter columns to only keep a specific set
  missing_cols <- setdiff(light_cols, colnames(db))
  if (length(missing_cols) == length(light_cols)) {
    cat('None of the requested columns were found in the light chain Change-O DB file:\n')
    cat(missing_cols, '\n\n')
    stop('Invalid light chain dataset.')
  } else if (length(missing_cols) > 0) {
    cat('WARNING! The following column(s) were not found in the light chain Change-O DB file:\n')
    cat(missing_cols, '\n\n')
  }
  db <- db[, colnames(db) %in% light_cols]
  
  # filter out light chains that have no matching heavy chain
  keep <- cell_ids_light %in% cell_ids_heavy
  db <- db[keep, ]
  cell_ids_light <- cell_ids_light[keep]
  
  # pair light chains with heavy chains based on cell IDs
  cell_index <- match(cell_ids_light, cell_ids_heavy)
  cell_indices <- unique(cell_index)
  
  # set a max number of paired light chains, to avoid creating too many columns
  max_pairs <- max(tabulate(cell_index))
  if (max_pairs > 5) {
    cat('WARNING! Found instances of > 5 light chains paired with a heavy chain.')
    cat('Only the first 5 paired light chains will be retained!\n\n')
    max_pairs <- 5
  }
  
  # initialize db_light data frame; label columns as "LC1_V_CALL", "LC1_J_CALL" ... "LC2_V_CALL", "LC2_J_CALL" ... etc.
  merged_light_colnames <- paste0('LC', unlist(lapply(1:max_pairs, function(x) rep(x, ncol(db)))), '_', colnames(db))
  db_light <- as.data.frame(matrix(NA, nrow=length(cell_ids_heavy), ncol=length(merged_light_colnames),
                                   dimnames=list(cell_ids_heavy, merged_light_colnames)))
  
  # populate db_light dataframe with paired light chain sequences
  cat('Pairing light and heavy chains ...\n')
  for (i in seq(length(cell_indices))) {
    db_rows <- db[cell_index == cell_indices[i], ]
    for (j in seq(nrow(db_rows))) {
      db_light[cell_indices[i], seq(ncol(db)) + ncol(db)*(j-1)] <- db_rows[j, ]
    }
  }
  
  # merge heavy and light chain db dataframes
  if (any(rownames(db_heavy) != rownames(db_light))) {
    # just a check, shouldn't get this error...
    stop('Heavy and light chain dataframes are not aligned!')
  }
  db_paired <- cbind(db_heavy, db_light)
}
#---------


############################################
### ADD DATA TO SEURAT OBJECT AND EXPORT ###
############################################
cat('Adding VDJ data to Seurat object metadata ...\n')
DATA <- AddMetaData(DATA, db_paired, col.name=colnames(db_paired))

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
png(file.path(opt$changeo_db_path, 'umap_VDJmutfreq_hc.png'), res=300, units='mm', width=120, height=100)
plotFeat(DATA, featName='HC_MU_FREQ_TOT', featMax=0.03, colorPalette=col_scale)
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











