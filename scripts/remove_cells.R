#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nRunning REMOVE CELLS with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object FILE."),
  make_option(c("-e", "--remove"),                type = "character",   metavar="character",   default='none',  help="Metadata column name(s) and value(s) corresponding to cells that should be removed. Input as 'column1, valA, valB; column2, valA'. For example, 'HC_12, 4, 8; louvain_0.8, 3'"),
  make_option(c("-k", "--keep"),                  type = "character",   metavar="character",   default='all',   help="Metadata column name(s) and value(s) corresponding to cells that should be kept. Can be used in combination with, or instead of, the --remove input."),
  make_option(c("-c", "--combine_method"),        type = "character",   metavar="character",   default='union', help="Options are 'union' or 'intersection'. Only relevant if specifying multiple metadata column names, and/or using both the '--remove' and '--keep' inputs. If 'union', then cells will be removed if they satisfy ANY of the removal criteria specified by the 'remove' and/or !('keep') inputs. If 'intersection', cells are removed only if they satisfy ALL removal critera."),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output DIRECTORY.")
) 
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)
#---------



##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
library(Seurat)
library(ggplot2)
#---------



#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\n### LOADING Seurat.v3 OBJECT ###\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------



#######################################
### REMOVE CELLS FROM SEURAT OBJECT ###
#######################################
cat("\n### REMOVING CELLS FROM THE DATA ###\n")

# determine which cells should be marked as "remove"
if (opt$remove != 'none') {
  remove_meta <- lapply(trimws(casefold(unlist(strsplit(opt$remove, ";")))), function(x) trimws(unlist(strsplit(x,','))))
  meta_cols <- unlist(lapply(remove_meta, function(x) x[1]))
  meta_cols <- colnames(DATA@meta.data)[match(meta_cols, casefold(colnames(DATA@meta.data)))]
  if (any(is.na(meta_cols))) {
    cat("Could not find the following columns in meta.data:\n")
    cat(unlist(lapply(remove_meta, function(x) x[1]))[is.na(meta_cols)])
    stop('Invalid meta data field.')
  }
  remove_meta <- lapply(remove_meta, function(x) x[-1])
  remove_cells <- matrix(FALSE, nrow=nrow(DATA@meta.data), ncol=length(remove_meta))
  for (i in 1:length(remove_meta)) {
    remove_cells[, i] <- casefold(DATA@meta.data[[meta_cols[i]]]) %in% casefold(remove_meta[[i]])
  }
} else {
  remove_cells <- matrix(nrow=nrow(DATA@meta.data), ncol=0)
}

# determine which cells should be marked as "keep"
if (opt$keep != 'all') {
  keep_meta <- lapply(trimws(casefold(unlist(strsplit(opt$keep, ";")))), function(x) trimws(unlist(strsplit(x,','))))
  meta_cols <- unlist(lapply(keep_meta, function(x) x[1]))
  meta_cols <- colnames(DATA@meta.data)[match(meta_cols, casefold(colnames(DATA@meta.data)))]
  if (any(is.na(meta_cols))) {
    cat("Could not find the following columns in meta.data:\n")
    cat(unlist(lapply(keep_meta, function(x) x[1]))[is.na(meta_cols)])
    stop('Invalid meta data field.')
  }
  keep_meta <- lapply(keep_meta, function(x) x[-1])
  keep_cells <- matrix(FALSE, nrow=nrow(DATA@meta.data), ncol=length(keep_meta))
  for (i in 1:length(keep_meta)) {
    keep_cells[, i] <- casefold(DATA@meta.data[[meta_cols[i]]]) %in% casefold(keep_meta[[i]])
  }
  remove_cells <- cbind(remove_cells, !keep_cells)
}

# combine results and decide which cells will be kept and removed
if (ncol(remove_cells) == 0) {
  stop('No cells were specified for removal.')
}
if (grepl('intersect', casefold(opt$combine_method))) {
  remove_cells <- rowSums(remove_cells) == ncol(remove_cells)
} else if ('union' %in% casefold(opt$combine_method)) {
  remove_cells <- rowSums(remove_cells) > 0
} else {
  stop('Invalid COMBINE_METHOD input. Valid options are "intersection" or "union".')
}
if (all(remove_cells)) {
  stop('All cells satisfied criteria for removal.')
}

# Plot UMAP showing cells to be removed
DATA@meta.data$remove_cells <- ifelse(remove_cells, 'Remove', 'Keep')
p <- DimPlot(object=DATA, pt.size=0.1, reduction='umap', group.by='remove_cells', cols=c('grey85','firebrick'))
source_dir <- sub('\\w+[.]rds$', '', opt$Seurat_object_path)
if (!dir.exists(paste0(source_dir,'umap_plots'))) { dir.create(paste0(source_dir,'umap_plots')) }
ggsave(p, filename='umap_removed_cells.png', path=paste0(source_dir,'umap_plots'), dpi=300, units="mm",
       width=190, height=150, limitsize=F)
DATA@meta.data$remove_cells <- NULL

# remove cells
DATA <- subset(DATA, cells=rownames(DATA@meta.data)[!remove_cells])



###################################
### SAVING RAW Seurat.v3 OBJECT ###
###################################
cat("\n### Saving the Seurat object ###\n")
saveRDS(DATA, file = paste0(opt$output_path,"/seurat_object.rds") )
#---------



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
cat("\n\n\n\n")
cat("\n##############################")
cat("\n### SCRIPT RAN SUCESSFULLY ###")
cat("\n##############################")
cat("\n\n\n\n")

cat("\n##########################")
cat("\n### SYSTEM INFORMATION ###")
cat("\n##########################")
Sys.info()
cat("\n\n\n\n")

cat("\n###########################")
cat("\n### SESSION INFORMATION ###")
cat("\n###########################")
sessionInfo()
#---------


