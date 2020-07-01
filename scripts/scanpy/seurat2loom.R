# Convert from seurat to scanpy anndata
# OBS! Requires loomR
#devtools::install_github(repo = 'mojaveazure/loomR', ref = 'develop')


library(optparse)

cat("\nRunning Seurat to loom object conversion with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--input"),    type = "character",   metavar="character",  help="Path to RDS file with Seurat object"),
  make_option(c("-o", "--output"),     type = "character",   metavar="character",  help="Output loom file"),
  make_option(c("-f", "--force"),     action = "store_true",   default=TRUE, help="Force overwrite existing output file, default True")
)
opt = parse_args(OptionParser(option_list=option_list))

#opt <- list()
#opt$input <- "../datasets/bone_marrow/sauron/merged/marrow_all_seurat_object.rds"
#opt$output <- "../datasets/bone_marrow/sauron/merged/marrow_all_seurat_object.loom"

if (file.exists(opt$output) && opt$force){
  print(sprintf("Removing file %s",opt$output))
  file.remove(opt$output)
}

library(loomR)
library(Seurat)

seurat <- readRDS(opt$input)

# OBS! Need to have  variable genes slot.
if (length(seurat@assays$RNA@var.features) == 0) {
  seurat <- FindVariableFeatures(object = seurat) 
}

# if any meta data column is empty or all NA - gives error.
# remove such columns.
nC <- apply(seurat@meta.data, 2, function(x) sum(na.omit(nchar(x)))>0)
nN <- apply(seurat@meta.data, 2, function(x) sum(!is.na(x))>0)
seurat@meta.data <- seurat@meta.data[,nC & nN]

# if character column has NAs - gives error
# replace with "NA
for (cn in colnames(seurat@meta.data)){
  cl <- class(seurat@meta.data[,cn])
  nn <- sum(is.na(seurat@meta.data[,cn]))
  if (cl == "character" & nn>0){
    print(cn)
    temp <- seurat@meta.data[,cn]
    temp[is.na(temp)] <- "NA"
    seurat@meta.data[,cn] <- temp
  } else if (cl == "factor" & nn>0){
    print(cn)
    temp <- seurat@meta.data[,cn]
    levels(temp) <- c(levels(temp), "NA")
    temp[is.na(temp)] <- "NA"
    seurat@meta.data[,cn] <- temp
  }
}

lc <- as.loom(seurat,  filename = opt$output)
lc$close_all()

# Adding: tissue
# Error in if (size == Inf) { : missing value where TRUE/FALSE needed
# Calls: as.loom ... getDtype -> <Anonymous> -> <Anonymous> -> <Anonymous>

#tmp <- ReadH5AD("../datasets/bone_marrow/sauron/merged/marrow_all_seurat_object.h5ad")
#Pulling expression matrices and metadata
#Data is unscaled
#Error in file[["obs"]][] : 
#  object of type 'environment' is not subsettable
# does not work with seurat  3.1.3 either.

