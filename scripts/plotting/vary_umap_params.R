#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
library(optparse)
cat("\nRunning UMAP parameter variation with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),  type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-c", "--plot_groups"),         type = "character",   metavar="character",   default='none',  help="Column name(s) of FACTORS in the Metadata matrix by which to color cells. Multiple names should be separated by a comma."),
  make_option(c("-f", "--plot_features"),       type = "character",   metavar="character",   default='none',  help="Column name(s) of NUMERIC fields in the Metadata matrix, or gene names, by which to color cells. Multiple names should be separated by a comma"),
  make_option(c("-k", "--pre_dim_reduct"),      type = "character",   metavar="character",   default='pca',   help="The reduced dimension embedding upon which the UMAP will be run, such as 'pca', 'mnn', etc."),
  make_option(c("-b", "--umap_base_params"),    type = "character",   metavar="character",   default='default',  help="UMAP parameters to be used as the 'base' case, upon which parameter variations will be performed. Separate parameters by commas: 'param1=value1, param2=value2, ...'. For example: 'n.neighbors=10, min.dist=0.1'."),
  make_option(c("-v", "--umap_vary_params"),    type = "character",   metavar="character",   default='default',  help="UMAP parameters to vary, and their values. Parameters are separated by semicolons, and values by commas. For the 'dims' parameter, enter only the total number, not a vector (i.e., 50 will be interpreted as 1:50). For example: 'n.neighbors, 5, 10, 15, 30; dims, 5, 10, 50; spread, 1, 3, 5'."),
  make_option(c("-a", "--assay"),               type = "character",   metavar="character",   default='RNA',   help="Assay to be used in the analysis."),
  make_option(c("-m", "--add_metadata"),        type = "character",   metavar="character",   default='none',  help="Path to csv file containing additional cell metadata to add to the Seurat object (for plotting)."),
  make_option(c("-o", "--output_path"),         type = "character",   metavar="character",   default='none',  help="Output directory")
)
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if (!dir.exists(opt$output_path)) {dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)

if ((casefold(opt$plot_groups) == 'none') & (casefold(opt$plot_features) == 'none')) {opt$plot_groups='assay'}
#---------


######################
### LOAD LIBRARIES ###
######################
cat("\nLoading libraries ...\n")
suppressMessages(suppressWarnings({
  library(Seurat)
  library(dplyr)
  library(cowplot)
  # library(RColorBrewer)
  library(parallel)
}))
#---------



#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\n### LOADING Seurat.v3 OBJECT ###\n")
DATA <- readRDS(opt$Seurat_object_path)
DATA@active.assay <- opt$assay
#---------


####################################################
### RETRIEVE AND INCORPORATE ADDITIONAL METADATA ###
####################################################
if (casefold(opt$add_metadata) != 'none') {
  cat("\n### LOADING ADDITIONAL METADATA ###\n")
  add_data <- read.csv2(opt$add_metadata, row.names=1)
  if (!any(rownames(add_data) %in% rownames(DATA@meta.data))) {
    add_data <- as.data.frame(t(add_data))
  }
  DATA <- AddMetaData(DATA, add_data)
}
#---------



################
### RUN UMAP ###
################

# specify UMAP default parameters
default_params <- list(reduction=casefold(opt$pre_dim_reduct),
                       dims=1:30,
                       n.neighbors=30L,
                       n.components=2L,
                       metric="correlation",
                       n.epochs=200,
                       learning.rate=1,
                       min.dist=0.3,
                       spread=1,
                       set.op.mix.ratio=1,
                       local.connectivity=1L,
                       repulsion.strength=1,
                       negative.sample.rate=5L,
                       seed.use=42L,
                       reduction.key="UMAP_",
                       verbose=F)

# parse "base" UMAP params and update umap_params
if (!casefold(opt$umap_base_params) %in% c('default', 'defaults', 'none')) {
  base_params <- eval(parse(text=paste0("list(", opt$umap_base_params, ")")))
  if (('dims' %in% names(base_params)) & (length(base_params$dims) == 1)) {
    base_params$dims <- 1:base_params$dims
  }
  default_params <- modifyList(default_params, base_params)
}

# parse varied UMAP params
vary_params <- paste0(sub(",", "=c(", casefold(unlist(strsplit(opt$umap_vary_params, ";")))), ")", collapse=",")
vary_params <- eval(parse(text=paste0("list(", vary_params, ")")))
cat('\nVarying the following parameters and values:\n')
print(vary_params); cat('\n')

# parse plotting parameters
if (casefold(opt$plot_groups) != 'none') {
  pg <- trimws(unlist(strsplit(casefold(opt$plot_groups), ',')))
  plot_groups <- colnames(DATA@meta.data)[casefold(colnames(DATA@meta.data)) %in% pg]
  if (length(pg) > length(plot_groups)) {
    cat('WARNING! The following plot_groups were not found in the metadata:', pg[!(pg %in% casefold(plot_groups))])
  }
} else {plot_groups <- NULL}
if (casefold(opt$plot_features) != 'none') {
  all_feats <- union(rownames(DATA), colnames(DATA@meta.data))
  pf <- trimws(unlist(strsplit(casefold(opt$plot_features), ',')))
  plot_feats <- all_feats[casefold(all_feats) %in% pf]
  if (length(pf) > length(plot_feats)) {
    cat('WARNING! The following plot_features were not found in the metadata:', pf[!(pf %in% casefold(plot_feats))])
  }
} else {plot_feats <- NULL}

# iterate through parameters and parameter values
for (vparam in names(vary_params)) {
  cat('\n\nVarying', vparam, '\n')
  param_vals <- vary_params[[vparam]]
  tmp <- DATA
  pgroup <- pfeat <- list()
  for (v in param_vals) {
    upar <- default_params
    if (vparam=='dims') {
      upar[[vparam]] <- 1:v
    } else {
      upar[[vparam]] <- v
    }
    cat('Generating UMAP with:', vparam, '=', v, '... \n')
    
    # long command because "do.call()" is extremely slow with large objects
    tmp <- RunUMAP(object=tmp, reduction.name='umap', dims=upar$dims,
                   reduction=upar$reduction, n.neighbors=upar$n.neighbors, n.components=upar$n.components,
                   metric=upar$metric, n.epochs=upar$n.epochs, learning.rate=upar$learning.rate, min.dist=upar$min.dist,
                   spread=upar$spread, set.op.mix.ratio=upar$set.op.mix.ratio, local.connectivity=upar$local.connectivity,
                   repulsion.strength=upar$repulsion.strength, negative.sample.rate=upar$negative.sample.rate,
                   seed.use=upar$seed.use, reduction.key=upar$reduction.key, verbose=upar$verbose)
    
    # generate plots
    col_scale <- c("grey85","navy")
    if (!is.null(plot_groups)) {
      pgroup <- c(pgroup, list(DimPlot(tmp, reduction='umap', group.by=plot_groups, ncol=length(plot_groups))))
    }
    if (!is.null(plot_feats)) {
      pfeat <- c(pfeat, list(FeaturePlot(tmp, reduction='umap', features=plot_feats, ncol=length(plot_feats), order=T, cols=col_scale, min.cutoff=0, max.cutoff=1)))
    }
  }

  # combine plots
  reduction_names <- paste(vparam, param_vals, sep='_')
  if (!is.null(plot_groups)) {
    png(filename=paste0(opt$output_path, '/umap_dimplots_', vparam, '.png'), units='mm', res=150, height=100*length(param_vals), width=120*length(plot_groups))
    print(plot_grid(plotlist=pgroup, ncol=1, labels=reduction_names, scale=0.9))
    invisible(dev.off())
  }
  if (!is.null(plot_feats)) {
    png(filename=paste0(opt$output_path, '/umap_featplots_', vparam, '.png'), units='mm', res=150, height=100*length(param_vals), width=120*length(plot_feats))
    print(plot_grid(plotlist=pfeat, ncol=1, labels=reduction_names, scale=0.9))
    invisible(dev.off())
  }
  invisible(gc())
}
cat('\nDone!\n\n')



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
print_session_info()
#---------







