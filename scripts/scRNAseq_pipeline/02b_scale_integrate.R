#!/usr/bin/env Rscript

#############################
### LOAD/INSTALL OPTPARSE ###
#############################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')};
library(optparse)
#---------

##################################
### DEFINE PATH TO LOCAL FILES ###
##################################
cat("\nRunning DATA SCALING AND INTEGRATION with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-j", "--sauron_script_path"),    type = "character",   metavar="character",   default='none',  help="Path to Sauron scripts"),
  make_option(c("-r", "--regress"),               type = "character",   metavar="character",   default='none',  help="Variables to be regressed out using linear modeling."),
  make_option(c("-b", "--integration_method"),    type = "character",   metavar="character",   default='cca,orig.ident',  help="Integration method to be used. 'cca', mmn', 'scale' and 'combat' are available at the moment. The batches (column names in the metadata matrix) to be removed should be provided as arguments comma separated. E.g.: 'Combat,sampling_day'. For MNN, an additional integer parameter is supplied as the k-nearest neighbour."),
  make_option(c("-v", "--var_genes"),             type = "character",   metavar="character",   default='scran,.2',  help="Whether use 'Seurat' or the 'Scran' method for variable genes identification. An additional value can be placed after a comma to define the level of dispersion wanted for variable gene selection. 'Seurat,2' will use the threshold 2 for gene dispersions. Defult is 'scran,0.001'. For Scran, the user should inpup the level of biological variance 'Scran,0.001'. An additional blocking parameter (a column from the metadata) can ba supplied to 'Scran' method block variation comming from uninteresting factors, which can be parsed as 'Scran,0.2,Batch'."),
  make_option(c("-N", "--nSeurat"),               type = "character",   metavar="character",   default='3000',  help="Number of var genes per dataset (as identified by seurat var gene selection) to use for the integration step."),
  make_option(c("-s", "--cluster_use"),           type = "character",   metavar="character",   default='all',  help="The cluster to be used for analysis."),
  make_option(c("-a", "--assay"),                 type = "character",   metavar="character",   default='RNA',  help="The default assay to use to integrate."),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output directory")
) 
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)

col_scale <- c("grey85","navy")
VAR_choice <- as.character(unlist(strsplit(opt$var_genes,",")))
#---------



###################################
### SETUP MULTICORE ENVIRONMENT ###
###################################
#options(future.globals.maxSize= 2048 * 1024^2)
#plan(strategy = "multicore", workers = nbrOfWorkers())
#---------



##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
initial.options <- commandArgs(trailingOnly = FALSE)
script_path <- opt$sauron_script_path
source( paste0(script_path,"/compute_hvgs.R") )
source( paste0(script_path,"/fast_ScaleData.R") )
source( paste0(script_path,"/inst_packages.R") )

suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(biomaRt)))
suppressMessages(suppressWarnings(library(igraph)))
suppressMessages(suppressWarnings(library(sva)))
suppressMessages(suppressWarnings(library(rafalib)))
suppressMessages(suppressWarnings(library(parallel)))
suppressMessages(suppressWarnings(library(scran)))
suppressMessages(suppressWarnings(library(scater)))
#---------



#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\n### LOADING Seurat.v3 OBJECT ###\n")
DATA <- readRDS(opt$Seurat_object_path)
DATA@active.assay <- opt$assay
#---------



###################################
### SELECT CELLS FROM A CLUSTER ###
###################################
cat("\n### SELECTING CELLS FROM A CLUSTER ###")
if (length(unlist(strsplit(opt$cluster_use,","))) >= 2 ){
  clustering_use <- as.character(unlist(strsplit(opt$cluster_use,",")))[1]
  clusters_to_select <- as.character(unlist(strsplit(opt$cluster_use,",")))[-1]
  if(!(clustering_use %in% colnames(DATA@meta.data))){
    cat("\nThe Clustering level was not found in the Metadata...\n")
    cat("\nAll cells will be used ...\n")
    cells_use <- rownames(DATA@meta.data)
  } else if(sum(clusters_to_select %in% unique(DATA@meta.data[,clustering_use])) == 0){
    cat("\nThe Cluster specifed was not found ...\n")
    cat("\nAll cells will be used ...\n")
    cells_use <- rownames(DATA@meta.data)
  } else {
    cells_use <- rownames(DATA@meta.data)[factor(DATA@meta.data[,clustering_use]) %in% clusters_to_select]}   #Filter out cells with no assigned clusters
  DATA <- SubsetData(DATA,assay = opt$assay,cells = cells_use)
} else {
  cat("\nThe name of the cluster or the cluster name were not found in your data.\n All cells will be used ...\n")
  cells_use <- rownames(DATA@meta.data)}
# sel <- rowSums(as.matrix(DATA@assays[[opt$assay]]@counts) >= 1) >= 1
# DATA <- CreateSeuratObject(as.matrix(DATA@assays[[opt$assay]]@counts[sel,cells_use]), meta.data = DATA@meta.data[cells_use,])
#---------



########################################
### SCALE DATA AND REGRESS VARIABLES ###
########################################
cat("\n### Scaling data and regressing uninteresting factors ###\n")
if (casefold(opt$regress) == 'none') {
  vars_regress <- NULL
} else {
  vars_regress <- unlist(strsplit(opt$regress,","))
  vars_missing <- !(vars_regress %in% colnames(DATA@meta.data))
  if (any(vars_missing)) { cat('WARNING! The following variables were not found in the metadata:', vars_regress[vars_missing])}
  if (all(vars_missing)) { vars_regress <- NULL}
}
DATA <- fast_ScaleData(DATA, vars.to.regress=vars_regress, assay=opt$assay)
# DATA <- ScaleData(DATA, vars.to.regress = vars_regress, assay=opt$assay)
#---------



#######################################
### INTEGRATE DATASETS USING COMBAT ###
#######################################
integration_method <- unlist(strsplit(opt$integration_method,","))

if ((length(integration_method) >= 2) & (casefold(integration_method[1]) == "combat") ){
  cat("\n### INTEGRATING DATASETS USING COMBAT ###\n")
  print(integration_method)
  
  #Defining batch variables
  batch <- factor(DATA@meta.data[,integration_method[2]])
  mod0 <- model.matrix(~1, data=as.data.frame(DATA@meta.data))
  
  #Get scaled count data from Seurat object
  scaled_data <- as.matrix(DATA@assays[[opt$assay]]@scale.data)[, rownames(DATA@meta.data)]
  
  #Applying Combat
  combat_data <- ComBat(dat=scaled_data, batch=batch, mod=mod0)
  
  #Add integrated data as new assay in Seurat object
  DATA@assays[["combat"]] <- CreateAssayObject(data=combat_data, min.cells=0, min.features=0)
  rm(combat_data,scaled_data,mod0);  invisible(gc())
}
#---------



####################################
### INTEGRATE DATASETS USING MNN ###
####################################
if ((length(integration_method) >= 2) & (casefold(integration_method[1]) == "mnn") ){
  cat("\n### INTEGRATING DATASETS USING MNN ###\n")
  print(integration_method)
  
  #Defining batch variables
  cat("\nCreating dataset list\n")
  batch <- as.character(factor(DATA@meta.data[,integration_method[2]]))
  DATA.list <- SplitObject(DATA, split.by = integration_method[2])
  
  if( (length(DATA.list) > 1) ){
    
    # define HVGs per dataset
    cat("\nComputing HVGs\n")
    for (i in 1:length(DATA.list)) {
      cat("\nProcessing dataset: ",names(DATA.list)[i]," \n")
      #DATA.list[[i]] <- NormalizeData(DATA.list[[i]], verbose = FALSE,scale.factor = 1000)
      if(file.exists(paste0(opt$output_path,"/variable_genes/var_genes_",names(DATA.list)[i],"/HVG_info_",VAR_choice[1],".csv"))){
        cat("\nVariable genes for this dataset found, and will be used. Skiping HVG calculation.\n")
        temp <- read.csv2(paste0(opt$output_path,"/variable_genes/var_genes_",names(DATA.list)[i],"/HVG_info_",VAR_choice[1],".csv"),row.names=1)
        DATA.list[[i]]@assays[[opt$assay]]@meta.features <- temp
        DATA.list[[i]]@assays[[opt$assay]]@var.features <- rownames(temp)[temp$use]
      } else {DATA.list[[i]] <- compute_hvgs(DATA.list[[i]],VAR_choice,paste0(opt$output_path,"/variable_genes/var_genes_",names(DATA.list)[i]),assay = opt$assay, nSeurat=as.numeric(opt$nSeurat))}
    }
    
    # select the most informative genes that are shared across all datasets:
    cat("\nComputing HVGs\n")
    universe <- unique(unlist(lapply(DATA.list,function(x){x@assays[[opt$assay]]@var.features})))
    head(universe, 50)
    cat("\n",length(universe)," genes found as variable within datasets\n")
    
    #Separating batch matricies
    DATA.list <- lapply(DATA.list, function(x){ x@assays[[opt$assay]]@scale.data[universe,]} )
    myinput <- list()
    
    if (is.na(integration_method[3])) {
      myinput[["k"]] <- 20
    } else {
      myinput[["k"]] <- as.numeric(integration_method[3])
    }
    myinput[["approximate"]] <-  TRUE
    myinput[["d"]] <-  51
    #myinput[["BPPARAM"]] <-  MulticoreParam(workers = detectCores()-1)
    
    #Applying MNN correction on scaled counts
    cat("\nApplying MNN correction on scaled counts\n")
    out <- do.call(fastMNN, args=c(DATA.list, myinput))
    cat("\nMNN computation done\n")
    out <- t(out$corrected)
    colnames(out) <- unlist(lapply(DATA.list,function(x){colnames(x)}))
    out <- out[,colnames(DATA)]
    rownames(out) <- paste0("dim",1:myinput$d)
    DATA@reductions[["mnn"]] <- CreateDimReducObject(embeddings = t(out),key = "MNN_",assay = opt$assay)
    rm(out, myinput);  invisible(gc())
  }
}
#---------





####################################
### INTEGRATE DATASETS USING CCA ###
####################################
if ((length(integration_method) >= 1) & (casefold(integration_method[1]) == "cca") ){
  cat("\n### INTEGRATING DATASETS USING CCA ###\n")
  print(integration_method)
  
  DATA.list <- SplitObject(DATA, split.by = integration_method[2])
  if( (length(DATA.list) > 1) ){
    
    for (i in 1:length(DATA.list)) {
      #DATA.list[[i]] <- NormalizeData(DATA.list[[i]], verbose = FALSE,scale.factor = 1000)
      DATA.list[[i]] <- compute_hvgs(DATA.list[[i]],VAR_choice,paste0(opt$output_path,"/var_genes_",names(DATA.list)[i]),assay = opt$assay, nSeurat=as.numeric(opt$nSeurat))
      gc()
    }
    universe <- unique(unlist(lapply(DATA.list,function(x){x@assays[[opt$assay]]@var.features})))
    
    # "scale" is set to FALSE in the FindIntegrationAnchors call below since we have already scaled the data above
    DATA.anchors <- FindIntegrationAnchors(object.list = DATA.list, dims = 1:30, anchor.features = universe, scale=FALSE)
    DATA <- IntegrateData(anchorset = DATA.anchors, dims = 1:30, new.assay.name = "cca")
    DATA@assays$cca@var.features <- rownames(DATA@assays$cca@data)
    rm(DATA.list,DATA.anchors); gc()
  }
}
#---------



###########################
### FIND VARIABLE GENES ###
###########################
if(DefaultAssay(DATA) == opt$assay){
  output_path <- paste0(opt$output_path,"/variable_genes/All_datasets_together")
  DATA <- compute_hvgs(DATA, VAR_choice, output_path, assay=opt$assay, nSeurat=as.numeric(opt$nSeurat))}
#---------



###################################
### SAVING Seurat.v3 OBJECT ###
###################################
cat("\n### Saving Seurat object ###\n")
saveRDS(DATA, file = paste0(opt$output_path,"/seurat_object.rds") )
#---------



#############################
### SYSTEM & SESSION INFO ###
#############################
print_session_info()
#---------

