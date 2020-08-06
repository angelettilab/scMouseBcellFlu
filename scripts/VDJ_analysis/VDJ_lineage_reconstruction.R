


library(alakazam)
library(igraph)
library(dplyr)


# specify directory containing genotyped VDJ files
proj_dir <- '/Users/jonrob/Documents/NBIS/LTS_projects/d_angeletti_1910/scMouseBcellFlu/'
geno_dir <- paste0(proj_dir, 'analysis/immcantation/genotyping/')

# get list of available files
geno_files <- dir(geno_dir, 'germ-pass[.]tab', full.names=T)

# specify which file to load
g_file <- geno_files[1]
  
# load the Change-O database file with germline sequence information (*_germ-pass.tab file)
db <- readChangeoDb(g_file)

# select a desired clone
head(sort(table(db$CLONE), decreasing=T))  # view most abundant clone IDs
sub_db <- subset(db, CLONE == 681)

# create ChangeOclone object for clone
clone <- makeChangeoClone(sub_db, text_fields=c('SEQUENCE_ID','C_CALL'), num_field='UMICOUNT')

# Run PHYLIP and parse output
# Need to download and install PHYLIP from http://evolution.genetics.washington.edu/phylip/getme-new1.html
phylip_exec <- "/Users/jonrob/Documents/NBIS/repos/phylip-3.695/exe/dnapars.app/Contents/MacOS/dnapars"  # path to dnapars executable within PHYLIP package
graph <- buildPhylipLineage(clone, phylip_exec, rm_temp=T, verbose=T)

# Modify graph and plot attributes
V(graph)$color <- "steelblue"
V(graph)$color[V(graph)$name == "Germline"] <- "black"
V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
V(graph)$label <- V(graph)$C_CALL
E(graph)$label <- ""

# Remove large default margins
par(mar=c(0, 0, 0, 0) + 0.1)

# Plot graph
plot(graph, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
     vertex.label.color="black", vertex.size=5, vertex.label.family='sans', vertex.label.cex=0.7)

# Add legend
legend("topleft", c("Germline", "Inferred", "Sample"), 
       fill=c("black", "white", "steelblue"), cex=0.75)




