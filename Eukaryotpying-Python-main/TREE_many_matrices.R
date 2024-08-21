#########LARGE DATASET ANALYSIS VERSION 2
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(dbscan)
library(gridExtra)  
library(stringr)
library(cluster)
library(phylogram)
library(msa)
library(ggtree)
library(colorspace)
library(ape)
library(seqinr)
library(colorspace)
library(phangorn)
library(randomcoloR)
library(dynamicTreeCut)

setwd("/Users/joelbarratt/Documents/P_VIVAX_PROJECT/GENOME_PIPELINE/CDC-Complete-Cyclospora-typing-workflow-ALPHA-TEST-master/Complete_Cyclospora_typing_workflow_MacOS_High_Sierra_ALPHA_TEST/Eukaryotpying-Python-main/TMP_matrix/")

all_files <- list.files("MATRICES")

#TABLE <- all_files[grepl("ALL_CHR_pairwisedistancematrix_H.csv", all_files)]

#TABLE <- tail(all_files, 1)




for( k in 1:length(all_files)){

matrix <- as.matrix(read.table(paste0("MATRICES/", all_files[k]), sep = ",", row.names=1, header=TRUE))
Ensemble_x <- as.hclust(agnes(x=matrix, diss = TRUE, stand = TRUE, metric = "manhattan", method = "ward"))
Ensemble_y <- as.phylo(Ensemble_x)

p <- ggtree(Ensemble_y, size = 0.9, layout = "circular") + geom_tiplab(color = "black", size = 2, offset = 0.05)

x <- p + ggplot2::xlim(0.0,2.0)

filename <- paste0(all_files[k],".pdf")

pdf(file = filename,   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
    
print(x)

dev.off()

}

















