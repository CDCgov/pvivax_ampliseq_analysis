library(cluster)
# library(tidyverse)
library(quantmod)
library(dplyr)
library(stringr)
library(filesstrings)
library(reshape)
# library(dbscan)
library(spaa)
require(parallel)
library(ape)
# library(utils)
# library(colorspace)
# library(DescTools)
#library(data.table)






#library(reticulate)

#UP November 19, 2020

date_today <- format(Sys.time(), '%Y-%m-%d_%H%M')

temporary_directory <- getwd()

args <- commandArgs(trailingOnly =T)


#### First you need to import the most recent ensemble matrix as a variable in R called "matrix" in R
# matrix_location <- readLines("MATRIX_LOCATION")

# setwd(matrix_location)

# myFiles <- list.files(pattern="*csv", all.files = FALSE)
# myMatrix <- (sort(myFiles, decreasing = TRUE)[1]) 
# matrix <- as.matrix(read.table(myMatrix, sep = ",", row.names=1, header=TRUE))

myMatrix <- args[1]
print(myMatrix)
matrix <- as.matrix(read.table(myMatrix, sep = ",", row.names=1, header=TRUE))


# setwd(temporary_directory)

##### Now you need to calculate the number of genetic clusters using the "gold standard" list. Use this to calculate the number of genetic clusters. "number_of_genetic_clusters" needs to be the name of this variable.


###set range of cluster numbers to test
min_cluster_number = as.numeric(readLines("CLUSTER_MIN"))
max_cluster_number = as.numeric(readLines("CLUSTER_MAX"))


### below here is where things change completely to the original:

folder <- readLines("SCRIPTS_LOC")
# print(folder)

source(paste(folder,"CLUSTER_COUNTER_V2.R", sep ="/"))
source(paste(folder,"CUTOFF_from_density.R", sep = "/"))


# source("../CLUSTER_COUNTER_V2.R")


# source("../CUTOFF_from_density.R") 




#NOW WE FIND CLUSTERS

#### First you need to import the most recent ensemble matrix as a variable in R called "matrix" in R


Ensemble_x <- as.hclust(agnes(x=matrix, diss = TRUE, stand = TRUE, method = "ward"))
Ensemble_y <- as.phylo(Ensemble_x)

write.tree(Ensemble_y, "tree.newick")

# run_in_system <- paste0("parnas -t tree.newick --cover --radius ", raw_threshold_for_PARNAS, " --clusters ",date_today,"_clusters.txt")

# system(run_in_system)



# cluster_location <- readLines("CLUSTER_LOCATION")

# clusters_to_move <- paste0(cluster_location,"/",date_today,"_clusters.txt")

# clusters_file <- paste0(date_today,"_clusters.txt")

# file.copy(clusters_file, clusters_to_move)

