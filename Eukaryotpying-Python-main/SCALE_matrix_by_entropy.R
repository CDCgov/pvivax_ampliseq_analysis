#########LARGE DATASET ANALYSIS VERSION 2
# library(ggplot2)
# library(RColorBrewer)
# library(reshape)
# library(dbscan)
# library(gridExtra)  
library(stringr)
# library(cluster)
# library(phylogram)
# library(msa)
# library(ggtree)
# library(colorspace)
# library(ape)
# library(seqinr)
# library(colorspace)
# library(phangorn)
# library(randomcoloR)
# library(dynamicTreeCut)


args <- commandArgs(trailing = T)
setwd("MATRICES")
matrices <- list.files()
df_list <- list()
print("start matrix loop")
for(k in matrices){
	print(k)
	newName <- str_split(k, "_", simplify= T)[3]
	print(newName)

	x <- read.table(k, sep = ",", row.names=1, header=TRUE)	

	df_list[[newName]] <- x
	
}
print(names(df_list))
# names(df_list) <- paste0("CHR",c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"))

### now we want to calculate the entropy of each chromosome to see which is more informative.

setwd("../STARTING_HAP_SHEET")


### get the set of haplotype sheets -- one for each chromosome.
hap_sheet_list <- list()
hapFiles <- list.files()
print("start hap sheet loop")
for(k in hapFiles){
	print(k)
	newName <- str_split(k, "_", simplify= T)[3]
	print(newName)
	
	data = read.csv(k, skip=0, stringsAsFactors = FALSE, sep = "\t")
		
	hap_sheet_list[[newName]] <- data
	
}

print(names(hap_sheet_list))
# names(hap_sheet_list) <- paste0("CHR",c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"))


# for(hap_sheets in names(df_list)){
# 	# setwd(hap_sheets)
	
# 	data = read.csv(list.files(), skip=0, stringsAsFactors = FALSE, sep = "\t")
	
# 	hap_sheet_list[[hap_sheets]] <- data
	
# 	setwd("..")
# }



#entropy in bans for each chromosome.
####################################################################################################################################

chromosome_entropies <- list()


for(a in names(hap_sheet_list)){


markers_captured <- unique(substr(colnames(hap_sheet_list[[a]]), 1,27)) #1 to 27 cuts the column headings down to their parts.
markers_captured <- markers_captured[2:length(markers_captured)] #this will contain a set of all markers on this chromosome


     markers_on_this_chromosome_entropies <- list()
     
            for(b in markers_captured){

this_marker <- dplyr::select(hap_sheet_list[[a]],contains(b)) ## lets find the columns for the first marker

this_marker_list <- list()

for(k in 1:length(colnames(this_marker))){
	this_marker_list[[k]] <- sum(str_count(this_marker[,k], pattern = "X"))
}

denominator <- Reduce(`+`, this_marker_list)


for(l in 1:length(colnames(this_marker))){
	
	this_marker_list[[l]] <- (sum(str_count(this_marker[,l], pattern = "X")))/denominator

	}


for(m in 1:length(colnames(this_marker))){
	
	this_marker_list[[m]] <- this_marker_list[[m]]*log(this_marker_list[[m]], 10)

	}

entropy_this_marker <- (Reduce(`+`, this_marker_list))*(-1)


markers_on_this_chromosome_entropies[[b]] <- entropy_this_marker

                                }                          
  sum_markers_on_this_chromosome_entropies <- Reduce(`+`, markers_on_this_chromosome_entropies)
  
  chromosome_entropies[[a]] <- sum_markers_on_this_chromosome_entropies                       

}

####################################################################################################################################


#now we multiple each distance matrix by its entropy.

#chromosome_entropies

for(z in names(df_list)){
	
	df_list[[z]] <- df_list[[z]]* chromosome_entropies[[z]]
}





matrix_sum <- Reduce(`+`, df_list)# takes the average of all the dataframes in a list.

normalized_matrix <- matrix_sum/max(matrix_sum)



todays_date <- Sys.Date()

# time_now <- format(Sys.time(), "%X")
# TIME  <- str_replace_all(time_now,":", "") 


time_now <- format(Sys.time(), "%H%M")
date_now <- Sys.Date()

name_current_distance_matrix <- paste("../",date_now,"_", time_now, "_ALL_CHR_pairwisedistancematrix_H.csv", sep="") ##modify so it puts in the right place.

write.csv(normalized_matrix, name_current_distance_matrix, quote = F)

setwd("..")

# this_name <- paste(todays_date,"_", TIME, "_ALL_CHR_pairwisedistancematrix_H.csv", sep="")

# file.copy(this_name, "..")
# setwd("..")
# file.copy(this_name, "../ensemble_matrices/")

# file.remove(this_name)

#matrix <- as.matrix(normalized_matrix)
#Ensemble_x <- as.hclust(agnes(x=matrix, diss = TRUE, stand = TRUE, method = "ward"))
#Ensemble_y <- as.phylo(Ensemble_x)

#p <- ggtree(Ensemble_y, size = 0.9, layout = "circular") + geom_tiplab(color = "black", size = 2, offset = 0.05)
#p + ggplot2::xlim(0.0,2.0)




























