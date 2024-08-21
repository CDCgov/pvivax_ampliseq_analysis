
library(ggplot2)
library(cluster)
library(ggtree)
library(ape)
library(RColorBrewer)


args <- commandArgs(trailingOnly =T)

#Read in matrix, clusters, and geographic prediction DB
matrix <- as.matrix(read.table(args[1], sep = ",", row.names=1, header=TRUE))
parnasClusters <- read.table(args[2], sep = "\t", header =F)
geographic <- read.table(args[3], sep ="\t", header =T)
travelHistory <- read.table(args[4], sep = "\t", header =T)
print("after read in")

stateID <- args[5]

geographic$Region_Prediction_1 <- gsub(" ", "", geographic$Region_Prediction_1)
colnames(parnasClusters) <- c("Seq_ID", "Cluster")

#parnas genetic clusters start at 0, make them start at 1 which is a little easier on the eye (just make sure you remember the numbers have increased by 1 when looking at other files)
parnasClusters$Cluster <- parnasClusters$Cluster +1

parnasClusters <- merge(parnasClusters, geographic, by.x = "Seq_ID", by.y = "Seq_ID", all.x =T)
parnasClusters <- merge(parnasClusters, travelHistory, by.x ="Seq_ID", by.y = "Seq_ID_final",  all.x = T)

print(parnasClusters)
#The number of unique clusters - used for highlighting each cluster in different colors
#Right now, it is random color assignment, but we could change so each numbered cluster is the same in each tree
correct_number_of_clusters <- length(unique(parnasClusters$Cluster))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
col=sample(color, correct_number_of_clusters)

#make sure that the matrix (used for making the tree) and the annotation text have the same order - vital for making sure tree is labeled correctly
Seq_ID <- as.data.frame(rownames(matrix))
Row <- as.data.frame(rownames(Seq_ID))
row_order <- cbind(Row, Seq_ID)
corrected_orderTest <- parnasClusters[order(match(parnasClusters[,1], row_order[,2])),]

#maket the tip label text
New_ID <- paste0(corrected_orderTest$Seq_ID,"_Cluster_",corrected_orderTest$Cluster, "_Region_", corrected_orderTest$Region_Prediction_1, "_Travel_", corrected_orderTest$Travel, sep ="")

#Get a list of all states represented in the dataset, this will be for generating reports
startFilter <- New_ID[!grepl("REF", New_ID)]
#old name format starts with the state
oldNames <- substr(startFilter,1,2)
oldNames <- oldNames[!grepl("^DM", oldNames)]
oldNames <- oldNames[!grepl("XX", oldNames)]

#new name format starts with DM
newNames <- startFilter[grepl("^DM", startFilter)]
newNames <- substr(newNames,6,7)

#combine all state names into one
uniqStates <- unique(c(oldNames,newNames))


#Generate phylo object for tree generation with updated labels
Ensemble <- data.frame(matrix)
rownames(Ensemble) <- c(New_ID)
Ensemble <- (as.matrix(Ensemble))
Ensemble_phylo <-  as.phylo(as.hclust(agnes(x=Ensemble, diss = TRUE, stand = TRUE, metric = "manhattan", method = "ward")))


#Start building the base tree. Get the region prediction for each sample, extracted as the tip label
tip <- Ensemble_phylo$tip.label
myRegions <- list(
    "AF" = tip[grepl("Region_AF", tip)],
    "LAM" = tip[grepl("Region_LAM", tip)],
    "WAS" = tip[grepl("Region_WAS", tip)],
    "WSEA" = tip[grepl("Region_WSEA", tip)],
    "ESEA" = tip[grepl("Region_ESEA", tip)],
    "MSEA" = tip[grepl("Region_MSEA", tip)],
    "EAS" = tip[grepl("Region_EAS", tip)],
    "OCE" = tip[grepl("Region_OCE", tip)]
)

# keepLabs <- c("DM23", "TX0002", )

#Basic tree
q <- ggtree(Ensemble_phylo, size = 0.9, layout = "circular")

#use groupOTU with each set of tip labels grouped by the region prediction - this will color the branches. I had to play around with making the tree to see what order i needed to specify the colors
q <- groupOTU(q, myRegions, "Regions") + aes(color=Regions) + scale_colour_manual(values = c("black","red", "blue", "green", "purple", "orange", "yellow", "pink", "darkgray")) + theme(legend.position = "top") 

#To add tip labels - as is - to the tree (for all samples)
# q <- q + geom_tiplab(color = "black", size = 1, offset = 0.60)

#To only add a label for the labels specified in the grep statements
q <- q + geom_tiplab(aes(subset = (grepl(stateID, label, fixed = T)==TRUE)), color = "black", size = 1, offset = 0.50)  
#Should be be to combine the above grepl statements into a single string but its not working for some reason
# keepLabs <- c("DM23", "TX0002")
# q <- q + geom_tiplab(aes(subset = (grepl(paste(keepLabs, collapse = "|"), label, fixed = T)==TRUE)), color = "black", size = 2, offset = 0.60) 

p <- q

#Put the tip labels into a dataframe - used for downstream annotation
tip_labels <- data.frame(Ensemble_phylo$tip.label)
tip_labels <- cbind(rownames(tip_labels), tip_labels)
colnames(tip_labels) <- c("Tips","Specimen")

#Loop through each cluster and assign a different color, based on the colors selected previously
for (i in 1: (correct_number_of_clusters)){

	#extract the tip numbers for each specimen in the given cluster
	clu_num_working <- paste("Cluster_",i,"_",sep="")
	foobar <- tip_labels[grep(clu_num_working, tip_labels$Specimen), ]
	foobar <- rownames(foobar)
	foobar <- as.numeric(foobar)

	#extract the color from col vector
	color_of_strips <- col[i]

	p <- p + geom_hilight(node=c(foobar), fill= color_of_strips, extend = 0.45, alpha = 1)

}

g <- p



time_now <- format(Sys.time(), "%H%M")
date_now <- Sys.Date()
name_tree <- paste(date_now,"_", time_now, "_", stateID, "_Pvivax_ampliseq_tree.pdf", sep="")
ggsave(g, file=name_tree, width=50, height=50, units = "cm")

# #This will be the loop for generating state specific trees
# n <- p + geom_tiplab(aes(subset = (grepl("TX", label, fixed = T)==TRUE)), color = "black", size = 2, offset = 0.60)
# name_tree <- paste(date_now,"_", time_now, "_Pvivax_ampliseq_clustersColored_TX_tree.pdf", sep="")
# ggsave(n, file=name_tree, width=50, height=50, units = "cm")


#################################################################################################################################
#if you want to make a super simple tree with no colors and just the specimen name on the tip labels

# Ensemble_phylo <-  as.phylo(as.hclust(agnes(x=matrix, diss = TRUE, stand = TRUE, metric = "manhattan", method = "ward")))
# pBasic <- ggtree(Ensemble_phylo, size = 0.9, layout = "circular") + geom_tiplab(color = "black", size = 2, offset = 0.05)
# name_treeBasic <- paste(date_now,"_", time_now, "_Pvivax_ampliseq_basic_tree.pdf", sep="")
# ggsave(pBasic, file=name_treeBasic, width=50, height=50, units = "cm")



















