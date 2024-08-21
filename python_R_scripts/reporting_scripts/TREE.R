
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

geographic$Region_Prediction_1 <- gsub(" ", "", geographic$Region_Prediction_1)
colnames(parnasClusters) <- c("Seq_ID", "Cluster")

#parnas genetic clusters start at 0, make them start at 1 which is a little easier on the eye (just make sure you remember the numbers have increased by 1 when looking at other files)
parnasClusters$Cluster <- parnasClusters$Cluster +1

parnasClusters <- merge(parnasClusters, geographic, by.x = "Seq_ID", by.y = "Seq_ID", all.x =T)
parnasClusters <- merge(parnasClusters, travelHistory, by.x ="Seq_ID", by.y = "Seq_ID_final",  all.x = T)

# print(parnasClusters)
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
# q <- q + geom_tiplab(aes(subset = (grepl("DM23_TX0008", label, fixed = T)==TRUE)), color = "black", size = 2.5, offset = 0.50)  + geom_tiplab(aes(subset = (grepl("TX0001_23", label, fixed = T)==TRUE)), color = "black", size = 2.5, offset = 0.50) + geom_tiplab(aes(subset = (grepl("TX0003_23_1PvD", label, fixed = T)==TRUE)), color = "black", size = 2.5, offset = 0.50)

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
#The next steps are to identify the node with a given label and then label those nodes with a different shape (or letter or whatever). 
#This is a manual process on an as needed basis, not typically done. I print the node to screen and then go to the geom_strip and assign manually
# AR_node <-tip_labels[grep("DM23_AR0001", tip_labels$Specimen),1]
# # print("TX_8node")
# # print(TX_8node)
# # TX_9node <-tip_labels[grep("DM23_TX0009", tip_labels$Specimen),1]
# # print("TX_9node")
# # print(TX_9node)
# TX_8node <-tip_labels[grep("DM23_TX0008", tip_labels$Specimen),1]
# print("TX_8node")
# print(TX_8node)
# # TXA_node <-tip_labels[grep("TX000A", tip_labels$Specimen),1]
# # print("TXA_node")
# # print(TXA_node)
# # TXB_node <-tip_labels[grep("TX000B", tip_labels$Specimen),1]
# # print("TXB_node")
# # print(TXB_node)
# TX1_node <-tip_labels[grep("TX0001_23", tip_labels$Specimen),1]
# print("TX1_node")
# print(TX1_node)
# # TX2_node <-tip_labels[grep("TX0002_23", tip_labels$Specimen),1]
# # print("TX2_node")
# # print(TX2_node)
# # TX4_node <-tip_labels[grep("TX0004_23", tip_labels$Specimen),1]
# # print("TX4_node")
# # print(TX4_node)
# TX3_node <-tip_labels[grep("TX0003_23_1PvD", tip_labels$Specimen),1]
# print("TX3_node")
# print(TX3_node)
# MN1_node <-tip_labels[grep("MN0001_23", tip_labels$Specimen),1]
# print("MN1_node")
# print(MN1_node)
# FL0010_node <-tip_labels[grep("FL0010", tip_labels$Specimen),1]
# print("FL0010_node")
# print(FL0010_node)
# FL4_node <-tip_labels[grep("FL0004_XX", tip_labels$Specimen),1]
# print("FL4_node")
# print(FL4_node)
# FL417_node <-tip_labels[grep("FL0004_17", tip_labels$Specimen),1]
# print("FL4_node")
# print(FL417_node)
# FL16_node <-tip_labels[grep("FL0016_XX", tip_labels$Specimen),1]
# print("FL16_node")
# print(FL16_node)
# DC_node <-tip_labels[grep("DC0001_23", tip_labels$Specimen),1]
# MN_node <-tip_labels[grep("MN0001_23", tip_labels$Specimen),1]
# AZ_node <-tip_labels[grep("AZ0001_23", tip_labels$Specimen),1]

# PAmarko <-tip_labels[grep("PA0001_20", tip_labels$Specimen),1]
# FLmarko <-tip_labels[grep("FL0001_13", tip_labels$Specimen),1]
# GA06 <-tip_labels[grep("GA0001_06", tip_labels$Specimen),1]
# # print("FL16_node")
# # print(FL16_node)
# # combineTX1 <- c(TXA_node, TXB_node, TX2_node, TX1_node)
# # combineTX2 <- c(TX3_node, TX4_node, MN_node)

# # combineTravel <- c(FL0010_node, FL4_node, FL16_node, FL417_node, DC_node, AZ_node)

# FLstart <-tip_labels[grep("FL00", tip_labels$Specimen),]
# # print(FLstart)
# FL2003 <- FLstart[grepl("_03", FLstart$Specimen),1]
# FL2023 <- FLstart[grepl("_23_1PvD", FLstart$Specimen),1]


# print("florida to add")
# print(FL2003)
# # print(FL2023)
# # fl2023_nodes <- 

# GAControl_1 <-tip_labels[grep("GA0001_15", tip_labels$Specimen),1]
# # print("Control 1")
# # print(GAControl_1)

# GAControl_2 <-tip_labels[grep("GA0001_Pv", tip_labels$Specimen),1]
# # print("Control 2")
# # print(GAControl_2)

# GAControl_3 <-tip_labels[grep("2415.2_Sal", tip_labels$Specimen),1]
# # print("Control 3")
# # print(GAControl_3)

# controls <- c(GAControl_1, GAControl_2, GAControl_3)

# g <- p +geom_strip(23,23,label = "1", offset = 0.51)
# g <- g +geom_strip(187,188,label = "2", offset = 0.51, offset.text = 0.05, barsize = 2, color = "blue")
# g <- g +geom_strip(181,182,label = "2", offset = 0.51, offset.text = 0.05, barsize = 2, color = "blue")
# g <- g +geom_strip(184,186,label = "3", offset = 0.51, offset.text = 0.05, barsize = 2, color = "green")
# g <- g +geom_strip(185,185,label = "3", offset = 0.51, color = "green")
# g <- g +geom_strip(106,107,label = "4", offset = 0.51, offset.text = 0.05, barsize = 2, color = "purple")
# g <- g +geom_strip(93,94,label = "5", offset = 0.51 , offset.text = 0.05, barsize = 2,, color = "red")
# g <- g +geom_strip(113,113,label = "5", offset = 0.51, color= "red")

# g <- p +geom_strip("GA0001_15_2PvX_Cluster_150_Region_LAM(1)","GA0001_15_1PvX_Cluster_150_Region_LAM(1)",label = "Control", offset = 0.51 , offset.text = 0.05, barsize = 2,, color = "black")
# g <- g +geom_strip("FL0008_23_2PvD_Cluster_151_Region_LAM(1)","FL0005_23_1PvD_Cluster_151_Region_LAM(1)",label = "Domestic 20023", offset = 0.51 , offset.text = 0.05, barsize = 2,, color = "Blue")

# g <- p + geom_cladelab(node = 131, label  ="Internal Control", textcolor = "black", offset = .55, angle = -15)
# g <- g + geom_cladelab(node = 92, label  ="FL Domestic 2023", textcolor = "black", offset = 0.55, angle = -12)
# g <- g + geom_cladelab(node = 116, label  ="FL Domestic 2003", textcolor = "black", offset = 0.55, angle = -75)
# g <- g +geom_strip(17,17,label = "Control", offset = 0.51 , offset.text = 0.05, barsize = 2,, color = "black")



#that completes the manual label (if wanted)

#Add strip of white on the tree, can be useful for separating different levels of a annotation but not always necessary
# p <- p + geom_hilight(node = c(as.numeric(tip_labels$Tips)),fill = "white", extend = 0.05, alpha = 1, linewidth = 1)

#Another way to distinguish samples on the tree is with a shape
# g <- g + geom_tippoint(aes(subset = (node %in% AR_node)), size = 4, shape = 21, fill = "sienna3", color="black")
# g <- g + geom_tippoint(aes(subset = (node %in% FL2003)), size = 4, shape = 21, fill = "green3", color="black")
# g <- g + geom_tippoint(aes(subset = (node %in% FL2023)), size = 4, shape = 21, fill = "purple2", color="black")
# g <- g + geom_tippoint(aes(subset = (node %in% TX1_node)), size = 4, shape = 21, fill = "blue", color="black")
# g <- g + geom_tippoint(aes(subset = (node %in% TX3_node)), size = 4, shape = 21, fill = "gold", color="black")
# g <- g + geom_tippoint(aes(subset = (node %in% TX_8node)), size = 4, shape = 21, fill = "red", color="black")
# g <- g + geom_tippoint(aes(subset = (node %in% controls)), size = 4, shape = 24, fill = "black", color="black")
# g <- g + geom_tippoint(aes(subset = (node %in% PAmarko)), size = 4, shape = 21, fill = "red", color="black")
# g <- g + geom_tippoint(aes(subset = (node %in% GA06)), size = 4, shape = 21, fill = "purple", color="black")
# g <- g + geom_tippoint(aes(subset = (node %in% FLmarko)), size = 4, shape = 21, fill = "blue", color="black")

#This draws a line between the specified nodes, could be useful in some instances
# g <- g + geom_taxalink(taxa1 = 99, taxa2 =100, outward = T, color = "red", size = 1)
# g <- g + geom_taxalink(taxa1 = 260, taxa2 =c(261,266,267), outward = T, color = "blue", size = 1)
# g <- g + geom_taxalink(taxa1 = 263, taxa2 =c(264,265,233), outward = T, color = "gold", size = 1)
# g <- g + geom_taxalink(taxa1 = 260, taxa2 =261, outward = F, color = "blue", size = 1)
# g <- g + geom_taxalink(taxa1 = 94, taxa2 =30, outward = F, color = "mediumorchid3", size = 1)
# g <- g + geom_taxalink(taxa1 = 95, taxa2 =42, outward = F, color = "orange1", size = 1)
# g <- g + geom_taxalink(taxa1 = 96, taxa2 =51, outward = F, color = "sienna2", size = 1)
# g <- g + geom_taxalink(taxa1 = 97, taxa2 =52, outward = F, color = "violetred2", size = 1)




time_now <- format(Sys.time(), "%H%M")
date_now <- Sys.Date()
name_tree <- paste(date_now,"_", time_now, "_Pvivax_allStates_tree.pdf", sep="")
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



















