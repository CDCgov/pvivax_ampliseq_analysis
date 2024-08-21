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


#/Users/joelbarratt/Documents/NEW_MODULE_3_WORK/ensemble_matrices

matrix <- as.matrix(read.table("2022-10-24Joel_haplotype_sheet_H_matrix.csv", sep = ",", row.names=1, header=TRUE))
Ensemble_x <- as.hclust(agnes(x=matrix, diss = TRUE, stand = TRUE, method = "ward"))



#Ensemble_x <- nj(matrix)

Ensemble_y <- as.phylo(Ensemble_x)

#ape::write.tree(Ensemble_y, file='cyclo_tree_newick.txt') 
#ape::write.nexus(Ensemble_y, file='cyclo_tree_nexus.nex')


tip_labels <- as.data.frame(Ensemble_y$tip.label)
colnames(tip_labels) <- "Seq_ID"
tip_labels$Tip <- rownames(tip_labels)



new_tips <- merge(cay_ash_clusters, tip_labels, by=c("Seq_ID"), all.x=TRUE)
#new_tips <- merge(dyanmic_cut, tip_labels, by=c("Seq_ID"), all.x=TRUE)
#two_species <- as.data.frame(melt(factor(cutree(Ensemble_x, k= 2))))
# two_species_fix <- cbind(rownames(two_species), two_species)
#rownames(two_species_fix) <- NULL
#colnames(two_species_fix) <- c("Seq_ID", "Cluster")
#two_species_fix_2 <- cbind(two_species_fix, paste0("Cluster_", two_species_fix$Cluster))
#colnames(two_species_fix_2) <- c("Seq_ID", "Cluster", "Species")
#new_tips <- merge(two_species_fix_2, tip_labels, by=c("Seq_ID"), all.x=TRUE)

sorted_new_tips <- new_tips[order(as.numeric(new_tips$Tip)),]
sorted_new_tips$NEW_NAME <- paste0(sorted_new_tips$Seq_ID," ", sorted_new_tips$Species)
Ensemble_y$tip.label <- sorted_new_tips$NEW_NAME

cols <- c("black","gray60", "black")
Ensemble_y <- groupClade(Ensemble_y, .node=c(2114, 2113)) ###These will be the 2 clusters 2114 is ashfordi, 2113 is cayetanensis




#outgroup <- "C_ChHenan"
 #rooted_tree <- root(Ensemble_y, resolve.root = TRUE, outgroup=outgroup)

p <- ggtree(Ensemble_y, size = 1.4, layout = "circular", aes(color = group)) + scale_color_manual(values = cols)

#p <- ggtree(Ensemble_y, size = 0.9, layout = "circular") + geom_tiplab(color = "black", size = 0.2, offset = 0.05)
#p <- ggtree(rooted_tree, size = 0.9, layout = "circular") + geom_tiplab(color = "black", size = 0.2, offset = 0.005)




#x <- p +  geom_tiplab2(color = "black", size = 0.6, offset = 0.05) #+ ####dont forget to add this plus symbol back if you want a node label.
x <- p #+  geom_tiplab(color = "black", size = 0.2, offset = 0.5) #+ ####dont forge
#x + geom_text2(aes(subset=!isTip, label=node), hjust=-0.05, size = 2, color = "red") ###

#length(unique(cay_ash_clusters$Species))






#palette <- distinctColorPalette(length(unique(cay_ash_clusters$Species)))
palette <- distinctColorPalette(length(unique(two_species_fix_2 $Species)))

unique_genotypes <- unique(cay_ash_clusters$Species)

#unique_genotypes <- 2
	
#TREE_tips = mclapply(1:length(unique(cay_ash_clusters$Species)), function (n) {
	TREE_tips = mclapply(1:length(unique(two_species_fix_2 $Species)), function (n) {

#these_tips <- filter(sorted_new_tips, Species == unique(cay_ash_clusters$Species)[n])
these_tips <- dplyr::filter(sorted_new_tips, Species == unique(two_species_fix_2 $Species)[n])

these_tips <- as.numeric(these_tips$Tip)

	}, mc.cores= number_of_threads)	
	
	
for(m in 1:length(TREE_tips)){
	x <- x + geom_hilight(node=c(TREE_tips[[m]]), fill=palette[m], extend = 0.45, alpha = 1)
}








##dynamic cuttree
	
#result_Dynamic_cut <- cutreeDynamic(Ensemble_x,  method = "tree") # -- USES ONLY THE TREE

#result_Dynamic_cut <- cutreeHybrid(Ensemble_x,  matrix) # USING HYBRID METHOD

#max(cutreeDynamic(Ensemble_x,  method = "tree"))
	
	
#dyanmic_cut <- as.data.frame(cbind(tip_labels, result_Dynamic_cut$labels))
#dyanmic_cut <- as.data.frame(cbind(tip_labels, result_Dynamic_cut))
#dyanmic_cut$Tip <- NULL	

#colnames(dyanmic_cut) <- c("Seq_ID", "Cluster")


#write.table(dyanmic_cut,  "DYNAMIC_CUT_HYBRID_TREE_ONLY.txt", col.names=T, quote=FALSE, sep="\t",row.names=FALSE)



#new_tips_DY <- merge(dyanmic_cut, tip_labels, by=c("Seq_ID"), all.x=TRUE)
#sorted_new_tips_DY <- new_tips_DY[order(as.numeric(new_tips_DY$Tip)),]
#sorted_new_tips_DY $NEW_NAME <- paste0(sorted_new_tips_DY $Seq_ID,"_Cluster_", sorted_new_tips_DY $Cluster,"_")
#Ensemble_y$tip.label <- sorted_new_tips_DY$NEW_NAME
#p <- ggtree(Ensemble_y, size = 1.4, layout = "circular") #+ scale_color_manual(values = cols)
#x <- p +  geom_tiplab(color = "black", size = 0.2, offset = 0.05) #+ ####dont forge
#length(unique(new_tips_DY$Cluster))

#palette <- distinctColorPalette(length(unique(new_tips_DY$Cluster)))
#unique_genotypes <- unique(new_tips_DY$Cluster)
	
	

#sorted_new_tips_DY$Cluster_next <- paste0("clus",sorted_new_tips_DY$Cluster,"_")

#TREE_tips = mclapply(1:26, function (n) {

#these_tips <- filter(sorted_new_tips_DY, Cluster_next == paste0("clus",n,"_"))
#these_tips <- as.numeric(these_tips$Tip)

#	}, mc.cores= number_of_threads)	
		
	
#for(m in 1:length(TREE_tips)){
#	x <- x + geom_hilight(node=c(TREE_tips[[m]]), fill=palette[m], extend = 0.45, alpha = 1)
#}




#tree cluster trees

number_of_threads <- 10

thresh <- 0.786
which_threshold <- (paste0("treecluster_", thresh,".out.txt"))

treeClusters <- read.table(which_threshold, header=T, sep="\t")

colnames(treeClusters) <- c("Seq_ID", "Cluster")


#/Users/joelbarratt/Documents/NEW_MODULE_3_WORK/ensemble_matrices

matrix <- as.matrix(read.table("/Users/joelbarratt/Documents/NEW_MODULE_3_WORK/ensemble_matrices/2022-10-24Joel_haplotype_sheet_H_matrix.csv", sep = ",", row.names=1, header=TRUE))
Ensemble_x <- as.hclust(agnes(x=matrix, diss = TRUE, stand = TRUE, metric = "manhattan", method = "ward"))
Ensemble_y <- as.phylo(Ensemble_x)




tip_labels <- as.data.frame(Ensemble_y$tip.label)
colnames(tip_labels) <- "Seq_ID"
tip_labels$Tip <- rownames(tip_labels)
new_tips <- merge(treeClusters, tip_labels, by=c("Seq_ID"), all.x=TRUE)
#new_tips <- merge(dyanmic_cut, tip_labels, by=c("Seq_ID"), all.x=TRUE)



sorted_new_tips <- new_tips[order(as.numeric(new_tips$Tip)),]
sorted_new_tips$NEW_NAME <- paste0(sorted_new_tips$Seq_ID," ", sorted_new_tips$Cluster)
Ensemble_y$tip.label <- sorted_new_tips$NEW_NAME

cols <- c("black","gray60", "black")
Ensemble_y <- groupClade(Ensemble_y, .node=c(2114, 2113)) ###These will be the 2 clusters 2114 is ashfordi, 2113 is cayetanensis

p <- ggtree(Ensemble_y, size = 1.4, layout = "circular", aes(color = group)) + scale_color_manual(values = cols)

#p <- ggtree(Ensemble_y, size = 1.4, layout = "circular") + geom_tiplab(color = "black", size = 0.4, offset = 0.05)


x <- p #+  geom_tiplab(color = "black", size = 0.2, offset = 0.5) #+ ####dont forge







palette <- distinctColorPalette(length(unique(treeClusters$Cluster)))
unique_genotypes <- unique(treeClusters$Cluster)
	

TREE_tips = mclapply(1:length(unique(treeClusters$Cluster)), function (n) {

these_tips <- dplyr::filter(sorted_new_tips, Cluster == unique(treeClusters$Cluster)[n])
these_tips <- as.numeric(these_tips$Tip)

	}, mc.cores= number_of_threads)	
	
	
for(m in 1:length(TREE_tips)){
	x <- x + geom_hilight(node=c(TREE_tips[[m]]), fill=palette[m], extend = 0.45, alpha = 1)
}



























