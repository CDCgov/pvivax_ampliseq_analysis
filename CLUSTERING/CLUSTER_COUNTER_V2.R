
#### THE FOLLOWING PART YOU REPEAT A THOUSAND TIMES SO THAT YOU CAN GET AN ACCURATE AVERAGE CUTOFF FROM MANY REPLICATES
#### otherwise you may get bouncy specimens again and the epis dont really like that. Paralellize it so it doesn't take a decade.

#number_of_species = as.numeric(readLines("SPECIES"))

## select the range of clusters to check for the maximum number of clusters that has at least 2 specimens in it (no less).
min_cluster_range <- 5
max_cluster_range <- 28

#bootstraps <- 1000 # number of times to run the code in order to calculate an average cutoff
bootstraps = as.numeric(readLines("BOOTSTRAPS"))

### first we need to select and import the most recent matrix -- DO WE STILL NEED THESE NEXT FEW LINES OF HASHED OUT CODE.
#matrix_location <- readLines("MATRIX_LOCATION")

# setwd(matrix_location)
#myFiles <- list.files(pattern="*csv", all.files = FALSE)
#myMatrix <- (sort(myFiles, decreasing = TRUE)[1])
#matrix <- as.matrix(read.table(myMatrix, sep = ",", row.names=1, header=TRUE))


# setwd("../CLUSTERING/TMP_REP/")
number_of_threads = as.numeric(readLines("THREADS"))



print(paste0("Now testing a range of cluster numbers from ", min_cluster_range, " to ", max_cluster_range, " to identify the maximum cluster number where no singletons exist in their own cluster."))




## Cluster_distance_matrix
Ensemble_x <- as.hclust(agnes(x=matrix, diss = TRUE, stand = TRUE, metric = "manhattan", method = "ward"))

#lets try nj
#Ensemble_x <- as.phylo(nj(matrix))

print("after ensemble")

pick_min_cluster_number <- c()

###just_for strongy data we can use max_cluster_range <- 200

for (k in min_cluster_range:max_cluster_range){
	print(k)
#for (k in 50){
print("start loop")
cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, k=k))))

#cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, k=48))))
#generate data frame showing cluster membership.

cluster_good <- cbind(rownames(cluster), cluster)
colnames(cluster_good) <- c("Seq_ID", "Cluster")
rownames(cluster_good) <- NULL
colnames(cluster_good) <- c("Seq_ID", "Cluster")
cluster_good <- cbind(cluster_good, paste0("Cluster_",cluster_good$Cluster,"_"))
colnames(cluster_good) <- c("Seq_ID","Cluster","Cluster_name")
cluster_good <- cluster_good[order(cluster_good$Cluster),]
print("after cluster good")
## generate list of data frames, where each data frame contains members of the same cluster
number_of_clusters <- max(as.numeric(cluster$value))
list_of_clusters <- list()

for(i in 1:number_of_clusters){
	print("within second loop")
	print(i)
	
foobar <- cluster_good %>% filter(str_detect(Cluster_name, paste0("Cluster_",i,"_")))
	
	if(length(foobar$Seq_ID) == 1){
			
		pick_min_cluster_number <- cbind(pick_min_cluster_number, print(k-1))
		
		break
	} else {
	
	list_of_clusters[[paste0("Cluster_",i,"_")]] <- foobar
	
### Now we need an "if" loop that checks if any given cluster has only one specimens. if a cluster has only 1 specimens the loop will re-cluster
### the matrix at the previous cluster number, and then it will save the list of clusters variable for the following steps.
      }	
   }
}

print("appropriate number of clusters")
number_of_clusters <- min(pick_min_cluster_number[,1])
print(number_of_clusters)



# now we need a loop to determine the number of specimens in the smallest cluster
# just in case we can select more than 2 specimens in the next step. If we can, then we should for the sake of diversity.
# This value will usually be 2 however.

cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, k= number_of_clusters))))
size_of_each_cluster <- vector()

for(p in 1:number_of_clusters){
print("check current clucster")
print(p)
current_cluster <- filter(cluster, value %in% p)
print(current_cluster)

size_of_each_cluster <- append(size_of_each_cluster,length(current_cluster$value))

}

smallest_cluster_size <- min(size_of_each_cluster)



NEW_clusters <- as.data.frame(melt(factor(cutree(Ensemble_x, k= number_of_clusters)))) # this is from the original large matrix - the Ensemble_x variable.
NEW_clusters <- cbind(NEW_clusters, c(rownames(NEW_clusters)))
colnames(NEW_clusters) <- c("value", "Seq_ID") ## this will be used later.
print("new Clusters")
print(NEW_clusters)

#extract patristic distances from your tree Ensemble_x
PatristicDistMatrix <-as.matrix(cophenetic(Ensemble_x))




