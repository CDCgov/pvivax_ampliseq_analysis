
print("Now computing the cutoff distance")

########################################################################################
# CAYETANENSIS - big loop below BIG LOOP MULTI-THREADED:
########################################################################################
print(number_of_clusters)
# number_of_clusters <- 5
GET_CUTOFF = mclapply(1:2, function (z) {
#CAYETANENSIS_get_cutoffs_under_first_peak = mclapply(1:10, function (z) {	

list_of_clusters <- list()
# 
for (m in 1:number_of_clusters){
	print("new Check cluster")
	print(m)

cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, k=m))))
print(cluster[,1])
#generate data frame showing cluster membership.

cluster_good <- cbind(rownames(cluster), cluster)
colnames(cluster_good) <- c("Seq_ID", "Cluster")
rownames(cluster_good) <- NULL
colnames(cluster_good) <- c("Seq_ID", "Cluster")
cluster_good <- cbind(cluster_good, paste0("Cluster_",cluster_good$Cluster,"_"))
colnames(cluster_good) <- c("Seq_ID","Cluster","Cluster_name")
cluster_good <- cluster_good[order(cluster_good$Cluster),]
print("numb clusters")
print(head(cluster_good, n = 20))
print(tail(cluster_good, n = 20))
# print(cluster_good)

## generate list of data frames, where each data frame contains members of the same cluster
number_of_clusters <- max(as.numeric(cluster$value))
list_of_clusters <- list()

for(i in 1:number_of_clusters){
	print("within list loop")
	print(i)
foobar <- cluster_good %>% filter(str_detect(Cluster_name, paste0("Cluster_",i,"_")))
	
	print(head(foobar))
	print(tail(foobar))
	list_of_clusters[[paste0("Cluster_",i,"_")]] <- foobar

      }	
   }

### randomly sample two specimens (or a number representing the smallest cluster size) from each cluster
### using previously created list of data frames - to generate a normalised set.

unbiased_list_of_specimens <- NULL

for(j in 1:number_of_clusters){
current_cluster <- paste0("Cluster_",j,"_")	
print(current_cluster)
current_data_frame <- list_of_clusters[[current_cluster]]
# random sampling of a number of specimens equal to the smallest cluster size, from each cluster.
random_ones <- sample_n(current_data_frame, smallest_cluster_size) 
unbiased_list_of_specimens <- rbind(unbiased_list_of_specimens, random_ones)
print(unbiased_list_of_specimens)
}


#d2 <- density(listed_no_self$value)
#plot(d2, lwd = 4, col = "black")
################################################################################

print("newest Clusters")
NEW_clusters <-  merge(unbiased_list_of_specimens, NEW_clusters, by.x = "Seq_ID", all.x=TRUE) 
# print(NEW_clusters)
NEW_clusters $Cluster <- NULL
NEW_clusters $Cluster_name <- NULL
# print(NEW_clusters)


## make a normalised patristic matrix
patristic_NEW_clusters_matrix <- matrix(nrow = length(NEW_clusters $Seq_ID), ncol = length(NEW_clusters $Seq_ID))
rownames(patristic_NEW_clusters_matrix) <-  NEW_clusters $Seq_ID
colnames(patristic_NEW_clusters_matrix) <- NEW_clusters $Seq_ID
# print(patristic_NEW_clusters_matrix)
for(p in rownames(patristic_NEW_clusters_matrix)){
	
	for(q in rownames(patristic_NEW_clusters_matrix)){
		
	patristic_NEW_clusters_matrix[p, q] <- PatristicDistMatrix[p, q]
		
	}
	
}
# print(patristic_NEW_clusters_matrix)

colnames(matrix) <- rownames(matrix)

## make a normalised raw genetic distance matrix

raw_normalized_matrix <- matrix(nrow = length(NEW_clusters $Seq_ID), ncol = length(NEW_clusters $Seq_ID))
rownames(raw_normalized_matrix) <-  NEW_clusters $Seq_ID
colnames(raw_normalized_matrix) <- NEW_clusters $Seq_ID

for(p in rownames(raw_normalized_matrix)){
	
	for(q in rownames(raw_normalized_matrix)){
		
	raw_normalized_matrix[p, q] <- matrix[p, q]
		
	}
	
}

print("raw matrix")
# print(raw_normalized_matrix)
#raw_matrix

raw_normalized_matrix <- as.dist(raw_normalized_matrix)
print("after as dist")
print("checking here")
raw_list_norm_matrix <- dist2list(raw_normalized_matrix)
print("After dist2list")
raw_list_norm_matrix_no_self2self <- filter(raw_list_norm_matrix, col != row) #remove self to self
print("after raw no self")

#sorted_raw_CAYETANENSIS_list_norm_matrix_no_self2self <- as.data.frame(raw_CAYETANENSIS_list_norm_matrix_no_self2self[order(raw_CAYETANENSIS_list_norm_matrix_no_self2self$value),])


#patristic_matrix
patristic_NEW_clusters_matrix <- as.dist(patristic_NEW_clusters_matrix)
patristic_list_norm_matrix <- dist2list(patristic_NEW_clusters_matrix)
patristic_norm_matrix_no_self2self <- filter(patristic_list_norm_matrix, col != row) #remove self to self
# print(patristic_norm_matrix_no_self2self)
# print("after no self")
# print(z/bootstraps)

D = sort(patristic_norm_matrix_no_self2self $value)[rank(raw_list_norm_matrix_no_self2self $value,ties.method="random")]
pat_plus_raw <- raw_list_norm_matrix_no_self2self
pat_plus_raw$D <- D

pat_plus_raw <- as.data.frame(pat_plus_raw)

sorted_pat_plus_raw <- pat_plus_raw[order(as.numeric(pat_plus_raw $value)),]


}, mc.cores= number_of_threads)

print("After bootstrap lloop")

########################################################################################
## FINAL CALCULATIONS
########################################################################################

foo <- GET_CUTOFF


for(z in 1:length(foo)){
	foo[[z]]$col <- NULL
	foo[[z]]$row <- NULL
	foo[[z]]$D <- NULL
	foo[[z]]$value <- as.numeric(foo[[z]]$value)
}

average_cutoff_under_first_peak_foo <- (Reduce("+", foo)/length(foo))




dens <- density(average_cutoff_under_first_peak_foo $value)
plot(dens, lwd = 4, col = "black") #+ title(main = TITLE) #+ abline(v=0.265, col="red", lty="dashed", lwd = 2)


max_first_step <- round(length(dens$x)*0.2)
first_peak_position <- which.max(dens$y[1:max_first_step])
last_peak_position <- which.max(dens$y)
dip_after_first_peak <- which.min(dens$y[first_peak_position:last_peak_position])

raw_threshold_for_PARNAS <- dens$x[(first_peak_position+dip_after_first_peak)]


write.table(raw_threshold_for_PARNAS, file = paste(date_today,".Threshold.csv", sep = ""), sep =",", quote = F, row.names = F)


#patristic_position <- (first_peak_position+dip_after_first_peak)

########################################################################################

#the part below does not do anything
keep_only_numbers <- GET_CUTOFF

for(z in 1:length(keep_only_numbers)){
	keep_only_numbers[[z]]$col <- NULL
	keep_only_numbers[[z]]$row <- NULL
	keep_only_numbers[[z]]$value <- NULL
	keep_only_numbers[[z]] <- as.numeric(keep_only_numbers[[z]]$D)
	#keep_only_numbers[[z]] <- as.numeric(keep_only_numbers[[z]]$value)
}

average_cutoff_under_first_peak <- (Reduce("+", keep_only_numbers)/length(keep_only_numbers))

dens_pat <- density(average_cutoff_under_first_peak)
#plot(dens_pat, lwd = 4, col = "black") 

#dens_pat$x[patristic_position]

########################################################################################










