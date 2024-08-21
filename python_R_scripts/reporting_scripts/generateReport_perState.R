library(knitr)
library(kableExtra)
library(dplyr)
library(readxl)

document_number = "Not Cleared"
version_number = " - Draft Document"

doc_control <- paste(document_number, version_number, sep="")

time_now <- format(Sys.time(), "%H%M")
date_now <- Sys.Date()
printDate <- paste("Date of Report Generation (yyyy/mm/dd): ", date_now)


args = commandArgs(trailingOnly = TRUE)


Abbreviations <- c("AF", "EAS", "ESEA", "LAM", "MSEA", "OCE", "WAS", "WSEA")
RegionNames <- c("Africa", 
"East Asia", 
"East Southeast Asia",
 "Central / South America", 
 "Malaysia Region Southeast Asia", 
 "Oceania", 
 "Western Asia",
 "Western Southeast Asia")
 geographicRegions <- data.frame(Abbreviations, RegionNames)
geoKable <- knitr::kable(geographicRegions, "html", padding = 40, line_sep = 2, align = "cc") %>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center", font_size = 12)



#Read in the latest metadata, clusters, geo prediction, and tree
fullMeta <- read_xlsx(args[1])
parnasClusters <- read.table(args[2], sep = "\t", header =F)
geographic <- read.table(args[3], sep ="\t", header =T)
stateID <- args[4]


#get the most recently generated tree for the state
treeList <- list.files(pattern = paste(stateID, "_Pvivax_ampliseq_tree.pdf", sep = ""))
treeFile <- c(treeList[-1])

#Some data cleaning
colnames(parnasClusters) <- c("Seq_ID", "Cluster")
parnasClusters$Cluster <- parnasClusters$Cluster +1
fullMeta$Travel <- gsub(" ", "", fullMeta$Travel)
colnames(fullMeta)[1] <- "Seq_ID"
fullMeta$Seq_ID <- gsub("-","_", fullMeta$Seq_ID)
fullMeta$Date_Collected <- as.Date(as.numeric(fullMeta$Date_Collected), origin = "1899-12-30")

#Extract the columns of interest from the metadata sheet

fullMeta <- fullMeta[,c(1,2,3,4,7,12,13)]

print(colnames(fullMeta))
print(colnames(parnasClusters))
merge1 <- merge(fullMeta, parnasClusters, by = "Seq_ID")
print(colnames(merge1))
print("after merge1")
print(colnames(geographic))
merge2 <- merge(merge1, geographic, by = "Seq_ID")
print("after merge2")

#more data cleaning so that the rows are in order and grouped together by cluster. Remove any duplicate rows (There are some floating around in different sheets)
merge2 <- merge2[order(merge2$Cluster),]
merge2 <- merge2[!duplicated(merge2),]

merge2 <- merge2[grepl(StateID, merge2$Seq_ID),]
# print(colnames(merge2))


#Make the rownames a new column so we can use the pack_rows command
rownames(merge2) <- 1:nrow(merge2)
merge2$name_of_rows <- rownames(merge2)


#make the kable table that will be in the rmarkdown report
kable_out <- knitr::kable(merge2, "html", padding = 40, line_sep = 2, align = "c") %>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "center", font_size = 12)

# print("after kable out")
# print(kable_out)
#Change the =9 if you update the number of columns in the final merge2 df
kable_out <- add_header_above(kable_out, c("Genetic clusters"= 10), italic = FALSE, font_size = 16, align = "left", line = FALSE, color = "#F5F5F5", background = "#9A9A9A", extra_css = "line-height: 20pt;")

print("second")
#loop throug the genetic clusters to find all specimens in each genetic cluster
list_of_genetic_clusters <- as.numeric(unique(merge2$Cluster))
for(l in list_of_genetic_clusters) {

	filtered_to_current_cluster <- filter(merge2, Cluster == l)
	specimens_to_look_at <- as.numeric(filtered_to_current_cluster$name_of_rows)

	#This is where the individual row numbers are extracted and then the range is used to pack rows of the same cluster into the same subtable
	left_value <- min(specimens_to_look_at)
	right_value <- max(specimens_to_look_at) 
	
	kable_out <- kable_out %>%
		pack_rows("Genetic cluster detected:", left_value, right_value, bold = TRUE, indent = FALSE,  label_row_css = "border-top: 2px solid; border-bottom: 2px solid; color:#704EA5")	
	
}
#Make column widths uniform (or can change if think it looks better)
kable_out <- column_spec(kable_out, 1, width = "7cm")
kable_out <- column_spec(kable_out, 2, width = "7cm")
kable_out <- column_spec(kable_out, 3, width = "7cm")
kable_out <- column_spec(kable_out, 4, width = "7cm")
kable_out <- column_spec(kable_out, 5, width = "7cm")
kable_out <- column_spec(kable_out, 6, width = "7cm")
kable_out <- column_spec(kable_out, 7, width = "7cm")
kable_out <- column_spec(kable_out, 8, width = "7cm")
kable_out <- column_spec(kable_out, 9, width = "7cm")
kable_out <- column_spec(kable_out, 10, width = "7cm")

#Now we're going to write the table to a csv file. This csv file will then be embedded in the rmarkdown report. I got rid of the row names columns because it isn't necessary in the csv file
merge2 <- merge2[,-9]
csvName <- paste(date_now, "_", time_now,"_", stateID, "_pvivax_clustering_report.csv", sep ="")

write.csv(merge2, file = csvName, quote = F, row.names =F)
file.copy(from = csvName, to = paste("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/P_vivax_fullAnalysis/python_R_scripts/reporting_scripts/",csvName, sep = ""))

#make the rmarkdown file
rmarkdown::render("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/P_vivax_fullAnalysis/python_R_scripts/reporting_scripts/Pvivax_markdown.Rmd", "html_document",  run_pandoc = TRUE)

htmlName <- paste(date_now, "_", time_now,"_", stateID, "_pvivax_clustering_report.html", sep = "")

file.copy(from = "/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/P_vivax_fullAnalysis/python_R_scripts/reporting_scripts/Pvivax_markdown.html", to = paste("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/P_vivax_fullAnalysis/reportFolder/",htmlName, sep = ""))
file.copy(from = csvName, to = paste("reportFolder/",csvName, sep = ""))
file.remove(csvName)