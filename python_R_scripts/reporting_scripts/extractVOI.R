library(knitr)
library(kableExtra)
library(dplyr)

document_number = "Not Cleared"
version_number = " - Draft Document"

doc_control <- paste(document_number, version_number, sep="")

time_now <- format(Sys.time(), "%H%M")
date_now <- Sys.Date()
printDate <- paste("Date of Report Generation (yyyy/mm/dd): ", date_now)
stateID <- "All_States"

# treeFile <- pdf("Rplots.pdf")



args = commandArgs(trailingOnly = TRUE)

# inVOI <- read.csv("Reportable_snps.csv", header = T, sep =",")
inVOI <- read.csv(args[2], header = T, sep =",")

# inGeoPrediction <- read.csv("../predictedOut/2023-10-18_1218_gatkRegion_prediction_withMetadata.txt", header = T, sep ="\t")
inGeoPrediction <- read.csv(args[3], header = T, sep ="\t")

# treeFile <- c("Rplots.pdf")

sampleList <- unique(inVOI$Sample_name) 

markers <- unique(inVOI$CHROM)


inVOI$Annotation <- gsub("missense_variant", "NonSynonymous", inVOI$Annotation)
inVOI$Annotation <- gsub("synonymous_variant", "Synonymous", inVOI$Annotation)
testDF <- as.data.frame(matrix(ncol = 5, nrow = length(sampleList)))
count = 0
for(i in sampleList){
	count <- count +1
	myTab <- inVOI[grep(i, inVOI$Sample_name),]
	reportable <- myTab[!grepl("Wildtype", myTab$Type),]
	reportable <- reportable[!grepl("No coverage", myTab$Type),]
	reportable <- reportable[grepl("WHO", myTab$SNP_REPORT),]
	if(length(reportable$Sample_name)>0){
	
	print(i)
	# print(head(myTab))
	# print(reportable)
	marker1 <- reportable[grep("Pvdhps", reportable$CHROM),c(1,2,13,5,6,15,10)]
	if(length(marker1$Sample_name)>0){
		final_marker1 <- c()
		marker1 <- marker1[order(marker1$AVG_COV, decreasing = T),]
		marker1$pasted <- paste(marker1$VOI, " [",  marker1$AVG_COV, ",", marker1$AVG_VAF, "]", sep ="")
		# print(marker1)
		for(j in 1:length(marker1$Sample_name)){
			final_marker1 <- c(final_marker1, marker1$pasted[j])
		}	
		final_marker1 <- paste(final_marker1, collapse ="; ")
		# print(paste(final_marker1, collapse ="; "))
	}

	else{
		final_marker1 <- "None"
		print("no variant marker1")
	}

	marker2 <- reportable[grep("Pvmdr1", reportable$CHROM), c(1,2,13,5,6,15,10)]
	if(length(marker2$Sample_name)>0){
		final_marker2 <- c()
		marker2 <- marker2[order(marker2$AVG_COV, decreasing = T),]
		marker2$pasted <- paste(marker2$VOI, " [", marker2$AVG_COV, ",", marker2$AVG_VAF, "]", sep ="")
		# print(marker2)
		for(j in 1:length(marker2$Sample_name)){
			final_marker2 <- c(final_marker2, marker2$pasted[j])
		}	
		final_marker2 <- paste(final_marker2, collapse ="; ")
		# print(paste(final_marker2, collapse ="; "))
	}

	else{
		final_marker2 <- "None"
		print("no variant marker2")
	}

	marker3 <- reportable[grep("Pvdhfr", reportable$CHROM), c(1,2,13,5,6,15,10)]
	if(length(marker3$Sample_name)>0){
		final_marker3 <- c()
		marker3 <- marker3[order(marker3$AVG_COV, decreasing = T),]
		marker3$pasted <- paste(marker3$VOI, " [",marker3$AVG_COV, ",", marker3$AVG_VAF, "]", sep ="")
		# print(marker3)
		for(j in 1:length(marker3$Sample_name)){
			final_marker3 <- c(final_marker3, marker3$pasted[j])
		}	
		final_marker3 <- paste(final_marker3, collapse ="; ")
		# print(paste(final_marker3, collapse ="; "))
	}

	else{
		final_marker3 <- "None"
		print("no variant marker3")
	}

	marker4 <- reportable[grep("Pvcrt", reportable$CHROM), c(1,2,13,5,6,15,10)]
	if(length(marker4$Sample_name)>0){
		final_marker4 <- c()
		marker4 <- marker4[order(marker4$AVG_COV, decreasing = T),]
		marker4$pasted <- paste(marker4$VOI, " [",  marker4$AVG_COV, ",", marker4$AVG_VAF, "]", sep ="")
		# print(marker4)
		for(j in 1:length(marker4$Sample_name)){
			final_marker4 <- c(final_marker4, marker4$pasted[j])
		}	
		final_marker4 <- paste(final_marker4, collapse ="; ")
		# print(paste(final_marker4, collapse ="; "))
	}

	else{
		final_marker4 <- "None"
		print("no variant marker4")
	}
	}
		else{
			print("no reportable")
		}
	testDF[count,1] <- i
	testDF[count,2] <- final_marker1
	testDF[count,3] <- final_marker2
	testDF[count,4] <- final_marker3
	testDF[count,5] <- final_marker4
	}
	# print(marker1)

inVOI$uniqueMutations <- paste(inVOI$CHROM, "_", inVOI$VOI, sep= "")
inVOI_nonWildtype <-  inVOI[!grepl("WHO", inVOI$SNP_REPORT),]
# print(table(inVOI_nonWildtype$uniqueMutations))

marker1_total <- inVOI_nonWildtype[grep("Pvdhps", inVOI_nonWildtype$CHROM),]
if(length(marker1_total$VOI > 0)){
marker1DF <- as.data.frame(table(marker1_total$VOI))
marker1DF[marker1DF ==0] <- NA
marker1DF_final <- marker1DF[complete.cases(marker1DF),]
marker1DF_final <- marker1DF_final[order(marker1DF_final$Freq, decreasing=T),]

rownames(marker1DF_final) <- marker1DF_final[,1]
marker1DF_final <- as.data.frame(marker1DF_final)
}else{
	marker1DF_final <- as.data.frame(matrix(ncol = 2))
}
colnames(marker1DF_final) <- c("Pvdhps_VOI", "Pvdhps_Frequency")
# print(marker1DF_final)

marker2_total <- inVOI_nonWildtype[grep("Pvmdr1", inVOI_nonWildtype$CHROM),]
if(length(marker2_total$VOI > 0)){
marker2DF <- as.data.frame(table(marker2_total$VOI))
marker2DF[marker2DF ==0] <- NA
marker2DF_final <- marker2DF[complete.cases(marker2DF),]
marker2DF_final <- marker2DF_final[order(marker2DF_final$Freq, decreasing=T),]

}else{
	marker2DF_final <- as.data.frame(matrix(ncol = 2))
}
colnames(marker2DF_final) <- c("Pvmdr1_VOI", "Pvmdr1_Frequency")
# print(marker2DF_final)

marker3_total <- inVOI_nonWildtype[grep("Pvdhfr", inVOI_nonWildtype$CHROM),]
if(length(marker3_total$VOI > 0)){
marker3DF <- as.data.frame(table(marker3_total$VOI))
marker3DF[marker3DF ==0] <- NA
marker3DF_final <- marker3DF[complete.cases(marker3DF),]
marker3DF_final <- marker3DF_final[order(marker3DF_final$Freq, decreasing=T),]

}else{
	marker3DF_final <- as.data.frame(matrix(ncol = 2))
}
colnames(marker3DF_final) <- c("Pvdhfr_VOI", "Pvdhfr_Frequency")
# print(marker3DF_final)

marker4_total <- inVOI_nonWildtype[grep("Pvcrt", inVOI_nonWildtype$CHROM),]
if(length(marker4_total$VOI > 0)){
marker4DF <- as.data.frame(table(marker4_total$VOI))
marker4DF[marker4DF ==0] <- NA
marker4DF_final <- marker4DF[complete.cases(marker4DF),]
marker4DF_final <- marker4DF_final[order(marker4DF_final$Freq, decreasing=T),]

}else{
	marker4DF_final <- as.data.frame(matrix(ncol = 2))
}
colnames(marker4DF_final) <- c("Pvcrt_VOI", "Pvcrt_VOI_Frequency")
# print(marker4DF_final)

# markerFreqFinal <- as.data.frame(as.matrix(nrow = 8))
# markerFreqFinal$
# markerFreqFinal <- cbind(marker1DF_final, marker2DF_final, marker3DF_final, marker4DF_final)
# print("freqFinal")
# print(markerFreqFinal)
inGeoTmp <- inGeoPrediction[,c(1,2,6,9,11)]
inGeoTmp$Miss <- (1 - inGeoTmp$Miss/113) * 100
colnames(inGeoTmp) <- c("CDC_Sample_ID", "Region Prediction (Probability)", "Perecent_Sites_Called", "Travel History", "Other ID")

colnames(testDF) <- c("Sample_ID","Pvdhps VOI","Pvmdr1 VOI", "Pvdhfr VOI", "Pvcrt VOI")
testDF2 <- merge(inGeoTmp, testDF, by.y = "Sample_ID", by.x = "CDC_Sample_ID", all.y = T)
# perSample_VOI <- testDF 
print(colnames(testDF2))
print("before kable")
perSample_VOI <- testDF2 %>%

# 			#create kable and align in center (not sure if align center is actually needed)
		  	kable(escape = F, align=rep('c', 5)) %>%

# 		  	#full width and boostrap options mean table doesn't take up full page and each cell has a border. I think it looks better this way
		  	kable_styling(position = "center", full_width = FALSE, bootstrap_options = "bordered") %>%

# 		  	#Color the rows depending on the validation decision. Could be changed based on how it looks
# 		  	row_spec(which(oneTable$Validation_Decision =="FAIL"), background = "#F88D76") %>%
# 		  	row_spec(which(oneTable$Validation_Decision =="CONDITIONAL_PASS"), background = "#76BAF8") %>%
# 		  	row_spec(which(oneTable$Validation_Decision =="COMPLETE_PASS"), background = "#F1F3F5") %>%

# 		  	#align everything in the center of the cells
		  	row_spec(1:length(sampleList), extra_css = "vertical-align:middle;") %>%
			
		  	column_spec(1:8,extra_css = "horizontal-align:middle") %>%

# 		  	#Give the Seq ID column a special dark grey color so it is easy to pick out
		    column_spec(1, bold = TRUE, border_right = TRUE, color = "black", background = "lightgrey") %>%
			column_spec(column = c(3,5,7,9), bold = TRUE, border_right = TRUE, color = "black", background = "#FFF5EE") %>%
			column_spec(column = c(2,4,6,8), bold = TRUE, border_right = TRUE, color = "black", background = "#FFFAFA") %>%

			# column_spec(column = c, bold = TRUE, border_right = TRUE, color = "black", background = "#FFF0F5")

			row_spec(which(testDF2$Perecent_Sites_Called < 50), background = "#F88D76")

print("after voi kable")

print(colnames(marker1DF_final))
kable1 <- kable(marker1DF_final, escape = F, align = rep('c', 5), row.names = F) %>% 
			kable_styling(position = "center", full_width = FALSE, bootstrap_options = "bordered") %>%
			row_spec(1:length(row.names(marker1DF_final)), extra_css = "vertical-align:middle;") %>%
		  	column_spec(1:2,extra_css = "horizontal-align:middle") %>%
	 		column_spec(1, bold = TRUE, border_right = TRUE, color = "black", background = "#8B8B7A") %>%
			column_spec(2, bold = TRUE, border_right = TRUE, color = "black", background = "#FFFFE0")

kable2 <- kable(marker2DF_final, escape = F, align = rep('c', 5), row.names =F) %>% 
			kable_styling(position = "center", full_width = FALSE, bootstrap_options = "bordered") %>%
			row_spec(1:length(row.names(marker2DF_final)), extra_css = "vertical-align:middle;") %>%
		  	column_spec(1:2,extra_css = "horizontal-align:middle") %>%
			column_spec(1, bold = TRUE, border_right = TRUE, color = "black", background = "#CDC5BF")%>%
			column_spec(2, bold = TRUE, border_right = TRUE, color = "black", background = "#FFF5EE")

kable3 <- kable(marker3DF_final, escape = F, align = rep('c', 5), row.names =F) %>% 
			kable_styling(position = "center", full_width = FALSE, bootstrap_options = "bordered") %>%
			row_spec(1:length(row.names(marker3DF_final)), extra_css = "vertical-align:middle;") %>%
		  	column_spec(1:2,extra_css = "horizontal-align:middle") %>%
			column_spec(1, bold = TRUE, border_right = TRUE, color = "black", background = "#9FB6CD")%>%
			column_spec(2, bold = TRUE, border_right = TRUE, color = "black", background = "#C6E2FF")

kable4 <- kable(marker4DF_final, escape = F, align = rep('c', 5), row.names =F) %>% 
			kable_styling(position = "center", full_width = FALSE, bootstrap_options = "bordered") %>%
			row_spec(1:length(row.names(marker4DF_final)), extra_css = "vertical-align:middle;") %>%
		  	column_spec(1:2,extra_css = "horizontal-align:middle") %>%
			column_spec(1, bold = TRUE, border_right = TRUE, color = "black", background = "#CDB38B")%>%
			column_spec(2, bold = TRUE, border_right = TRUE, color = "black", background = "#FFDEAD")
# kable3
# kable4
print("before combine")
# finalOutput <- knitr::kable(list(marker1_kable, marker2_kable, marker3_kable, marker4_kable))
print("final output")
# print(finalOutput)
rmarkdown::render(args[1])

rename_html <- paste0(date_now,"_", time_now, "_", stateID,"_Pvivax_MaRS_genotyping_report.html")
file.copy(from="Pvivax_markdown_VOI.html", to= rename_html, overwrite = TRUE)
# #This text goes in the header of the report


# #Get the latest validation files, will need to have wrapper script that pulls only the latest files. Something to grab each states most recent file
# myFiles <- list.files(pattern = "*cyclosporaValidations.csv")

# read.table

# #Go through each states validation csv file and create an html report
# for(i in 1:length(myFiles)){

# 	#most partners have a two letter abbreviation, internal testing (CDC) and canada (CAN) have three letters. Get the correct name for the state ID in either case
# 	stateID <- substr(myFiles[i], 12,14)
# 	if(stateID == "CDC" || stateID == "CAN"){
# 		stateID <- stateID
# 	} else{
# 		stateID <- substr(myFiles[i], 12,13)
# 	}


# 	#replace a extra characters included in the python output. Could be done in python script but it works fine here
# 	oneTable <- read.csv(myFiles[i])

# 	#The files will be empty if no validations have been submitted yet, so only do the full processing for partners that have submitted validations
# 	if(length(row.names(oneTable))>0){
# 		oneTable[] <- lapply(oneTable, gsub, pattern="\\[", replacement='')
# 		oneTable[] <- lapply(oneTable, gsub, pattern="\\]", replacement='')
# 		oneTable[] <- lapply(oneTable, gsub, pattern="\\{", replacement='')
# 		oneTable[] <- lapply(oneTable, gsub, pattern="\\}", replacement='')
# 		oneTable[] <- lapply(oneTable, gsub, pattern="_}", replacement='')
# 		oneTable[] <- lapply(oneTable, gsub, pattern="\'Mt_C\'", replacement='Mt_Junction')
# 		oneTable[] <- lapply(oneTable, gsub, pattern="'", replacement='')
# 		oneTable[] <- lapply(oneTable, gsub, pattern="_$", replacement='')
# 		oneTable[] <- lapply(oneTable, gsub, pattern="_,", replacement=',')
		
# 		# oneTable$MissingMarkers	<- gsub("Mt_C", "Mt_Junction",oneTable$MissingMarkers)

# 		val1 <- oneTable[grep("001[12345]_", oneTable$Seq_ID),2]
# 		val2 <- oneTable[grep("002[12345]_", oneTable$Seq_ID),2]
# 		val3 <- oneTable[grep("003[12345]_", oneTable$Seq_ID),2]
# 		val4 <- oneTable[grep("004[12345]_", oneTable$Seq_ID),2]
# 		val5 <- oneTable[grep("005[12345]_", oneTable$Seq_ID),2]
# 		val6 <- oneTable[grep("006[12345]_", oneTable$Seq_ID),2]
# 		val7 <- oneTable[grep("007[12345]_", oneTable$Seq_ID),2]
# 		val8 <- oneTable[grep("008[12345]_", oneTable$Seq_ID),2]
# 		val9 <- oneTable[grep("009[12345]_", oneTable$Seq_ID),2]
# 		val10 <- oneTable[grep("010[12345]_", oneTable$Seq_ID),2]
# 		val11 <- oneTable[grep("011[12345]_", oneTable$Seq_ID),2]
# 		val12 <- oneTable[grep("012[12345]_", oneTable$Seq_ID),2]
# # 
# 		# print(oneTable[grep("001[12345]_", oneTable$Seq_ID),])
# 		# print(oneTable[grep("011[12345]_", oneTable$Seq_ID),])
# 		# print(oneTable)
# 		# resultVec <- c()
# 		# testTab <- dataFrame()
# 		# for(i in 1:length(row.names(oneTable))){
# 			# result <- oneTable[i,2]
# 			# resultVec <- c(result, resultVec)
# 			# name <- paste("val", i, sep = "")

# 		# 
# 		# }
# 		# allVals <- c(val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12)
# 		allVals <- data.frame("Result" = 12)
# 		allVals[1,1] <- paste(val1, collapse = "_")
# 		allVals[2,1] <- paste(val2, collapse = "_")
# 		allVals[3,1] <- paste(val3, collapse = "_")
# 		allVals[4,1] <- paste(val4, collapse = "_") 
# 		allVals[5,1] <- paste(val5, collapse = "_")
# 		allVals[6,1] <- paste(val6, collapse = "_")
# 		allVals[7,1] <- paste(val7, collapse = "_")
# 		allVals[8,1] <- paste(val8, collapse = "_")
# 		allVals[9,1] <- paste(val9, collapse = "_")
# 		allVals[10,1] <- paste(val10, collapse = "_")
# 		allVals[11,1] <- paste(val11, collapse = "_")
# 		allVals[12,1] <- paste(val12, collapse = "_")
# 		print(stateID)
# 		print(allVals)
# 		allVals2 <- allVals[!apply(allVals == "", 1, all), ]
# 		print(allVals2)
# 		# print(oneTable[1,2])
# 		# print(stateID)
# 		# print(resultVec)
# 		# print(name)
# 		val1 <- allVals2[1]
# 		val2 <- allVals2[2]
# 		val3 <- allVals2[3]
# 		val4 <- allVals2[4]
# 		val5 <- allVals2[5]
# 		val6 <- allVals2[6]

# 		# print(oneTable)
# 		# print(val1)
# 		# print(val2)
# 		# print(val3)
# 		# print(val4)
# 		# print(val5)
# 		# print(val6)
		
# 		numberValidations <- length(row.names(oneTable))

# 		if(numberValidations < 6){
# 			output <- c("Not all validation specimens have been submitted to CDC. Please submit all 6 validation specimens to CDC. The results below apply to the validations you have submited; however, you will not be cleared to process clinical specimens until all 6 validations specimens have a PASS designation.")
# 		}	
# 		#Determine the status for the lab's validation specimens. If all pass validation, inform the lab. The 'output' variable is the text that goes into the html report under the 'Actions to Take' section.
# 		else if("COMPLETE_PASS" %in% val1 && "COMPLETE_PASS" %in% val2 && "COMPLETE_PASS" %in% val3 && "COMPLETE_PASS" %in% val4 && "COMPLETE_PASS" %in% val5 && "COMPLETE_PASS" %in% val6){
# 			output <- c("All validation specimens PASS. The lab may proceed with sequencing clinical specimens and submitting data to CDC")
		
# 		#Else, look at the most recent specimen for each validation and see if it says Fail. If any of the validations say fail, then inform the lab.
# 		} else if("FAIL" %in% tail(val1, n = 1) || "FAIL" %in% tail(val2, n =1) || "FAIL" %in% tail(val3, n =1) || "FAIL" %in% tail(val4, n =1) || "FAIL" %in% tail(val5, n =1) || "FAIL" %in% tail(val6, n =1)){
# 			output <- c("At least one validation specimen is classified as FAIL. Please refer to the Results Table to see the status of the validation specimens. CDC will not process any clinical specimen sequencing data until all validation specimens are either COMPLETE_PASS or CONDITIONAL_PASS.")
		
# 		#If not all specimens are PASS and none of the specimens are considered FAIL, then one or more must be a CONDITIONAL_PASS specimen. Inform the lab of this
# 		} 
# 		else{
# 			output <- c("At least one validation specimens is classified as CONDITIONAL_PASS. This means at least one validation specimen is missing 1 or more of the 8 genotyping markers but this vaildation specimen still belongs to the same genetic cluster as the original specimen. CDC will accept clinical specimen sequencing data generated by the lab. CDC recommends re-processing the validation specimen(s) classified as CONDITIONAL_PASS to see if all 8 markers can be amplified for each vaildation specimen; however, re-processing is not required.")

# 		}

# 		#Parse through the table to create the stateVD table. This table will go straight into the final html report.
# 		stateVD <- oneTable %>%

# 			#create kable and align in center (not sure if align center is actually needed)
# 		  	kable(escape = F, align=rep('c', 5)) %>%

# 		  	#full width and boostrap options mean table doesn't take up full page and each cell has a border. I think it looks better this way
# 		  	kable_styling(position = "center", full_width = FALSE, bootstrap_options = "bordered") %>%

# 		  	#Color the rows depending on the validation decision. Could be changed based on how it looks
# 		  	row_spec(which(oneTable$Validation_Decision =="FAIL"), background = "#F88D76") %>%
# 		  	row_spec(which(oneTable$Validation_Decision =="CONDITIONAL_PASS"), background = "#76BAF8") %>%
# 		  	row_spec(which(oneTable$Validation_Decision =="COMPLETE_PASS"), background = "#F1F3F5") %>%

# 		  	#align everything in the center of the cells
# 		  	row_spec(1:numberValidations, extra_css = "vertical-align:middle;") %>%
# 		  	column_spec(1:7,extra_css = "horizontal-align:middle") %>%

# 		  	#Give the Seq ID column a special dark grey color so it is easy to pick out
# 		    column_spec(1, bold = TRUE, border_right = TRUE, color = "black", background = "lightgrey")


# 		#run the rmarkdown script to print html to file  
# 		rmarkdown::render(args[1])
		
# 	#kableextra commands get angry if the kable is empty (when a partner has yet to submit validation specimens so their reports are empty). So create placeholder text for the 'output' and 'stateVD' variables
# 	} else{
# 		output <- c("No validations submitted")
# 		stateVD <- c("No validations submitted")

# 		#run the rmarkdown script to print html to file  
# 		rmarkdown::render(args[1])

# 	#Correct the name of the output files so that have the partner name and date
# 	rename_html <- paste0(date_now,"_",stateID,"_validation_specimensReport.html")
# 	file.copy(from="../../validation_reports.html", to= rename_html, overwrite = TRUE)

# }

