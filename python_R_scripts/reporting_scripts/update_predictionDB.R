library(dplyr)
library(readxl)

args <- commandArgs(trailingOnly = T)


newFile<- tail(list.files("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/P_vivax_fullAnalysis/GeoPrediction_output/predictedOut/", pattern ="gatkRegion"),1)
# print(new)
time_now <- format(Sys.time(), "%H%M")
date_now <- Sys.Date()
name_output <- paste(date_now,"_", time_now,"_pvivax_geoPrediction_simple.txt", sep ="")
print(name_output)


existingPredictions <- read.table(args[1], header = T, sep = "\t")

# newPredictions <- read.table(args[2], header = T, sep = "\t")
newPredictions <- read.table(paste("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/P_vivax_fullAnalysis/GeoPrediction_output/predictedOut/",newFile, sep = ""), header = T, sep = "\t")

print("read new predictions")
newPredictions_clean <- newPredictions[,c(1,2)]

newPredictions_clean$Region_Prediction_1 <- gsub(" ", "_", newPredictions_clean$Region_Prediction_1)
newPredictions_clean$Sample <- gsub("-", "_", newPredictions_clean$Sample)


colnames(newPredictions_clean)[1] <- "Seq_ID"

print(colnames(existingPredictions))
print(colnames(newPredictions_clean))
combinedFile <- rbind(existingPredictions, newPredictions_clean)

print("combined file")

combinedFile_final <- combinedFile[!duplicated(combinedFile$Seq_ID), ]

write.table(combinedFile_final, file = name_output, quote = F, sep  ="\t", row.names =F)
# print(combinedFile_final)

fullMeta <- read_xlsx("/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/P_vivax_fullAnalysis/pvivax_allBatches_metadata_update.xlsx")
travelInfo <- fullMeta[,c(1,3)]
travelInfo$Seq_ID_final <- gsub("-","_",travelInfo$Seq_ID_final)
print("before last")
name_outputTravel <- paste(date_now,"_", time_now,"_pvivax_travel.txt", sep ="")
write.table(travelInfo, file = name_outputTravel, quote = F, sep  ="\t", row.names =F)



