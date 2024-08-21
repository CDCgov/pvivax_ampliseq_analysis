require(stringr)
require(gtools)
require(parallel)
library(tidyverse)

#UPDATED August 5, 2021

todays_date <- Sys.Date()


threads <- readLines("THREADS")
cyclone_dir <- readLines("TOP_DIR") 
chromosomes <- readLines("CHROMOSOMES")

dir.create("STARTING_HAP_SHEET")


hap_sheet_loc <- paste0(cyclone_dir,"/haplotype_sheets/")
newest_hap <- tail(list.files(hap_sheet_loc),1)
get_new_hap <- paste0(hap_sheet_loc, newest_hap)

file.copy(get_new_hap, "./STARTING_HAP_SHEET/")

setwd("./STARTING_HAP_SHEET/")


hap_sheet = read.csv(newest_hap, skip=0, stringsAsFactors = FALSE, sep = "\t")


for(k in 1:chromosomes){
	
	if(k < 10){
		
	dir.create(paste0("CHR0",k))
	this_chrom <- paste0("CHR0",k)	
	extract_these <- c(1, which(str_detect(colnames(hap_sheet), this_chrom), arr.ind = TRUE)) # find index of hap sheet for this chromosome only
    hap_sheet_this_chrom <- hap_sheet[, extract_these] # generate hap sheet for this single chromosome
    write_loc <- paste0("./", this_chrom,"/",todays_date,"_",this_chrom,"_haplotype_data_sheet.txt")
    write.table(hap_sheet_this_chrom, write_loc, quote=FALSE, sep='\t', row.names=F)
	
	} else {
	
	dir.create(paste0("CHR",k))
	this_chrom <- paste0("CHR",k)	
	extract_these <- c(1, which(str_detect(colnames(hap_sheet), this_chrom), arr.ind = TRUE)) # find index of hap sheet for this chromosome only
    hap_sheet_this_chrom <- hap_sheet[, extract_these] # generate hap sheet for this single chromosome
    write_loc <- paste0("./", this_chrom,"/",todays_date,"_",this_chrom,"_haplotype_data_sheet.txt")
    write.table(hap_sheet_this_chrom, write_loc, quote=FALSE, sep='\t', row.names=F)
		
	}
	
}


setwd("..")

