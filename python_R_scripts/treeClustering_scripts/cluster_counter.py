#!/usr/bin/python3
#/Users/brycematlock/Documents/UPDx_Workflow/FINAL_2016_2244_1_S1_L001_..00000000.run_summary/2016_2244_1_S1_L001_..clusters.updated.clstr
#/Users/brycematlock/Documents/UPDx_Workflow/FINAL_2016_2244_1_S1_L001_..00000000.run_summary/count_summary_2016_2244_1_S1_L001_..csv

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-x', type=argparse.FileType('r'), help='The updated clstr file to perform the count from')
parser.add_argument('-y', type=argparse.FileType('w'), help='Name of the summary file')
args = parser.parse_args()

# with open("HAPLOTYPE_CALLER_CYCLO_V2/TMP/STX11691_20.clean_merged_CLUSTERS.fq.clstr") as file:
# 	for 

#read in the augmented clstr file (can probably change this when playing around)
inFile = pd.read_csv(args.x, sep = ">")
# print(inFile)
#get the lines that have Cluster in them
test = inFile.filter(like = "Cluster")
# print(test)
#dropna will remove lines without Cluster, 
test = test.dropna()
#reset the actul index of the df and also get the index for each line (the index is used to calculate size later on)
test1 = test.reset_index()


# #put the first cluster into a row (right now it is the column header)f
newRow = pd.DataFrame({'index':0, 'Cluster0':"Cluster0"}, index = [0])
test1 = pd.concat([newRow, test1]).reset_index(drop = True)
# print(test1)
# #add the final row with the size of the original file so that you can calculate size for last cluster
test1.loc[len(test1.index)] = [len(inFile), "Final"]
# print(test1)
test1 = test1[test1["Cluster0"].str.contains("Read") == False]
# print(test1)
myDict = {}
for i in range(0,len(test1.index)-1):
#   #get the index of the cluster and the cluster coming after it
    cluster = test1.iloc[i,1]
    end = test1.iloc[i+1,0]
    start = test1.iloc[i,0]
# 	#cluster 0 has correct size, all other clusters need to subtract 1
    if cluster == "Cluster0":
        size = end - start
    else:
        size = end - start -1
    myDict[cluster] = size

# #conver to datafrme
final = pd.DataFrame.from_dict(myDict, orient = "index")
final.to_csv(args.y, header=False)