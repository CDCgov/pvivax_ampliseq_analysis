#!/usr/bin/python3

import os
import glob
import shutil
import re
import datetime
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from os import path
from collections import defaultdict


#read in files and create temporary directory
parser = ArgumentParser()
parser.add_argument("-g", "--newGenotypes")
parser.add_argument("-r", "--refGenotypes")
parser.add_argument("-t", "--tempPath")
parser.add_argument("-o", "--haplotypeSheet_output")
# parser.add_argument("-f", "--failedSpecimens")
parser.add_argument("-H", "--haplotypeDir")
# parser.add_argument("-Q", "--qcLogs")
args = parser.parse_args()

newPath = args.newGenotypes
refPath = args.refGenotypes
tempPath = args.tempPath
hapFoldPath = args.haplotypeDir



try:
	shutil.rmtree(tempPath)
except OSError:
	pass
os.mkdir(tempPath)

#get names of new genotypes and references, then copy to temp directory
genotypes = [os.path.basename(x) for x in glob.glob(os.path.join(args.newGenotypes, "*"))]
print(genotypes)
refs = [os.path.basename(x) for x in glob.glob(os.path.join(args.refGenotypes, "*"))]
print(refs)

for f in genotypes:
	shutil.copy(path.join(newPath,f),tempPath)

for f in refs:
	shutil.copy(path.join(refPath,f), tempPath)

print("after copy")
#combine reference and new genotype names into 1 variable, remove name if file is empty
Seq_ID = refs + genotypes
Seq_ID_full = [s for s in Seq_ID if os.stat(path.join(tempPath,s)).st_size >= 200000 ] #get list of files that are not empty

#next step is to read in blast files 
blast = ""
for fname in Seq_ID_full:
	with open(path.join(tempPath,fname)) as infile:
		blast += infile.read()

#remove empty line from initialization
blast = os.linesep.join([s for s in blast.splitlines() if s])

#get list of unique markers, convert to pandas dataframe
markers = []
for line in blast.splitlines():
	x = line.split("\t")[0]
	markers.append(x)
a = np.array(markers)
uniques = np.unique(a)
df = pd.DataFrame(index = Seq_ID_full, columns = uniques)

#loop through markers and add an X if specimen has that marker
#d dictionary is to support report generation
d = defaultdict(list)
for add_x in Seq_ID_full:
	with open(path.join(tempPath, add_x)) as infile:
		data = infile.read()
		for currentMarker in uniques:
			if re.search(currentMarker, data):
				df.at[add_x, currentMarker] = "X"
				d[add_x].append(currentMarker)



#list of all unique haplotypes
# allHaps = uniques.tolist()
#get unique base names with PART for each marker
#second split is because junction haplotypes are named slightly different than the rest, so takes extra step
# base1 = [i.split("_Hap")[0] for i in allHaps]
# base2 = [i.split("mt")[0] for i in base1]
# uniq1 = np.unique(base1).tolist()
# print(uniq1)

### The following step is specific to Cyclospora and the acceptance/rejection criteria using the 8 markers published in XXXX### 
### Please be thorough if editing this step to use more markers
### Result of script is to get list of specimens that will be filtered out in module 2. This is used to create reports. No specimens are actually removed in this script, all removing of specimens occurs in module 2
# filteredOut = []
# qcOut = []
# for key in d:
# 	#steps to get basename of markers for each specimen
# 	test = d.get(key)
# 	result1 = [i.split("_Hap")[0] for i in test] 
# 	result2  = [i.split("_PART")[0] for i in result1]
# 	result3 = [i.split("mt")[0] for i in result2]
# 	# Commented out - use filtering criteria for cyclospora. Specimens filtered out if it has less than 4 markers, or if it has exactly 4 markers and 2 or more are the Nu_CDS markers.
# 	#Current approach is to classify specimens as fail if they have less than 5 markers.
# 	if len(np.unique(result3)) < 5:
# 		result4 = [i.split("DS")[0] for i in result3]
# 		list1 = np.unique(result4).tolist()
# 		string1 = " ".join(list1)
# 		filteredOut.append(key)
# 		# if string1.count("Nu_C") > 1:
# 		# 	filteredOut.append(key)
# 		# elif len(np.unique(result4)) < 4:
# 		# 	filteredOut.append(key)
# 	if "NG" in key:
# 		# print("negatives")
# 		# print(key)
# 		if len(np.unique(test)) > 1:
# 			# print("fail negative")
# 			# print(key)
# 			qcOut.append(key)
# 	elif "PS" in key:
# 		# print("positives")
# 		if len(np.unique(result3)) < 5:
# 			# print("fail positive")
# 			# print(key)
# 			qcOut.append(key)

#Get genotype of the specimens that are filtered out. These are the markers that did amplify for those specimens 
# finalFilters = defaultdict(list)
# for m in filteredOut:
# 	tempFilter = []
# 	with open(path.join(tempPath, m)) as inFile:
# 		for row in inFile:
# 			tempFilter.append(row.split()[0])
# 	marker1 = [i.split("_Hap")[0] for i in tempFilter] 
# 	marker2 = [i.split("mt")[0] for i in marker1]
# 	marker3 = np.unique(marker2).tolist()
# 	finalFilters[m].append(marker3)

# #Get markers that did not amplify for the filtered out specimens. 
# filteredNotSeq = defaultdict(list)
# for key in finalFilters:
# 	a = finalFilters.get(key)
# 	noSeq = list(set(uniq1).difference(a[0]))
# 	filteredNotSeq[key].append(noSeq)

# #Convert each of those filtered out specimen dictionaries to a pandas dataframe
# if len(filteredOut) == 0:
# 	notSeq_df = pd.DataFrame()
# 	df2 = pd.DataFrame()
# else:
# 	notSeq_df = pd.DataFrame([[k] + v[0] for k,v in filteredNotSeq.items()])
# 	notSeq_df = notSeq_df.set_index(list(notSeq_df)[0])
# 	df2 = pd.DataFrame([[k] + v[0] for k,v in finalFilters.items()])
# 	df2 = df2.set_index(list(df2)[0])


#Add specimens with 0 markers to bottom of the haplotype sheet, as well as the two dataframes of filtered out specimens. For report generation
# Seq_ID_empty = [s for s in Seq_ID if os.stat(path.join(tempPath,s)).st_size == 0 ]
# for empty in Seq_ID_empty:
	# df = df.append(pd.Series(name = empty, dtype = "string"))
	# df = pd.concat([df, pd.Series(name = empty)])
	# df = df.append(pd.Series(name = empty))
	# df2 = df2.append(pd.Series(name = empty, dtype = "string"))
	# notSeq_df = notSeq_df.append(pd.Series(data = uniq1, name = empty))

#fix naming, then write to files
# df2 = df2.replace(to_replace = ["Mt_C"], value = ["Mt_Cmt"])
# notSeq_df = notSeq_df.replace(to_replace = ["Mt_C"], value = ["Mt_Cmt"])



#make sure cells without X in haplotype sheet are left empty. Set column/index name as Seq_ID
df = df.fillna('')
df.columns.name = 'Seq_ID'
df = df[~df.index.duplicated(keep="first")]

#get date to add to file name and write to file			
getDay = datetime.datetime.now()
today = str(getDay.strftime("%Y") + "-" + getDay.strftime("%m") + "-" + getDay.strftime("%d") + "-"+  getDay.strftime("%H") +  getDay.strftime("%M") + "_")
haplotypeFolder = hapFoldPath 
folderDay = os.path.join(haplotypeFolder, today)
fileName = folderDay + str(args.haplotypeSheet_output)

fileNameLocal = today +str(args.haplotypeSheet_output) 
df.to_csv(fileNameLocal, index = True, header = True, index_label = 'Seq_ID', sep =  "\t")

# failedPath = args.failedSpecimens
# failedSpecimen_markersPresent = today + "failedSpecimens_markersPresent.csv"
# failedSpecimen_markersAbsent = today + "failedSpecimens_markersAbsent.csv"
# presentName = os.path.join(failedPath, failedSpecimen_markersPresent)
# absentName = os.path.join(failedPath, failedSpecimen_markersAbsent)

# PS_NG_QC = today + "PS_NG_qcFail.txt"
# PS_NG_QC_name = os.path.join(args.qcLogs, PS_NG_QC)
# psNG_df = pd.DataFrame(qcOut)
# psNG_df.columns.name = "Controls FAIL QC"

# psNG_df.to_csv(PS_NG_QC_name, header = True)

# df2.to_csv(presentName, index =True, header = False)
# notSeq_df.to_csv(absentName, index = True, header = False)
#remove temp directory
shutil.rmtree(tempPath)
