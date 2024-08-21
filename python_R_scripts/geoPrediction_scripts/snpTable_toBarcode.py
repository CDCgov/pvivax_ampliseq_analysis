#!/usr/bin/python3

import pandas as pd
import time
import argparse
import os
import json
from datetime import datetime
now = datetime.now()

#This script will take the output of variants to table and extract the markers of interest (geoCombined) and format a csv file for BALK classifier

parser = argparse.ArgumentParser()
parser.add_argument("--snpFile")
parser.add_argument("--balkKey")
parser.add_argument("--sampleList")
parser.add_argument("--outfile")
parser.add_argument("--refBarcode")
# parser.add_argument("--vcfType")

args = parser.parse_args()


#The balk key has the ampliseq amplicon name, the SNP position in PvP01 and the SNP locus in bp from the ampliseq amplicon (it is +100bp of the actual SNP in the amplicon, beause the reference sequence has +100bp of the amplicon)
newKey = pd.read_csv(args.balkKey, sep = "\t")
uniqueLoci = newKey.snpLocation_rawChromosome_balkFormat.unique()

#vcf table. Replace single quotes on specimen names (if necessary)
vcfTable = pd.read_csv(args.snpFile, sep = "\t")
vcfTable.columns = vcfTable.columns.str.replace("'", "")

#list of all samples in the vcf table
mySamps = pd.read_csv(args.sampleList, sep = "\t", header = None)

#Get the ampliseq amplicon name and position in a single column
newKey['exactName'] = newKey['joelRef_name'].astype(str) + "_" + newKey['snpLocus_distStartJoelRef'].astype(str) 
vcfTable['exactName'] = vcfTable['CHROM'].astype(str)  + "_" + vcfTable['POS'].astype(str) 

#Keep only the loci in both files - this will keep only the ampliseq loci in the balkKey file, because the vcf table has every single variant detect
result = newKey.merge(vcfTable, left_on = "exactName", right_on = "exactName", how = "inner")

#Dictionary of ampliseq amplicon and PvP01 locus
renameDict = dict(zip(newKey.exactName, newKey.snpLocation_rawChromosome_balkFormat))


#New df without the genotype calls. Pretty sure this isn't necessary
# getDepth = result.loc[:,~result.columns.str.endswith('GT')]
result = result.fillna(0)
sampGT_dict = {}
for i in mySamps[0]:
	print(i)
	#get the genotype and depth at each locus for each sample. 
	sampGT_list = []
	sampDF = result.loc[:,result.columns.str.startswith(i)]
	sampDF['joelRef_name']  = result['joelRef_name']
	sampDF['snpLocation_rawChromosome_balkFormat'] = result['snpLocation_rawChromosome_balkFormat']
	for j in uniqueLoci:
		# print(j)
		#loop through each locus and find the genotype calls for that sample at the locus
		try:
			lociDF = sampDF[sampDF.snpLocation_rawChromosome_balkFormat == j]
			# print(lociDF)
			#since there are often two amplicons that cover the same SNP, keep the genotype call that has the greater depth of coverage.
			lociDF2=lociDF.loc[lociDF[".".join([i,'DP'])].idxmax()]
			# print(lociDF2)
			# print("before myGT")
			myGT = lociDF.iloc[0,0]
			#Convert the GT call into nucleotide/X (for missing)/N (for heterozygous)
			# print(myGT)
			if myGT == "./.":
				finalGT = "X"
			elif myGT == ".|.":
				finalGT = "X"
			elif myGT[0] != myGT[2]:
				finalGT = "N"
			else:
				finalGT = myGT[0]
		except:
			print(j)
			print("skip locus")
			finalGT = "X"
		sampGT_list.append(finalGT)
	sampGT_dict[i] = sampGT_list

#format into CSV barcode for BALK classifier
toTranspose = pd.DataFrame.from_dict(sampGT_dict)
toTranspose['Locus'] = uniqueLoci
toTranspose.set_index('Locus', inplace =True)
mySampsFinal = toTranspose.T

#put the columns in the same order as the geoCombined classifier. Not sure if necessary but not a big lift. Just make sure the correct file is used
df = pd.read_csv(args.refBarcode)
lociList = list(df.columns)
newDF = mySampsFinal.reindex(columns=lociList)
newDF['Sample'] = newDF.index

date_time = now.strftime("%Y-%m-%d_%H%M")
newDF.to_csv('_'.join([date_time, args.outfile]), index = False)
vcfTable.to_csv('_'.join([date_time,args.snpFile]), index = False)





