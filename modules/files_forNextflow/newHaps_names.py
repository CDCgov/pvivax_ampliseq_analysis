#!/usr/bin/python3

import pandas as pd
from argparse import ArgumentParser
import os

#this script will take in all potential new haplotypes and name them in a way that is consistent with the original naming convention.
#But does so outside of the find SNPs loop
#this means some new haps maybe be identified twice but the python script will handle those cases and keep only a single new haplotype per sequence

#

parser = ArgumentParser()
parser.add_argument("-n", "--newHaplotypes")
parser.add_argument("-e", "--existingHaplotypes")
parser.add_argument("-o", "--newNameDir")
args = parser.parse_args()


inFile = pd.read_csv(args.newHaplotypes, sep = " ", header = None)
existing = pd.read_csv(args.existingHaplotypes, sep = " ", header = None) 

inFile.rename({0:'specimen', 1:'locus', 2:'sequence'}, axis = 1, inplace = True)

#Keep only one record of each new sequence
uniqNew = inFile.drop_duplicates(subset='sequence', keep = 'first')


for i in uniqNew['locus'].unique():
	#count how many times the locus appears in the existing haplotypes, necessary for naming
	matches = existing[0].str.count(i)
	existingHaps = sum(matches)

	#get the matches between the new locus and existing haplotypes for naming
	myMatch = uniqNew['locus'].str.match(i)
	m = myMatch.index[myMatch == True].tolist()
	for j in m:
		#make a counter for  the new haplotypes
		existingHaps += 1

		# the below steps are all necessary to print the files correctly. Probably a better way to do this
		locus1 = uniqNew.loc[[j],['locus']].values
		locus2 = ''.join(map(str,locus1))
		locus3 = locus2.replace('[', '')
		locus4 = locus3.replace(']', '')
		locus5 = locus4.replace("'", "")

		sequence1 = uniqNew.loc[[j],['sequence']].values
		sequence2 = ''.join(map(str,sequence1))
		sequence3 = sequence2.replace('[', '')
		sequence4 = sequence3.replace(']', '')
		sequence5 = sequence4.replace("'", "")

		specimen1 = uniqNew.loc[[j],['specimen']].values
		specimen2 = ''.join(map(str,specimen1))
		specimen3 = specimen2.replace('[', '')
		specimen4 = specimen3.replace(']', '')
		specimen5 = specimen4.replace("'", "")

		#write fasta to directory/file following existing naming convention
		finalSeq = ">" + locus5 + str(existingHaps) +  "\n" + sequence5 + "\n"
		finalName = specimen5 + "." + locus5 + str(existingHaps) + ".fasta"
		with open(os.path.join(args.newNameDir,finalName), "w") as outFasta:
			outFasta.write(finalSeq)
