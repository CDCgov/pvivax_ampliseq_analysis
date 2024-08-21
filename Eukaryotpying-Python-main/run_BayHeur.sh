#!/bin/bash -l




euktypDir=/Users/joelbarratt/Documents/P_VIVAX_PROJECT/GENOME_PIPELINE/CDC-Complete-Cyclospora-typing-workflow-ALPHA-TEST-master/Complete_Cyclospora_typing_workflow_MacOS_High_Sierra_ALPHA_TEST/Eukaryotpying-Python-main/DISTCOMP # modify as required


#Define input and output dir and file names
pycodeDir=$euktypDir
haplosheetDir="$(dirname "$pycodeDir")"/haplotype_sheets
haplosheetFileName=$(ls $haplosheetDir | sort | tail -1)
haplosheetPathFN=$haplosheetDir/$haplosheetFileName
echo "The newest haplotype file is found:"
echo $haplosheetPathFN
outputDir="$(dirname "$pycodeDir")"/ensemble_matrices
outputFNBase=${haplosheetFileName/.txt/""}
outputhh="${outputDir}/${outputFNBase}_hh_matrix.csv" ## Outdated heuristic normalisation method -- can be ignored.
outputbh="${outputDir}/${outputFNBase}_bh_matrix.csv" ## Ensemble matrix
outputB="${outputDir}/${outputFNBase}_B_matrix.csv" ## Bayesian matrix
outputH="${outputDir}/${outputFNBase}_H_matrix.csv" ## Heuristic matrix



#Run the command -- first modify "markerList.txt to tell the code which markers are essential" 
python3 $euktypDir/Pycode_distcomp/BayesianHeuristic.py -infile $haplosheetPathFN -outbh $outputbh -outhh $outputhh \
-outB $outputB -outH $outputH -sampmeetlocirequire samplesMeetCutoff.txt \
-expectlocifile $euktypDir/Pycode_distcomp/markerList.txt -expectlocinumber 5 # edit the "expectlociflag" to tell the code how many of the markers within your essential "markerList.txt" are required for distance computation

#Immediately the genotypes meeting the minimum data requirements will be printed (i.e., genotypes with any 5 markers (-expectlocinumner 5) within the file "markerList.txt").

echo "Genetic distance computation Complete. "
