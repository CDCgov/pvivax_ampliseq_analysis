Plasmodium_Ampliseq_Nextflow#!/bin/bash

#UPDATED August 5, 2021

#####################################################################################################################################################################
SOFTWARE_DIRECTORY=nextflow_testing  # do not modify this line   ####################################################
#####################################################################################################################################################################   


#Tell me the folder where you pasted and extracted the Complete genotyping workflow zip file.
working_directory=/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/P_vivax_fullAnalysis


## How many threads would you like to use.
number_of_threads=10


### Epsilon. This must be any number between 0 and 1.
#epsilon=0.3072
numer_of_chromosomes=14










#########################################################   DO NOT MODIFY BELOW THIS POINT


cd $working_directory/$SOFTWARE_DIRECTORY/Eukaryotpying-Python-main

rm -rf TMP_matrix
mkdir TMP_matrix

cd TMP_matrix


echo $working_directory/$SOFTWARE_DIRECTORY > TOP_DIR
echo $number_of_threads > THREADS
echo $numer_of_chromosomes > CHROMOSOMES
#echo $epsilon > EPSILON


Rscript ../chrome_by_chrome.r # this script splits the haplotype sheet up, so you have a single haplotype sheet for each chromosome.

cd STARTING_HAP_SHEET


#THE LOOP BELOW WILL GENERATE A GENETIC DISTANCE MATRIX FOR EACH CHROMOSOME SEPARATELY.

ls -d */ |  while read folders

do

cd $folders

cp *.txt $working_directory/$SOFTWARE_DIRECTORY/Eukaryotpying-Python-main/haplotype_sheets/

cd $working_directory/$SOFTWARE_DIRECTORY/Eukaryotpying-Python-main

bash run_BayHeur.sh

rm $working_directory/$SOFTWARE_DIRECTORY/Eukaryotpying-Python-main/haplotype_sheets/*.txt

cd $working_directory/$SOFTWARE_DIRECTORY/Eukaryotpying-Python-main/TMP_matrix/STARTING_HAP_SHEET

done


cp -r  $working_directory/$SOFTWARE_DIRECTORY/Eukaryotpying-Python-main/ensemble_matrices/*sheet_H_matrix.csv \
$working_directory/$SOFTWARE_DIRECTORY/Eukaryotpying-Python-main/TMP_matrix/
rm -rf $working_directory/$SOFTWARE_DIRECTORY/Eukaryotpying-Python-main/ensemble_matrices/*matrix.csv

cd $working_directory/$SOFTWARE_DIRECTORY/Eukaryotpying-Python-main/TMP_matrix

mkdir MATRICES

mv *.csv MATRICES

cd MATRICES



## THIS R SCRIPT BELOW WILL TAKE THE DISTANCE MATRICES GENERATED FOR EACH CHROMOSOME, AND WILL SCALE EACH DISTANCE MATRIX BY ITS ENTROPY
## ONCE SCALING IS COMPLETE, IT WILL GENERATE A SINGLE NORMALIZED DISTANCE MATRIX THAT WE CLUSTER TO PRODUCE A TREE.
## THAT DISTANCE MATRIX WILL BE PRINTED TO THE "ensemble_matrices" folder in the main ALPHA_TEST directory.
Rscript $working_directory/$SOFTWARE_DIRECTORY/Eukaryotpying-Python-main/SCALE_matrix_by_entropy.R




echo "EUKARYOTYPING COMPLETE!"

cd $working_directory/$SOFTWARE_DIRECTORY/Eukaryotpying-Python-main

#rm -rf TMP_matrix

#rm THREADS
#rm EPSILON

cd $working_directory/$SOFTWARE_DIRECTORY/

Rscript TREE.r
echo "Tree generated!"


