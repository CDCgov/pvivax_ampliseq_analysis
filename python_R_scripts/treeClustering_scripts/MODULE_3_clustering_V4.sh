
#!/bin/bash

#UP November 19, 2020

###################################################################################################################################
# do not modify the next 2 lines.                                                           #######################################
                                                           #######################################
###################################################################################################################################

###USER MUST MODIFY THE FOLLOWING LINES:



# TELL ME THE DIRECTORY WHERE YOU PASTED THE CYCLONE FOLDER. LITERALLY WHERE YOU UNZIPPED CYCLONE AND INTEND TO RUN/INSTALL IT.
cyclone_location=/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/P_vivax_fullAnalysis

# CYCLONE=nextflow_testing            
CLUSTERING=$cyclone_location/CLUSTERING   



# What number of clusters do you want to start at?
cluster_min=5


# What number of clusters do you want to finish at?
cluster_max=50


# Tell me the number of threads you would like to use
number_of_threads=11



## bootstraps - number of random selections of unbiased matrices do you want to make?
bootstraps=1000

matrixUse=`echo $1 | awk -F "/" '{print$NF}'`

















#########################################################   DO NOT MODIFY BELOW THIS POINT


LOC=$CLUSTERING

matrix_folder=$cyclone_location/ensemble_matrices
clusters_folder=$cyclone_location/clusters_detected


###########################################################################################################################################################
######                                                  ###################################################################################################
######  WRITE DIRECTORIES FOR TMP FILES & VARIABLES     ###################################################################################################
######                                                  ###################################################################################################
###########################################################################################################################################################


# cd $LOC

# rm -rf TMP_REP

# mkdir TMP_REP
# cd TMP_REP

echo $cluster_min > CLUSTER_MIN
echo $cluster_max > CLUSTER_MAX
echo $matrix_folder > MATRIX_LOCATION
echo $clusters_folder > CLUSTER_LOCATION
echo $number_of_threads > THREADS
echo $bootstraps > BOOTSTRAPS
echo $LOC > SCRIPTS_LOC



#Rscript $LOC/CLUSTER_FINDER.R

# Rscript $LOC/CLUSTER_FINDER_V2.R $cyclone_location/$CYCLONE/newTest.csv
Rscript $LOC/CLUSTER_FINDER_V2.R $matrixUse

# cd $LOC

# rm -rf TMP_REP

# cd $cyclone_location/$CYCLONE

# echo "CLUSTERING COMPLETE!"




