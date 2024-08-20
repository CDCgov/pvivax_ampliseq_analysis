#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process runModule3_Pvivax {
  label 'ggtreeEnv'

  input:
  path latestMatrix

  output:
  path "*tree.newick", emit: tree
  path "*Threshold.csv", emit: threshold

  script:
  """

  bash ${params.pyscripts}/treeClustering_scripts/MODULE_3_clustering_V4.sh ${latestMatrix}

  """ 
}

process runModule3_parnas{
  label 'parnasEnv'

  publishDir "${params.clusterFolder}", pattern: "*clusters.txt", mode: "copy"

  input:
  path tree
  path cutoff

  output:
  path "*clusters.txt", emit: clustersOutput

  script:
  """
  cutoffVar=`tail -1 ${cutoff}`

  echo \$cutoffVar > cutoffExample.txt

  parnas -t ${tree} --cover --radius \$cutoffVar  --clusters  ${cutoff.simpleName}_clusters.txt

  """ 
}

// process updateGeoPredictionDB {
//   label 'reporting'

//   publishDir "${params.treeOut}", pattern: "*txt", mode: "copy"

//   input:
//   path existing
//   // path newPredictions

//   output:
//   path "*_geoPrediction_simple.txt", emit: geoPredictionUpdate
//   path "*_pvivax_travel.txt", emit: travelUpdate


//   script:
//   """
//   Rscript ${params.pyscripts}/reporting_scripts/update_predictionDB.R $existing
//   """
// }

// process makeTree_full {
//   label 'ggtreeEnv'

//   publishDir "${params.treeOut}", pattern: "*pdf", mode: "copy"

//   input:
//   path matrix
//   path parnasClusters
//   path simplePrediction
//   path travelInfo
  
//   output:
//   path "*tree.pdf", emit: treeOut

//   script:
//   """
//   Rscript ${params.pyscripts}/reporting_scripts/TREE.R $matrix $parnasClusters $simplePrediction $travelInfo
//   """

// }

// process makeTree_states {
//   label 'ggtreeEnv'

//   publishDir "${params.treeOut}", pattern: "*pdf", mode: "copy"

//   input:
//   path matrix
//   path parnasClusters
//   path simplePrediction
//   path travelInfo
  
//   output:
//   path "*Pvivax_ampliseq_tree.pdf", emit: stateTree

//   script:
//   """
//   cp  ${params.pyscripts}/reporting_scripts/Pvivax_markdown.Rmd .
//   Rscript ${params.pyscripts}/reporting_scripts/TREE_perState.R $matrix $parnasClusters $simplePrediction $travelInfo ${params.reportState}
//   """

// }


// process runReportGen_allStates {
//   label 'reporting'
//   publishDir "${params.reportOut}", pattern: "*pvivax_clustering_report*", mode: "copy"

//   input:
//   path metaIn
//   path clustersIn
//   path geoIn
//   path treeIn

//   output:
//   path "*_pvivax_clustering_report.csv"
//   path "*_pvivax_clustering_report.html"

//   script:
//   """
//   cp  ${params.pyscripts}/reporting_scripts/Pvivax_markdown.Rmd .
//   Rscript ${params.pyscripts}/reporting_scripts/generateReport.R $metaIn $clustersIn $geoIn $treeIn

//   """

// }

// process runReportGen_singleState {
//   //take tree frmo makeTree_states

//   label 'reporting'
//   publishDir "${params.reportOut}", pattern: "*pvivax_clustering_report*", mode: "copy"

//   input:
//   path metaIn
//   path clustersIn
//   path geoIn
//   path treeIn

//   output:
//   path "*_pvivax_clustering_report.csv"
//   path "*_pvivax_clustering_report.html"

//   script:
//   """
//   Rscript ${params.pyscripts}/reporting_scripts/generateReport_perState.R $metaIn $clustersIn $geoIn $treeIn ${params.reportState}

//   """

// }
