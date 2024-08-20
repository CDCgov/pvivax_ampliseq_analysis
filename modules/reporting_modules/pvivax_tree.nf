#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process updateGeoPredictionDB {
  label 'reporting'

  publishDir "${params.treeOut}", pattern: "*txt", mode: "copy"

  input:
  path existing
  // path newPredictions

  output:
  path "*_geoPrediction_simple.txt", emit: geoPredictionUpdate
  path "*_pvivax_travel.txt", emit: travelUpdate


  script:
  """
  Rscript ${params.pyscripts}/reporting_scripts/update_predictionDB.R $existing
  """
}

process makeTree_full {
  label 'ggtreeEnv'

  publishDir "${params.treeOut}", pattern: "*pdf", mode: "copy"

  input:
  path matrix
  path parnasClusters
  path simplePrediction
  path travelInfo
  
  output:
  path "*tree.pdf", emit: treeOut

  script:
  """
  Rscript ${params.pyscripts}/reporting_scripts/TREE.R $matrix $parnasClusters $simplePrediction $travelInfo
  """

}

process makeTree_states {
  label 'ggtreeEnv'

  publishDir "${params.treeOut}", pattern: "*pdf", mode: "copy"

  input:
  path matrix
  path parnasClusters
  path simplePrediction
  path travelInfo
  
  output:
  path "*Pvivax_ampliseq_tree.pdf", emit: stateTree

  script:
  """
  cp  ${params.pyscripts}/reporting_scripts/Pvivax_markdown.Rmd .
  Rscript ${params.pyscripts}/reporting_scripts/TREE_perState.R $matrix $parnasClusters $simplePrediction $travelInfo ${params.reportState}
  """

}
