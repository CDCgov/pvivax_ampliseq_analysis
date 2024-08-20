#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process individualMats {
  label 'module2'
  publishDir "${params.individualMatrixFolder}", pattern: "*H_matrix.csv", mode: "copy"

  input:
  path latestSheet

  output:
  path '*H_matrix.csv', emit: matrixOut
  // path '*.csv'

  script:
  """
  python3 $baseDir/Eukaryotpying-Python-main/DISTCOMP/Pycode_distcomp/BayesianHeuristic.py -infile ${latestSheet} -outbh ${latestSheet.simpleName}_bh_matrix.csv -outhh ${latestSheet.simpleName}_hh_matrix.csv -outB ${latestSheet.simpleName}_B_matrix.csv -outH ${latestSheet.simpleName}_H_matrix.csv -sampmeetlocirequire samplesMeetCutoff.txt -expectlocifile ${params.wf2MarkerList} -expectlocinumber ${params.wf2MinLociNum}
  """
}

process individualChroms {
  label 'module2'

  input:
  path latestSheet
  each myChrom

  output:
  path "*Nu*haplotype_data_sheet.txt", emit: chromHapSheet

  script:
  """

  python3 ${params.pyscripts}/treeClustering_scripts/extractChroms.py -i ${latestSheet} -c ${myChrom}


  """
}

process scaleMatrices {
  label 'rScaling'
  publishDir "${params.matrixFolder}", pattern: "*ALL_CHR_pairwisedistancematrix_H.csv", mode: "copy"

  input:
  path matrixList
  path haplotypeList

  output:
  path "*_ALL_CHR_pairwisedistancematrix_H.csv", emit: finalMatrix

  script:
  """
  mkdir MATRICES
  cp ${matrixList} MATRICES

  mkdir STARTING_HAP_SHEET
  cp ${haplotypeList} STARTING_HAP_SHEET

  Rscript $baseDir/Eukaryotpying-Python-main/SCALE_matrix_by_entropy.R

  """
}
