#!/usr/bin/env nextflow
nextflow.enable.dsl=2



process marsVOI_allStates {
  label 'reporting'

  publishDir "${params.reportOut}", pattern: "*Pvivax_MaRS_genotyping_report.html", mode: "copy"

  input:
  path gatkCountry
  path reportableSNPs

  output:
  path "*Pvivax_MaRS_genotyping_report.html"

  script:
  """
  cp ${params.pyscripts}/reporting_scripts/extractVOI.R .
  cp ${params.pyscripts}/reporting_scripts/Pvivax_markdown_VOI.Rmd .
  Rscript extractVOI.R  Pvivax_markdown_VOI.Rmd  ${reportableSNPs} ${gatkCountry}
  """
}

process marsVOI_singleState {
  label 'reporting'

  publishDir "${params.reportOut}", pattern: "*Pvivax_MaRS_genotyping_report.html", mode: "copy"

  input:
  path gatkCountry
  path reportableSNPs

  output:
  path "*Pvivax_MaRS_genotyping_report.html"

  script:
  """
  cp ${params.pyscripts}/reporting_scripts/extractVOI.R .
  cp ${params.pyscripts}/reporting_scripts/Pvivax_markdown_VOI.Rmd .
  Rscript extractVOI.R  Pvivax_markdown_VOI.Rmd  ${reportableSNPs} ${gatkCountry} ${params.reportState}
  """
}



process runReportGen_allStates {
  label 'reporting'
  publishDir "${params.reportOut}", pattern: "*pvivax_clustering_report*", mode: "copy"

  input:
  path metaIn
  path clustersIn
  path geoIn
  path treeIn

  output:
  path "*_pvivax_clustering_report.csv"
  path "*_pvivax_clustering_report.html"

  script:
  """
  cp  ${params.pyscripts}/reporting_scripts/Pvivax_markdown.Rmd .
  Rscript ${params.pyscripts}/reporting_scripts/generateReport.R $metaIn $clustersIn $geoIn $treeIn

  """

}

process runReportGen_singleState {
  //take tree frmo makeTree_states

  label 'reporting'
  publishDir "${params.reportOut}", pattern: "*pvivax_clustering_report*", mode: "copy"

  input:
  path metaIn
  path clustersIn
  path geoIn
  path treeIn

  output:
  path "*_pvivax_clustering_report.csv"
  path "*_pvivax_clustering_report.html"

  script:
  """
  Rscript ${params.pyscripts}/reporting_scripts/generateReport_perState.R $metaIn $clustersIn $geoIn $treeIn ${params.reportState}

  """

}
