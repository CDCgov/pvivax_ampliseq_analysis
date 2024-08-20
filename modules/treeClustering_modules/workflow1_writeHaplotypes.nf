#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Give path to folders used in this module
NEW_HAPS_DIR_Junc = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/JUNC_New_haps")
nonJunc_NEWHAPS = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/NON_JUNC_MAPPING/nonJunc_NEW_HAPS")
tmpGenotypes = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/TMP_GENOTYPES")
NEW_HAPS_DIR = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/REF_SEQS/BLASTING/NEW_HAPS")
finalGenotypes = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/SPECIMEN_GENOTYPES")

process combineGenotypes {
  publishDir "${tmpGenotypes}", pattern: "*GENOTYPE_final", mode: "copy"
  
  input:
  each inFASTQ
  path genotypes1
  path genotypes2


  output:
  path "*GENOTYPE_final", emit: finalGenotypes
  path "*GENOTYPE_final"

  script:
  """
  cat ${tmpGenotypes}/${inFASTQ.simpleName}.GENOTYPE_noJunction ${tmpGenotypes}/${inFASTQ.simpleName}.GENOTYPE_withJunction >> ${inFASTQ.simpleName}.GENOTYPE_final
  """
}

process copyGenotypes{
  input:
  path endGenotype

  output:
  path 'endMarker', emit: endFile

  script:
  """
  cp "${tmpGenotypes}/${endGenotype.simpleName}.GENOTYPE_noJunction" "${finalGenotypes}/${endGenotype.simpleName}"
  echo "done" >> endMarker
  """
}


// Call python script to make haplotype sheets. 
// Then copy newly identified haplotypes to the new haplotypes folder
process haplotypeSheet {
 label 'bioinformaticProcessing'

 publishDir "$baseDir/haplotype_sheets", pattern: "*haplotype_sheet*.txt", mode: "copy"

  input:
  path endFile

  output:
  path "*haplotype_sheet*.txt", emit: finalHapSheet
  path "*haplotype_sheet_largeFiles.txt", emit: finalHapSheetLarge

  script:
  """
  python3 ${params.pyscripts}/treeClustering_scripts/haplotypeSheet_gen_forNextflow_vivax.py -g $baseDir/HAPLOTYPE_CALLER_PVIVAX/SPECIMEN_GENOTYPES -r ${params.hapSheet_specifcReferences} -t $baseDir/HAPLOTYPE_CALLER_PVIVAX/tempDir -o haplotype_sheet.txt  -H $baseDir/haplotype_sheets
  python3 ${params.pyscripts}/treeClustering_scripts/haplotypeSheet_gen_forNextflow_vivaxOnlyLargeFiles.py -g $baseDir/HAPLOTYPE_CALLER_PVIVAX/SPECIMEN_GENOTYPES -r ${params.hapSheet_specifcReferences} -t $baseDir/HAPLOTYPE_CALLER_PVIVAX/tempDir -o haplotype_sheet_largeFiles.txt  -H $baseDir/haplotype_sheets 


  newNonJuncs=`ls ${nonJunc_NEWHAPS}/*.fasta | wc -l`
  newJuncs=`ls ${NEW_HAPS_DIR_Junc}/*.fasta  | wc -l`

  if [ \$newNonJuncs -gt 0 ]
  then
  ls ${nonJunc_NEWHAPS}/*.fasta | awk -F "/" '{print\$NF}' | while read line; do cp ${nonJunc_NEWHAPS}/\$line ${NEW_HAPS_DIR}/\$line
  done
  fi

  if [ \$newJuncs -gt 0 ]
  then
  ls ${NEW_HAPS_DIR_Junc}/*.fasta | awk -F "/" '{print\$NF}' | while read line; do cp ${NEW_HAPS_DIR_Junc}/\$line ${NEW_HAPS_DIR}/\$line
  done
  fi


  """
}

process haplotypeSheet_specific {
 label 'bioinformaticProcessing'

 publishDir "$baseDir/haplotype_sheets", pattern: "*haplotype_shee*.txt", mode: "copy"

  input:
  // path endFile

  output:
  path "*haplotype_sheet.txt", emit: finalHapSheet
  path "*haplotype_sheet_largeFiles.txt", emit: finalHapSheetLarge

  script:
  """
  python3 ${params.pyscripts}/treeClustering_scripts/haplotypeSheet_gen_forNextflow_vivax.py -g ${params.hapSheet_specifcGenotypes} -r ${params.hapSheet_specifcReferences} -t $baseDir/HAPLOTYPE_CALLER_PVIVAX/tempDir -o haplotype_sheet.txt  -H $baseDir/haplotype_sheets 
  python3 ${params.pyscripts}/treeClustering_scripts/haplotypeSheet_gen_forNextflow_vivaxOnlyLargeFiles.py -g ${params.hapSheet_specifcGenotypes} -r ${params.hapSheet_specifcReferences} -t $baseDir/HAPLOTYPE_CALLER_PVIVAX/tempDir -o haplotype_sheet_largeFiles.txt  -H $baseDir/haplotype_sheets 


  """
}

