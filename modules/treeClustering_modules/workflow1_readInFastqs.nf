#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Give path to folders used in this module
processReads = file("$projectDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/PROCESSED_READS")
tmpGenotypes = file("$projectDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/TMP_GENOTYPES")
sanitizedDir = file("$projectDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/SANITIZED_READS")
qcFolder = file("$projectDir/HAPLOTYPE_CALLER_PVIVAX/QClogs")


process bbduk {
  label 'bioinformaticProcessing'

  publishDir "${processReads}", pattern: "*clean_merged.fastq", mode: "copy"
  publishDir "${tmpGenotypes}", pattern: "*GENOTYPE*", mode: "copy"
  input:
  tuple val(sample_id), path(sample_files)

  output:
  tuple val(sample_id), path("*.clean_merged.fastq")
  path "*.clean_merged.fastq", emit: readyReads
  path "*GENOTYPE*"
  
  script:
  """
  bbduk.sh -h > bbdukWorking.file   
  which samtools > samtoolsWorking.file
  echo "testing" > justTesting.file

  echo ${sample_files[0]} > file1.txt
  echo ${sample_files[1]} > file2.txt

  bbduk.sh in1=${sample_files[0]} in2=${sample_files[1]} ref=${params.adapters} out1=${sample_id}.clean1.bbduk.fastq out2=${sample_id}.clean2.bbduk.fastq \
  minlen=50 ktrim=r k=23 mink=11 hdist=1 lhist=lhist.txt -Xmx20g
  
  newName="\$(echo ${sample_id} | cut -c 1-${params.nameLength} | sed 's/-/_/g')"

  cat ${sample_id}.clean1.bbduk.fastq ${sample_id}.clean2.bbduk.fastq > \$newName.clean_merged.fastq

  touch \$newName.GENOTYPE_noJunction
  touch \$newName.GENOTYPE_withJunction
  touch \$newName.GENOTYPE_final
  """
}


process fastqc_beforeMerge {
  label 'bioinformaticProcessing'

  publishDir "${qcFolder}", pattern: "*qualityCheck.txt", mode: "copy"

  input:
  tuple val(sample_id), path(sample_files)

  output:
  path "*qualityCheck.txt"

  shell:
  '''
  mkdir fastqc_out
  fastqc !{sample_files} -o fastqc_out
  for x in fastqc_out/*.zip; do unzip "$x"; done
  echo !{sample_files[0].simpleName}_fastqc > myName

  echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' !{sample_files[0].simpleName}_fastqc/fastqc_data.txt | grep -v "^>" | grep -v "^#B" |  awk '{print$2}' | awk '{s+=$1}END{print s/NR}'` > !{sample_id}.clean1_qualityCheck.txt   
  grep "Total Sequences" !{sample_files[0].simpleName}_fastqc/fastqc_data.txt >> !{sample_id}.clean1_qualityCheck.txt
  grep "Sequence length" !{sample_files[0].simpleName}_fastqc/fastqc_data.txt >> !{sample_id}.clean1_qualityCheck.txt
  
  echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' !{sample_files[1].simpleName}_fastqc/fastqc_data.txt | grep -v "^>" | grep -v "^#B" |  awk '{print$2}' | awk '{s+=$1}END{print s/NR}'` > !{sample_id}.clean2_qualityCheck.txt  
  grep "Total Sequences" !{sample_files[1].simpleName}_fastqc/fastqc_data.txt >> !{sample_id}.clean2_qualityCheck.txt
  grep "Sequence length" !{sample_files[1].simpleName}_fastqc/fastqc_data.txt >> !{sample_id}.clean2_qualityCheck.txt

  '''

}

process mapPhix {
  label 'bioinformaticProcessing'

 // publishDir "${qcFolder}", pattern: "*phix.sam", mode: "copy"
  publishDir "${qcFolder}", pattern: "*phix.out", mode: "copy"

  input:
  tuple val(sample_id), path(sample_files)

  output:
  path "*phix.out"
  //path "*phix.sam"

  script:
  """
  bowtie2 -p 5 -x $projectDir/modules/files_forNextflow/illumina_phiX.bt_index -1 ${sample_files[0]} -2 ${sample_files[1]} -S ${sample_id}.phix.sam &> ${sample_id}.phix.out
  """

}

process fastqc_afterMerge {
  label 'bioinformaticProcessing'

  publishDir "${qcFolder}", pattern: "*qualityCheck_postBBDUK.txt", mode: "copy"

  input:
  path sample_files

  output:
  path "*qualityCheck_postBBDUK.txt"
  path "*passQC.fastq", emit: qcFASTQs, optional: true

  shell:
  '''
  mkdir fastqc_out
  fastqc !{sample_files} -o fastqc_out
  for x in fastqc_out/*.zip; do unzip "$x"; done
  echo !{sample_files.baseName} > myName

  echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' !{sample_files.baseName}_fastqc/fastqc_data.txt | grep -v "^>" | grep -v "^#B" |  awk '{print$2}' | awk '{s+=$1}END{print s/NR}'` > !{sample_files.simpleName}.qualityCheck_postBBDUK.txt   
  grep "Total Sequences" !{sample_files.baseName}_fastqc/fastqc_data.txt >> !{sample_files.simpleName}.qualityCheck_postBBDUK.txt
  grep "Sequence length" !{sample_files.baseName}_fastqc/fastqc_data.txt >> !{sample_files.simpleName}.qualityCheck_postBBDUK.txt


  passReads=`grep "Total Sequences" !{sample_files.baseName}_fastqc/fastqc_data.txt | awk -F "\t" '{print$2}'`
  passQC=`grep "Average Quality" !{sample_files.simpleName}.qualityCheck_postBBDUK.txt | awk -F " " '{print$3}' | awk -F "." '{print$1}'`

  if [ $passReads -gt 50000 ] && [ $passQC -gt 29 ]
  then
  cp !{sample_files} !{sample_files.simpleName}.passQC.fastq
  fi
  
  '''

}