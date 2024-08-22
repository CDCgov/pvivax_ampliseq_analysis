#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Give path to folders used in this module
renamedReads = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/RENAMED_READS")
refHapsDir = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/REF_HAPS")
nonJunc_mapping = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/NON_JUNC_MAPPING")
nonJunc_NEWHAPS = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/NON_JUNC_MAPPING/nonJunc_NEW_HAPS")
nonJunc_clippedReads = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/NON_JUNC_MAPPING/CLIPPED_READS")
NEW_HAPS_DIR = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/REF_SEQS/BLASTING/NEW_HAPS")
tmpGenotypes = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/TMP_GENOTYPES")

process makeRefs {
  input:
  each(newHaps)
  
  output:
  path 'refHaps.fasta', emit: completeFile

  script:
  """
  cat ${newHaps} > all.fasta
  sed '/Junction/,+1 d' all.fasta > refHaps.fasta
  """
}

process makeRefs_v2 {
  input:
  // each(newHaps)
  
  output:
  path 'refHaps.fasta', emit: completeFile

  script:
  """
  cat ${NEW_HAPS_DIR}/*.fasta > all.fasta
  sed '/Junction/,+1 d' all.fasta > refHaps.fasta
  """
}

process finalNewRefs {
  input:
  path(allNew)

  output:
  path 'fullRefs.fasta', emit: outFinalRefs
  path "shortHapRefs.fasta", emit: outShortRefs

  script:
  """
  cat ${allNew} ${params.fullNonJunction_refs} > fullRefs.fasta
  cat ${allNew} ${params.originalHaplotypes_refs} > shortHapRefs.fasta
  """
}

process finalNewRefs_v2 {
  input:
  path(allNew)

  output:
  path 'fullRefs.fasta', emit: outFinalRefs
  path "shortHapRefs.fasta", emit: outShortRefs

  script:
  """
  cat ${allNew} ${params.fullNonJunction_refs} > fullRefs.fasta
  cat ${allNew} ${params.originalHaplotypes_refs} > shortHapRefs.fasta
  """
}

process mapping {
   label 'mappingSoftware'

  input:
  each(refHaps)
  // each(sample_files)
  tuple val(sampleID), path(sample_files)

  output:
  path '*.mapped_only.sam', emit: mappedOnlySam

  script:
  """

  which samtools > samtoolsTestingMapping.file
  which bwa > bwaTstingMapping.file

  bwa index ${refHaps}
  
  bwa mem -t ${params.threads} ${refHaps} ${sample_files[0]} ${sample_files[1]} > ${sampleID}.alignment.sam
  samtools view -h -F 4 ${sampleID}.alignment.sam | samtools view -bS > ${sampleID}.mapped_only.sam 

  samtools fastq ${sampleID}.mapped_only.sam > ${sampleID}.mapped_only.fastq
  """
}

// bwa mem -t ${params.threads} ${refHaps} ${sample_files} > ${sample_files.simpleName}.alignment.sam

// This sam to fastq step is rquired with renaming the sequences. If you just do use the samtools fastq convert, you will get incorrect haplotypes for some reason
process mapping2fastq {
   label 'mappingSoftware'
 
  publishDir "${nonJunc_mapping}", pattern: "*mapped_only.fastq", mode: "copy"

  input:
  path(samMapped)

  output:
  path "*mapped_only.fastq", emit: mappedFastq

  shell:
  '''
  samtools view !{samMapped} | awk '{print("@"\$1"\\n"\$10"\\n+\\n"\$11)}' | awk '{if(NR%4==1) \$0=sprintf("@Read_%d",(1+i++)); print;}' > !{samMapped.simpleName}.mapped_only.fastq
  '''
}

process bowtie2 {
  label 'mappingSoftware'

  publishDir "${nonJunc_mapping}", pattern: "*.sorted.alignment.*", mode: "copy"

  input:
  path(mappedFastq)

  output:
  path("*.alignment.bam")
  path("*.alignment.bam.bai")
  path "*.alignment.bam", emit: finishSorting

  script:
  """
  bowtie2-build "${params.fullNonJunction_refs}" originalHaps_bt_index
  bowtie2 -x originalHaps_bt_index -U ${mappedFastq.simpleName}.mapped_only.fastq --threads ${params.threads} --local > ${mappedFastq.simpleName}.alignment.bam
  samtools sort ${mappedFastq.simpleName}.alignment.bam > ${mappedFastq.simpleName}.sorted.alignment.bam
  samtools index ${mappedFastq.simpleName}.sorted.alignment.bam
  """
}
//code that works...
// java -jar ${params.biostar84452}  ${sorted.simpleName[0]}.\$myName.softClipped.bam > ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only.bam
// this step takes a while because it repeats this process for each part of each marker for each sample. But it is necessary to do it this way
process bedExtraction {
  label 'mappingSoftware'

  errorStrategy 'ignore'



  publishDir "${nonJunc_mapping}", pattern: "*coverag*", mode: "copy"

  publishDir "${nonJunc_mapping}", pattern: "*.f*", mode: "copy"



  input:

  val(marker) 

  each(sorted)



  output:

  path("*coverag*")

  path("*.f*")

  path "*multi-seq",  emit: multiOut, optional: true

  path("*mapped_only.fastq"), emit: finalClipped_fastq



  script:

  """

  myName="\$(echo $marker)"



  shortName="\$(echo \$myName | sed 's/_PART_*.//g')"

  grep \$myName ${params.bedFile} > \$myName.myMatch.bed

  left_trim="\$(cat \$myName.myMatch.bed | cut -f2)"

  right_trim="\$(cat \$myName.myMatch.bed | cut -f3)"



  echo \$left_trim \$right_trim > \$myName.myMatch.trim

  echo \$shortName > \$shortName.markerNames



  



  samtools view -bu -F 4 "${nonJunc_mapping}/${sorted.simpleName[0]}.sorted.alignment.bam" \$shortName | java -jar ${params.samjs} -e "record.alignmentStart <= \$left_trim && record.alignmentEnd >= \$right_trim" --samoutputformat BAM -o ${sorted.simpleName[0]}.\$myName.bam



  java -jar ${params.pcrclipreads} --bed \$myName.myMatch.bed ${sorted.simpleName[0]}.\$myName.bam --samoutputformat BAM -o ${sorted.simpleName[0]}.\$myName.softClipped.bam


  java -jar ${params.biostar84452}  ${sorted.simpleName[0]}.\$myName.softClipped.bam > ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only.bam



  samtools view -h -F 4 ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only.bam| samtools view -bS > ${sorted.simpleName[0]}.\$myName.mapped_only.sam



  samtools fastq ${sorted.simpleName[0]}.\$myName.mapped_only.sam > ${sorted.simpleName[0]}.\$myName.mapped_only.fastq





  samtools sort -T -n ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only.bam -o ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only_aln_sorted.bam



  samtools index ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only_aln_sorted.bam



  samtools depth -aa -d 0 ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only_aln_sorted.bam | awk "/^\$shortName/" | head -\$right_trim | awk "NR > \$left_trim" > ${sorted.simpleName[0]}.\$myName.coverage



  avgCov="\$(awk '{ total += \$3; count++ } END { print total/count }' ${sorted.simpleName[0]}.\$myName.coverage)"



  echo "scale=0; (\$avgCov * 0.10)/1" | bc > ${sorted.simpleName[0]}.\$myName.assignHaps_coverageCutoff

  echo "scale=0; (\$avgCov * 0.25)/1" | bc > ${sorted.simpleName[0]}.\$myName.newHaps_coverageCutoff

  newHaps_cutoff=`echo "scale=0; (\$avgCov * 0.25)/1" | bc`



  newHaps2_cutoff='\$(cat ${sorted.simpleName[0]}.\$myName.newHaps_coverageCutoff)'



  echo \$newHaps_cutoff > ${sorted.simpleName[0]}.\$myName.newHapsTest_coverageCutoff



  cd-hit-est -i ${sorted.simpleName[0]}.\$myName.mapped_only.fastq -o ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fq -c 1 -g 1 -d 0 -T ${params.threads} -s 1 -M ${params.RAM}





  clus="\$(wc -l ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fq | awk '{print\$1}')"



  if [ \$clus -gt 0 ]

  then

  seqret -sequence ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fq -outseq ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fasta



  seqret -sequence ${sorted.simpleName[0]}.\$myName.mapped_only.fastq -outseq ${sorted.simpleName[0]}.\$myName.mapped_only.fasta



  perl ${params.multiSeq_pl} ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fasta ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fq.clstr ${sorted.simpleName[0]}.\$myName.newHaps_multi-seq \$newHaps_cutoff

  

  fi

  """

}




process bedExtraction_v2 {
  label 'mappingSoftware'

  errorStrategy 'ignore'

  publishDir "${nonJunc_mapping}", pattern: "*coverag*", mode: "copy"
  publishDir "${nonJunc_mapping}", pattern: "*.f*", mode: "copy"

  input:
  val(marker) 
  each(sorted)

  output:
  path("*coverag*")
  path("*.f*")
  path "*multi-seq",  emit: multiOut, optional: true
  path("*mapped_only.fastq"), emit: finalClipped_fastq
  tuple val("${sorted.simpleName[0]}.${marker}"), path("*clusterCounted.csv"), path("*multi-seq"), emit:bedTuple, optional: true


  script:
  """
  myName="\$(echo $marker)"

  shortName="\$(echo \$myName | sed 's/_PART_*.//g')"
  grep \$myName ${params.bedFile} > \$myName.myMatch.bed
  left_trim="\$(cat \$myName.myMatch.bed | cut -f2)"
  right_trim="\$(cat \$myName.myMatch.bed | cut -f3)"

  echo \$left_trim \$right_trim > \$myName.myMatch.trim
  echo \$shortName > \$shortName.markerNames

  

  samtools view -bu -F 4 "${nonJunc_mapping}/${sorted.simpleName[0]}.sorted.alignment.bam" \$shortName | java -jar ${params.samjs} -e "record.alignmentStart <= \$left_trim && record.alignmentEnd >= \$right_trim" --samoutputformat BAM -o ${sorted.simpleName[0]}.\$myName.bam

  java -jar ${params.pcrclipreads} --bed \$myName.myMatch.bed ${sorted.simpleName[0]}.\$myName.bam --samoutputformat BAM -o ${sorted.simpleName[0]}.\$myName.softClipped.bam

  java -jar ${params.biostar84452}  ${sorted.simpleName[0]}.\$myName.softClipped.bam > ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only.bam

  samtools view -h -F 4 ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only.bam| samtools view -bS > ${sorted.simpleName[0]}.\$myName.mapped_only.sam

  samtools fastq ${sorted.simpleName[0]}.\$myName.mapped_only.sam > ${sorted.simpleName[0]}.\$myName.mapped_only.fastq


  samtools sort -T -n ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only.bam -o ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only_aln_sorted.bam

  samtools index ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only_aln_sorted.bam

  samtools depth -aa -d 0 ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only_aln_sorted.bam | awk "/^\$shortName/" | head -\$right_trim | awk "NR > \$left_trim" > ${sorted.simpleName[0]}.\$myName.coverage

  avgCov="\$(awk '{ total += \$3; count++ } END { print total/count }' ${sorted.simpleName[0]}.\$myName.coverage)"

  echo "scale=0; (\$avgCov * 0.10)/1" | bc > ${sorted.simpleName[0]}.\$myName.assignHaps_coverageCutoff
  echo "scale=0; (\$avgCov * 0.25)/1" | bc > ${sorted.simpleName[0]}.\$myName.newHaps_coverageCutoff
  newHaps_cutoff=`echo "scale=0; (\$avgCov * 0.25)/1" | bc`

  newHaps2_cutoff='\$(cat ${sorted.simpleName[0]}.\$myName.newHaps_coverageCutoff)'

  echo \$newHaps_cutoff > ${sorted.simpleName[0]}.\$myName.newHapsTest_coverageCutoff

  cd-hit-est -i ${sorted.simpleName[0]}.\$myName.mapped_only.fastq -o ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fq -c 1 -g 1 -d 0 -T ${params.threads} -s 1 -M ${params.RAM}


  clus="\$(wc -l ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fq | awk '{print\$1}')"

  if [ \$clus -gt 0 ]
  then
  seqret -sequence ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fq -outseq ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fasta

  seqret -sequence ${sorted.simpleName[0]}.\$myName.mapped_only.fastq -outseq ${sorted.simpleName[0]}.\$myName.mapped_only.fasta

  perl ${params.multiSeq_pl} ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fasta ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fq.clstr ${sorted.simpleName[0]}.\$myName.newHaps_multi-seq \$newHaps_cutoff

  cat ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fq.clstr  | sed 's/Cluster /Cluster/g' > ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters_toCount.fq.clstr
  python3 ${params.pyscripts}/treeClustering_scripts/cluster_counter.py -x ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters_toCount.fq.clstr -y ${sorted.simpleName[0]}.\$myName.clusterCounted.csv
  
  fi
  """
}


process bedExtraction_plusNewHaps {
  label 'mappingSoftware'
 

  errorStrategy 'ignore'

  publishDir "${nonJunc_mapping}", pattern: "*coverag*", mode: "copy"
  publishDir "${nonJunc_mapping}", pattern: "*.f*", mode: "copy"

  input:
  
  val(marker) 
  each(sorted)
  each(refDB)

  output:
  path("*coverag*")
  path("*.f*")
  path "*multi-seq",  emit: multiOut, optional: true
  path("*mapped_only.fastq"), emit: finalClipped_fastq
  tuple val("${sorted.simpleName[0]}.${marker}"), path("*clusterCounted.csv"), path("*multi-seq"), emit:bedTuple, optional: true
  tuple val("${sorted.simpleName[0]}"), path("*mapped_only.fastq"), emit: tryBedFastq, optional:true
  path "allNewHapsOut.file", emit: allNewHaps, optional: true


  script:
  """
  myName="\$(echo $marker)"

  shortName="\$(echo \$myName | sed 's/_PART_*.//g')"
  grep \$myName ${params.bedFile} > \$myName.myMatch.bed
  left_trim="\$(cat \$myName.myMatch.bed | cut -f2)"
  right_trim="\$(cat \$myName.myMatch.bed | cut -f3)"

  echo \$left_trim \$right_trim > \$myName.myMatch.trim
  echo \$shortName > \$shortName.markerNames

  

  samtools view -bu -F 4 "${nonJunc_mapping}/${sorted.simpleName[0]}.sorted.alignment.bam" \$shortName | java -jar ${params.samjs} -e "record.alignmentStart <= \$left_trim && record.alignmentEnd >= \$right_trim" --samoutputformat BAM -o ${sorted.simpleName[0]}.\$myName.bam

  java -jar ${params.pcrclipreads} --bed \$myName.myMatch.bed ${sorted.simpleName[0]}.\$myName.bam --samoutputformat BAM -o ${sorted.simpleName[0]}.\$myName.softClipped.bam

  java -jar ${params.biostar84452}  ${sorted.simpleName[0]}.\$myName.softClipped.bam > ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only.bam

  samtools view -h -F 4 ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only.bam| samtools view -bS > ${sorted.simpleName[0]}.\$myName.mapped_only.sam

  samtools fastq ${sorted.simpleName[0]}.\$myName.mapped_only.sam > ${sorted.simpleName[0]}.\$myName.mapped_only.fastq


  samtools sort -T -n ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only.bam -o ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only_aln_sorted.bam

  samtools index ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only_aln_sorted.bam

  samtools depth -aa -d 0 ${sorted.simpleName[0]}.\$myName.finalClipped.map_to_marker_only_aln_sorted.bam | awk "/^\$shortName/" | head -\$right_trim | awk "NR > \$left_trim" > ${sorted.simpleName[0]}.\$myName.coverage

  avgCov="\$(awk '{ total += \$3; count++ } END { print total/count }' ${sorted.simpleName[0]}.\$myName.coverage)"

  echo "scale=0; (\$avgCov * 0.10)/1" | bc > ${sorted.simpleName[0]}.\$myName.assignHaps_coverageCutoff
  echo "scale=0; (\$avgCov * 0.25)/1" | bc > ${sorted.simpleName[0]}.\$myName.newHaps_coverageCutoff
  newHaps_cutoff=`echo "scale=0; (\$avgCov * 0.25)/1" | bc`

  newHaps2_cutoff='\$(cat ${sorted.simpleName[0]}.\$myName.newHaps_coverageCutoff)'

  echo \$newHaps_cutoff > ${sorted.simpleName[0]}.\$myName.newHapsTest_coverageCutoff

  cd-hit-est -i ${sorted.simpleName[0]}.\$myName.mapped_only.fastq -o ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fq -c 1 -g 1 -d 0 -T ${params.threads} -s 1 -M ${params.RAM}


  clus="\$(wc -l ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fq | awk '{print\$1}')"

  if [ \$clus -gt 0 ]
  then
  seqret -sequence ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fq -outseq ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fasta

  seqret -sequence ${sorted.simpleName[0]}.\$myName.mapped_only.fastq -outseq ${sorted.simpleName[0]}.\$myName.mapped_only.fasta

  perl ${params.multiSeq_pl} ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fasta ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fq.clstr ${sorted.simpleName[0]}.\$myName.newHaps_multi-seq \$newHaps_cutoff

  cat ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters.fq.clstr  | sed 's/Cluster /Cluster/g' > ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters_toCount.fq.clstr
  python3 ${params.pyscripts}/treeClustering_scripts/cluster_counter.py -x ${sorted.simpleName[0]}.\$myName.clean_merged_Clusters_toCount.fq.clstr -y ${sorted.simpleName[0]}.\$myName.clusterCounted.csv
  
  fi


  ls ${sorted.simpleName[0]}.\$myName.newHaps_multi-seq  | grep -v ".fasta" > ${sorted.simpleName[0]}.\$myName.newHaps_multi-seq.list


  cat ${sorted.simpleName[0]}.\$myName.newHaps_multi-seq.list | while read i
  do

  clusDepth="\$(grep "Cluster\$i"  ${sorted.simpleName[0]}.\$myName.clusterCounted.csv | awk -F ',' '{print\$2}')"  

  echo \$clusDepth > myDepth_\$i


  cat ${sorted.simpleName[0]}.\$myName.newHaps_multi-seq/\$i > ${sorted.simpleName[0]}.\$myName.newHaps_multi-seq/\$i.fasta


  makeblastdb -in ${sorted.simpleName[0]}.\$myName.newHaps_multi-seq/\$i.fasta -input_type fasta -dbtype nucl

  blastn -db ${sorted.simpleName[0]}.\$myName.newHaps_multi-seq/\$i.fasta -query ${refDB} -evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -num_threads ${params.threads} -out ${sorted.simpleName[0]}.\$myName.\$i.singleMatch_blastOut -max_target_seqs 1 -word_size 7 -dust no -soft_masking false -outfmt "6 qseqid"

  match_present="\$(wc -l ${sorted.simpleName[0]}.\$myName.\$i.singleMatch_blastOut | awk '{print\$1}')"


  covCutoff=`cat "${sorted.simpleName[0]}.\$myName.newHaps_coverageCutoff" | awk '{print\$1}'`

  if [ \$match_present == 0  ] && [ \$clusDepth -gt ${params.depthNewHap_min} ] && [ \$clusDepth -gt \$covCutoff ]
  
  then 

  blastn -db ${sorted.simpleName[0]}.\$myName.newHaps_multi-seq/\$i.fasta -query ${refDB} -evalue 0.001 -perc_identity 75 -qcov_hsp_perc 85 -num_threads ${params.threads} -out ${sorted.simpleName[0]}.\$myName.\$i.newHaps_blast_RED_SIM -max_target_seqs 1 -word_size 7 -dust no -soft_masking false -outfmt "6 qseqid"
  
  theMarker=`head -1 ${sorted.simpleName[0]}.\$myName.\$i.newHaps_blast_RED_SIM | sed 's/_Hap_.*/_Hap_/'`
  theSequence=`cat ${sorted.simpleName[0]}.\$myName.newHaps_multi-seq/\$i | tail -3 | tr -d '\n' | sed 's/^>Read_[0-9]*//g'`
  theSequence_FDA=`cat ${sorted.simpleName[0]}.\$myName.newHaps_multi-seq/\$i | tail -1`
  theSpecimen=`echo ${sorted.simpleName[0]}`

  echo \$theSpecimen \$theMarker \$theSequence > allNewHapsOut.file

  fi
  done 
  """
}




process makeInitialReads {
  input:
  path(reads)

  output:
  path "*Initial_clipped_for_hap_calling.fastq", emit: readsOut

  script:
  """
  cat ${reads} > ${reads.simpleName}.Initial_clipped_for_hap_calling.fastq
  """
}

process makeFinalReads {
  publishDir "${nonJunc_clippedReads}", pattern: "*Final_clipped_for_hap_calling.fastq", mode: "copy"

  input:
  path(reads2)

  output:
  path("*Final_clipped_for_hap_calling.fastq")
  path("*Final_clipped_for_hap_calling.fastq"), emit: finalClipped, optional: true

  script:
  """
  cat ${reads2} > ${reads2.simpleName}.Final_clipped_for_hap_calling.fastq
  """
}

process makeFinalReads_v2 {
  publishDir "${nonJunc_clippedReads}", pattern: "*Final_clipped_for_hap_calling.fastq", mode: "copy"

  input:
  // path(reads2)
  tuple val(sample), path(reads)

  output:
  path("*Final_clipped_for_hap_calling.fastq")
  path("*Final_clipped_for_hap_calling.fastq"), emit: finalClipped, optional: true

  script:
  """
  echo ${reads} > myFile
  cat ${reads} > ${sample}.Final_clipped_for_hap_calling.fastq
  """
}

process newHaps1 {
  label 'mappingSoftware'

  errorStrategy 'ignore'

  input:
  each(newMulti)
  path(refDB)

  output:
  path "allNewHapsOut.file", emit: allNewHaps, optional: true

  script:
  """
  ls ${newMulti} | grep -v ".fasta" > ${newMulti}.list


  cat ${newMulti}.list | while read i
  do
  cat ${newMulti}/\$i > ${newMulti}/\$i.fasta

  makeblastdb -in ${newMulti}/\$i.fasta -input_type fasta -dbtype nucl

  blastn -db ${newMulti}/\$i.fasta -query ${refDB} -evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -num_threads ${params.threads} -out ${newMulti.baseName}.\$i.singleMatch_blastOut -max_target_seqs 1 -word_size 7 -dust no -soft_masking false -outfmt "6 qseqid"

  blastn -db ${newMulti}/\$i.fasta -query "${nonJunc_mapping}/${newMulti.baseName}.mapped_only.fasta" -evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -num_threads ${params.threads} -out ${newMulti.baseName}.\$i.newHaps_blast_coverage_guess -max_target_seqs 1 -word_size 7 -dust no -soft_masking false -outfmt "6 qseqid"

  match_present="\$(wc -l ${newMulti.baseName}.\$i.singleMatch_blastOut | awk '{print\$1}')"

  coverage_estimate_by_blast="\$(wc -l ${newMulti.baseName}.\$i.newHaps_blast_coverage_guess | awk '{print\$1}')"

  covCutoff=`cat "${nonJunc_mapping}/${newMulti.baseName}.newHaps_coverageCutoff" | awk '{print\$1}'`

  if [ \$match_present == 0  ] && [ \$coverage_estimate_by_blast -gt ${params.depthNewHap_min} ] && [ \$coverage_estimate_by_blast -gt \$covCutoff ]
  
  then 

  blastn -db ${newMulti}/\$i.fasta -query ${refDB} -evalue 0.001 -perc_identity 75 -qcov_hsp_perc 100 -num_threads ${params.threads} -out ${newMulti.baseName}.\$i.newHaps_blast_RED_SIM -max_target_seqs 1 -word_size 7 -dust no -soft_masking false -outfmt "6 qseqid"
  
  theMarker=`head -1 ${newMulti.baseName}.\$i.newHaps_blast_RED_SIM | sed 's/_Hap_.*/_Hap_/'`
  theSequence=`cat ${newMulti}/\$i | tail -3 | tr -d '\n' | sed 's/^>Read_[0-9]*//g'`
  theSequence_FDA=`cat ${newMulti}/\$i | tail -1`
  theSpecimen=`echo ${newMulti.simpleName}`

  echo \$theSpecimen \$theMarker \$theSequence > allNewHapsOut.file

  fi
  done 
  """
}

process newHaps1_v2{
  label 'mappingSoftware'

  errorStrategy 'ignore'

  input:
  // each(newMultiOrig)
  each(refDB)
  tuple val(sampleMarkerPart), path(counted), path(newMulti)


  output:
  path "allNewHapsOut.file", emit: allNewHaps, optional: true

  script:
  """
  ls ${newMulti} | grep -v ".fasta" > ${newMulti}.list


  cat ${newMulti}.list | while read i
  do

  clusDepth="\$(grep "Cluster\$i" ${counted} | awk -F ',' '{print\$2}')"  

  echo \$clusDepth > myDepth_\$i


  cat ${newMulti}/\$i > ${newMulti}/\$i.fasta


  makeblastdb -in ${newMulti}/\$i.fasta -input_type fasta -dbtype nucl

  blastn -db ${newMulti}/\$i.fasta -query ${refDB} -evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -num_threads ${params.threads} -out ${newMulti.baseName}.\$i.singleMatch_blastOut -max_target_seqs 1 -word_size 7 -dust no -soft_masking false -outfmt "6 qseqid"

  match_present="\$(wc -l ${newMulti.baseName}.\$i.singleMatch_blastOut | awk '{print\$1}')"


  covCutoff=`cat "${nonJunc_mapping}/${newMulti.baseName}.newHaps_coverageCutoff" | awk '{print\$1}'`

  if [ \$match_present == 0  ] && [ \$clusDepth -gt ${params.depthNewHap_min} ] && [ \$clusDepth -gt \$covCutoff ]
  
  then 

  blastn -db ${newMulti}/\$i.fasta -query ${refDB} -evalue 0.001 -perc_identity 75 -qcov_hsp_perc 85 -num_threads ${params.threads} -out ${newMulti.baseName}.\$i.newHaps_blast_RED_SIM -max_target_seqs 1 -word_size 7 -dust no -soft_masking false -outfmt "6 qseqid"
  
  theMarker=`head -1 ${newMulti.baseName}.\$i.newHaps_blast_RED_SIM | sed 's/_Hap_.*/_Hap_/'`
  theSequence=`cat ${newMulti}/\$i | tail -3 | tr -d '\n' | sed 's/^>Read_[0-9]*//g'`
  theSequence_FDA=`cat ${newMulti}/\$i | tail -1`
  theSpecimen=`echo ${newMulti.simpleName}`

  echo \$theSpecimen \$theMarker \$theSequence > allNewHapsOut.file

  fi
  done 
  """

}

process nameNewHaps {
  label 'bioinformaticProcessing'

  input:
  path(allNewHaps)
  path(existingHaps)

  output:
  path "newHapsDone.txt", emit: nonJunc_hapsFinished
  path "newFastas.fasta", emit: newHaps_singleFile, optional: true

  script:
  """
  cat ${allNewHaps} > concatHaps.txt
  grep "^>" ${existingHaps} > existingHapNames.txt
  python3 ${params.pyscripts}/treeClustering_scripts/newHaps_trial.py -n concatHaps.txt  -e existingHapNames.txt -o "${nonJunc_NEWHAPS}"

  touch newFastas.fasta
  newHapFiles="\$(ls ${nonJunc_NEWHAPS}/*.fasta | wc -l | awk '{print\$1}')"
  echo \$newHapFiles > mySize
  if [ \$newHapFiles -gt 0 ]
  then
  cat ${nonJunc_NEWHAPS}/*.fasta >> newFastas.fasta
  fi
  echo "done" > newHapsDone.txt

  """
}



process nameNewHaps_sameLength {
  label 'bioinformaticProcessing'

  input:
  path(allNewHaps)
  path(existingHaps)

  output:
  path "newHapsDone.txt", emit: nonJunc_hapsFinished
  path "newFastas.fasta", emit: newHaps_singleFile, optional: true

  script:
  """
  cat ${allNewHaps} > concatHaps.txt
  grep "^>" ${existingHaps} > existingHapNames.txt
  python3 ${params.pyscripts}/treeClustering_scripts/newHaps_requireMatchHapLength.py" -n concatHaps.txt -l "$baseDir/modules/files_forNextflow/CDC_8markers_hapLength.txt" -e existingHapNames.txt -o "${nonJunc_NEWHAPS}"

  touch newFastas.fasta
  cat ${nonJunc_NEWHAPS}/*.fasta >> newFastas.fasta
  echo "done" > newHapsDone.txt

  """
}


process makeRefs_assignHaps2 {
  input:
  path(newHaps)

  output:
  path 'refHaps.fasta', emit: completeFile

  script:
  """
  cat ${newHaps} > refHaps.fasta
  """
}

process makeRefs_assignHaps2_v2 {
  input:
  // path(newHaps)

  output:
  path 'refHaps.fasta', emit: completeFile

  script:
  """
  cat ${NEW_HAPS_DIR}/*.fasta > refHaps.fasta
  """
}

process finalNewRefs_assignHaps2 {
  publishDir "${refHapsDir}", pattern: "hapRefs.fasta", mode: "copy"

  input:
  path(originalNew)
  path(allNew)

  output:
  path('hapRefs.fasta')
  path 'hapRefs.fasta', emit: outFinalRefs

  script:
  """
  cat ${allNew} ${originalNew} ${params.originalHaplotypes_refs} > hapRefs.fasta
  """
}


process blastDB {
  label 'mappingSoftware'

  publishDir "${refHapsDir}", pattern: "*.fast*", mode: "copy"
 
  input:
  path(haps)

  output:
  path("*fast*")
  path "*fast*", emit: outFinalRefs_db

  script:
  """
  makeblastdb -in ${haps} -input_type fasta -dbtype nucl
  """
}


/*
The firstClustering and blastMulti steps have been updated to v2 to reduce processing times. 
The original blastMulti process is a major bottleneck when there many markers being analyzed and/or when the sample FASTQ has many reads (e.g., >250k reads).
The orignial blastMulti process used simple loops to search through every part of every marker and perform 2 blast steps. 
The first blast step takes the cluster sequences from 'firstClustering' and blasts those against all reference haplotypes. If there is a hit to any ref haplotype, then proceed to the second blast step
The second blast step is to blast all of the specimen's reads against the sequence that was used as the query in the first step. This will find the depth of the specimen reads mapping to the clustered read that hit a reference.
The depth of the hit is used to see if it passes the dynamic cutoff calculated in the BED extraction steps.

However, the depth determined by the second blast step is actually already calculated (albeit somewhat hidden).
The CD-HIT step in the firstClustering step clusters all identical reads together from all of that specimen's reads. Then the .clstr file has how many reads fell into that 100% identity cluster.
We just need a simple python script to extract the number of reads within each cluster in the clstr file.
This depth is the same as the depth calculated in the 2nd blast step, because the second blast step is a blast of all reads against each CD-HIT cluster at 100% identity. So that 2nd blast step is basically just redoing the clustering step.

How to fix this in the v2 version of firstClustering and blastMulti

1) In the firstClustering_v2 process - perform the same cd-hit clustering and perl make multiseq
2) quickly reformat clstr file and run through python script to count number of reads within each cluster 
3) Go through the multiseq folder and rename each sequence so they can be more easily pulled out to match cluster
4) Put depth for each cluster into a single file per cluster (rather than all depths into a single file)
5) Output 2 tuples with sampleName: sampleName.depth, and sampleName: sampleName.multiSeq. Putting the sampleName as the key of the tuple allows tracking outputs by sample and guarentees all outputs within a sample are linked

6) Concatanate all multiseq files per sample into a single file per sample. Can do this using the tuple input variable rather than loops
7) Blast the all multiseq file per sample against references in a single blast, rather than individual blast for each multiseq file
8) Format blast hits, multiseq, and depth files so that they are merged by Cluster ID, then renamed so that the Seq_ID is prepended to the haplotype name that the clusters hit to.
  Once again, there is a single file with all hits per sample, rather than processes each haplotype per sample individually in loops
9)Then look at the depth for the sample + marker dynamic cutoff that was calculated in a previous step and see if the depth from the .clstr file exceeds the cutoffs. If it exceeds cutoffs, then multiseq sequence for that cluster to a fasta file
10) This fasta file is the same ultimate output from the orginal blastMulti step and it goes in to the hapsBlast process
11) The hapsBlast process is necessary becaues that is a blast of the reference haplotypes against the sample fasta. So this guarantees we have full length haplotypes for our new specimens


*/

process firstClustering {
  label 'mappingSoftware'

  publishDir "${renamedReads}", pattern: "RENAMED*fasta", mode: "copy"

  input:
  path(clippedReads)

  output:
  path("RENAMED*fasta"), optional: true
  path '*multi-seq', emit: multiSeq_out, optional: true

  shell:
  '''
  cat !{clippedReads} | awk '{if(NR%4==1) $0=sprintf("@Read_%d",(1+i++)); print;}' > RENAMED.!{clippedReads.simpleName}.clipped_for_hap_calling.fastq
  
  clipLen="\$(wc -l RENAMED.!{clippedReads.simpleName}.clipped_for_hap_calling.fastq | awk '{print$1}')"
  if [ $clipLen -gt 0 ]
  then
  seqret -sequence RENAMED.!{clippedReads.simpleName}.clipped_for_hap_calling.fastq -outseq RENAMED.!{clippedReads.simpleName}.clipped_for_hap_calling.fasta

  cd-hit-est -i RENAMED.!{clippedReads.simpleName}.clipped_for_hap_calling.fastq -o !{clippedReads.simpleName}.clean_merged_Clusters.fq -c 1 -g 1 -d 0 -T !{params.threads} -s 1 -M !{params.RAM}

  seqret -sequence !{clippedReads.simpleName}.clean_merged_Clusters.fq -outseq !{clippedReads.simpleName}.clean_merged_Clusters.fasta

  perl !{params.multiSeq_pl} !{clippedReads.simpleName}.clean_merged_Clusters.fasta !{clippedReads.simpleName}.clean_merged_Clusters.fq.clstr !{clippedReads.simpleName}.multi-seq !{params.assignHap_min}
  else
  touch RENAMED.!{clippedReads.simpleName}.clipped_for_hap_calling.fasta
  mkdir !{clippedReads.simpleName}.multi-seq
  fi
  '''
}

process firstClustering_v2{

  label 'mappingSoftware'

  publishDir "${renamedReads}", pattern: "RENAMED*fasta", mode: "copy"

  input:
  path(clippedReads)

  output:
  path("RENAMED*fasta"), emit: sampleFasta ,  optional: true
  // path '*multi-seq', emit: multiSeq_out, optional: true
  // path '${clippedReads.simpleName}.*.multiOut', emit: multiSeq_Ind, optional: true
  // path '*toCount.fq.clstr', emit:clusterDepth, optional: true
  // path '${clippedReads.simpleName}.*.depth', emit:clusterDepthSingle, optional: true
  tuple val("${clippedReads.simpleName}"), path("${clippedReads.simpleName}.*.depth"), emit: depthTuple, optional: true
  tuple val("${clippedReads.simpleName}"), path("${clippedReads.simpleName}.*.multiOut"), emit: multiTuple, optional: true


  shell:
  '''
  cat !{clippedReads} | awk '{if(NR%4==1) $0=sprintf("@Read_%d",(1+i++)); print;}' > RENAMED.!{clippedReads.simpleName}.clipped_for_hap_calling.fastq
  
  clipLen="\$(wc -l RENAMED.!{clippedReads.simpleName}.clipped_for_hap_calling.fastq | awk '{print$1}')"
  if [ $clipLen -gt 0 ]
  then
  seqret -sequence RENAMED.!{clippedReads.simpleName}.clipped_for_hap_calling.fastq -outseq RENAMED.!{clippedReads.simpleName}.clipped_for_hap_calling.fasta

  cd-hit-est -i RENAMED.!{clippedReads.simpleName}.clipped_for_hap_calling.fastq -o !{clippedReads.simpleName}.clean_merged_Clusters.fq -c 1 -g 1 -d 0 -T !{params.threads} -s 1 -M !{params.RAM}

  seqret -sequence !{clippedReads.simpleName}.clean_merged_Clusters.fq -outseq !{clippedReads.simpleName}.clean_merged_Clusters.fasta

  perl !{params.multiSeq_pl} !{clippedReads.simpleName}.clean_merged_Clusters.fasta !{clippedReads.simpleName}.clean_merged_Clusters.fq.clstr !{clippedReads.simpleName}.multi-seq !{params.assignHap_min}
  
  else
  touch RENAMED.!{clippedReads.simpleName}.clipped_for_hap_calling.fasta
  mkdir !{clippedReads.simpleName}.multi-seq
  touch !{clippedReads.simpleName}.clean_merged_Clusters.fq.clstr 
  fi

  ls !{clippedReads.simpleName}.multi-seq > !{clippedReads.simpleName}.list

  cat !{clippedReads.simpleName}.clean_merged_Clusters.fq.clstr | sed 's/Cluster /Cluster/g' > !{clippedReads.simpleName}.toCount.fq.clstr

  clusLen="\$(wc -l !{clippedReads.simpleName}.toCount.fq.clstr | awk '{print$1}')"

  if [ $clusLen -gt 0 ]
  then
  python3 !{params.pyscripts}/treeClustering_scripts/cluster_counter.py -x !{clippedReads.simpleName}.toCount.fq.clstr -y !{clippedReads.simpleName}.clusterCounted.csv

  cat !{clippedReads.simpleName}.list | while read j
  do
  cp !{clippedReads.simpleName}.multi-seq/$j !{clippedReads.simpleName}.Cluster$j.trial
  cat !{clippedReads.simpleName}.Cluster$j.trial | sed "s/>Read/>Cluster${j}-Read/g" > !{clippedReads.simpleName}.Cluster$j.multiOut
  grep "Cluster$j," !{clippedReads.simpleName}.clusterCounted.csv > !{clippedReads.simpleName}.Cluster$j.depth
  done
  fi

  '''

}

// Dynamic cutoff comes from the finding new Hap SNPS script.
// To assign haps, we need depth above the set parameter valued AND the dynamic cutoff 
// This process also has a couple of bash loops. There's probably an elegant nextflow way to write this process but I couldn't get it running correctl using purely nextflow channels, etc., so i went with bash loops instead
// I had to add the path(refDB) to take in the blastDB process output channel otherwise the script would start before the ref db was created
process blastMulti {
  label 'mappingSoftware'

  publishDir "${nonJunc_mapping}", pattern: "*hapCov", mode: "copy"

  input:
  each(myMulti)
  path(refDB)

  output:
  path "*finalBlast.fasta", emit: finalHaps_blast, optional: true
  path "*hapCov", optional: true

  script:
  """
  touch ${myMulti.simpleName}.finalBlast.fasta

  ls ${myMulti} > ${myMulti}.list

  cat ${myMulti}.list | while read i
  do 
  blastn -db "${refHapsDir}/hapRefs.fasta" -query ${myMulti}/\$i -evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -num_threads ${params.threads} -max_target_seqs 1 -word_size 7 -dust no -soft_masking false  -outfmt "6 sseqid" -out ${myMulti.simpleName}.\$i.blastout
  
  hit="\$(wc -l ${myMulti.simpleName}.\$i.blastout | awk '{print\$1}')"
  if [ \$hit -gt 0 ]
  then
  tmp_var="\$(cat ${myMulti.simpleName}.\$i.blastout | sed 's/_Hap_.*//')"
  echo \$tmp_var >> ${myMulti.simpleName}.\$i.testHit
  cat ${myMulti}/\$i | sed "1s/.*/>${myMulti.simpleName}.\$tmp_var/" > ${myMulti.simpleName}.\$tmp_var.forBlast.fasta

  makeblastdb -in ${myMulti.simpleName}.\$tmp_var.forBlast.fasta -input_type fasta -dbtype nucl
  
  blastn -db ${myMulti.simpleName}.\$tmp_var.forBlast.fasta -query "${renamedReads}/RENAMED.${myMulti.simpleName}.clipped_for_hap_calling.fasta" -evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -num_threads ${params.threads} -out ${myMulti.simpleName}.\$tmp_var.estimatedCov -max_target_seqs 1 -word_size 7 -dust no -soft_masking false -outfmt "6 qseqid"
  
  cov="\$(wc -l ${myMulti.simpleName}.\$tmp_var.estimatedCov | awk '{print\$1}')"
  tmp_var2="\$(cat ${myMulti.simpleName}.\$i.blastout)"
  echo \$cov > ${myMulti.simpleName}.\$tmp_var2.hapCov

  dynamicCutoff=`cat ${nonJunc_mapping}/${myMulti.simpleName}.\$tmp_var.assignHaps_coverageCutoff`

  if [ \$cov -gt ${params.assignHap_min} ] && [ \$cov -gt \$dynamicCutoff ]
  then 
  cat ${myMulti.simpleName}.\$tmp_var.forBlast.fasta >> ${myMulti.simpleName}.finalBlast.fasta
  fi
  fi
  done
  """
}

process blastMulti_v2{

  label 'mappingSoftware'

  // publishDir "${newMethGenotypes}", pattern: "*haplotypes", mode: "copy"

  input:

  tuple val(seq_ID), path(multiFile)
  tuple val(seq_ID2), path(depthFile)
  each(refDB)
  // path(mySamp)

  output:
  // path("*haplotypes")
  path("${seq_ID}.finalBlast.fasta"), emit: blastFasta, optional: true
  // tuple val("${indMulti.baseName}"), path("${indMulti.baseName}.Num.fasta"), emit:testTuple

  shell:
  '''
  
  echo !{multiFile} >> myFiles

  cat myFiles | tr " " "\n" | awk -F "." '{print$1}'  > test

  cat !{multiFile} >> !{seq_ID}_hapsBlast.fasta

  blastn -db !{refHapsDir}/hapRefs.fasta -query !{seq_ID}_hapsBlast.fasta -evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -num_threads !{params.threads} -max_target_seqs 1 -word_size 7 -dust no -soft_masking false  -outfmt "6 sseqid qseqid" -out !{seq_ID}.blastout

  cat !{seq_ID}.blastout | sed 's/_Read/-Read/g' | awk -F "-" '{print$1"\t"$2}' |awk -F "\t" '{print$1"\t"$2"\t"$3}' | sort -k 2 > !{seq_ID}.blastFormat

  cat !{depthFile} | awk -F "," '{print$1"\t"$2}' | sort -k 1 >> !{seq_ID}.myDepth

  awk -F "\t" '{print$2}' !{seq_ID}.blastFormat | while read line; do grep $line -w -h !{seq_ID}.myDepth !{seq_ID}.blastFormat | tr "\n" "\t" ; echo ; done | awk -F "\t" '{print$3"\t"$1"\t"$2"\t"$4"\t"$5}' > !{seq_ID}.hapDepth  

  mySamp=!{seq_ID}

  awk -v var=!{seq_ID} -F "_" '{print$0, "\t"var"."$1"_"$2"_"$3"_"$4"_"$5}'  !{seq_ID}.hapDepth > tester

  awk -F "\t" '{print$6"\t"$3"\t"$1"\t"$4"\t"$4"-"$5}' tester | while read search actDepth hap clusterName readName
  do
  FILE=!{nonJunc_mapping}/$search.assignHaps_coverageCutoff
  echo "$FILE"
  if test -f "$FILE";
  then
  dynamicCutoff=`cat !{nonJunc_mapping}/$search.assignHaps_coverageCutoff`
  echo $dynamicCutoff >> sampleMarkerDepth
  echo !{seq_ID}.$clusterName.multiOut >> multiNameTesting
  else
  dynamicCutoff=!{params.assignHap_min}
  echo "no file exists"
  fi
  if [ $actDepth -gt $dynamicCutoff ] && [ $actDepth -gt !{params.assignHap_min} ]
  then 
  cat !{seq_ID}.$clusterName.multiOut >> !{seq_ID}.finalBlast.fasta
  echo $hap >> !{seq_ID}.haplotypes
  fi
  done



  '''


}
process hapsBlast {
  label 'mappingSoftware'
  
  publishDir "${tmpGenotypes}", pattern: "*GENOTYPE_noJunction", mode: "copy"

  input:
  each(myHaps)
  path(refDB)

  output:
  path "*GENOTYPE_noJunction", emit: genotypeNoJunction
  path("*GENOTYPE_noJunction")
  path("endMarker"), emit: SNP_hapsDone

  script:
  """
  hapLen="\$(wc -l ${myHaps} | awk '{print\$1}')"

  if [ \$hapLen -gt 0 ]
  then

  makeblastdb -in ${myHaps} -input_type fasta -dbtype nucl

  blastn -db ${myHaps} -query "${refHapsDir}/hapRefs.fasta" -evalue 0.001 -perc_identity 100 -qcov_hsp_perc 100 -num_threads ${params.threads} -out ${myHaps.simpleName}.GENOTYPE_noJunction -max_target_seqs 1 -word_size 7 -dust no -soft_masking false -outfmt "6 qseqid pident mismatch gapopen gaps sseq evalue bitscore"
  else
  touch ${myHaps.simpleName}.GENOTYPE_noJunction
  fi

  echo "done" >> endMarker
  """
}
