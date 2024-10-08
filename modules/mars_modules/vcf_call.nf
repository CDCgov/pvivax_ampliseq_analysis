process Sam_sort {
  label 'nfNest'

    tag { "Sam_sort ${sample_id}"}
    // publishDir "${params.out}/bam_out", mode:'copy'

    // input is reads and refenace
    input:
    tuple val(sample_id), path(sam)

    // output is samfile
    output:
    tuple val (sample_id) ,path("${sample_id}_sorted.bam"), emit: Bam_file


    //script //  bwa mem ref.fa read1.fq read2.fq > aln-pe.sam
    script:
    """

    samtools view -S -b ${sample_id}_align.sam > ${sample_id}_align.bam
    samtools sort -o "${sample_id}_sorted.bam" ${sample_id}_align.bam
    """

}


process Picard_add_read{
   label 'nfNest'

    tag { "Picard_add_read${sample_id}"}
    // publishDir "${params.out}/picard_out", mode:'copy'


    input:
    tuple val (sample_id) ,path("${sample_id}_align.sam")

    output:
    tuple val (sample_id) ,path("${sample_id}_picard_readgroup.bam"), emit: Picard_out_bam


    script:

    """

    #picard AddOrReplaceReadGroups I=${sample_id}_align.sam O=${sample_id}_picard_readgroup.bam \\
    SORT_ORDER=coordinate RGLB=ExomeSeq  RGPL=Illumina RGPU=NextSeq  RGSM=${sample_id} CREATE_INDEX=True

    java -jar /usr/local/bin/picard/picard.jar AddOrReplaceReadGroups I=${sample_id}_align.sam O=${sample_id}_picard_readgroup.bam \\
    SORT_ORDER=coordinate RGLB=ExomeSeq  RGPL=Illumina RGPU=NextSeq  RGSM=${sample_id} CREATE_INDEX=True


    """

}

process Get_Bed {
   label 'nfNest'
publishDir "${params.out}/BED", mode:'copy'

 input:
  path(gff)

 output:
  file ('mars_genes.bed')


 script:

  """

   awk 'BEGIN { OFS="\t" } {if (\$3=="gene") {print \$1,\$4-1,\$5,\$10,\$16,\$7}}' ${gff}  > mars_genes.bed

  """
}

process VCF_call {
   label 'nfNest'
  tag  { "VCF_call ${sample_id }"}
  publishDir "${params.out}/vcf_files", mode:'copy'
  cpus params.threads
  memory params.memory


  input:
    tuple path(ref), path("*"), path ("mars_genes.bed"), val (sample_id), path("${sample_id}_picard_readgroup.bam")

  output:
    //tuple val (sample_id) ,path("${sample_id}.mpileup"),      emit: BCF_out
    //tuple val (sample_id) ,path("${sample_id}_bcf.vcf"),      emit: Variants_out
    tuple val (sample_id) ,path("${sample_id}_samtools.vcf"),path("${sample_id}_Freebayes.vcf"),  path("${sample_id}_gatk.vcf"), path("${sample_id}_vardict.vcf") , emit: variants



  script:

    """

    samtools index ${sample_id}_picard_readgroup.bam


    # Samtools
    bcftools mpileup -O b -o ${sample_id}.mpileup -f ${ref} ${sample_id}_picard_readgroup.bam
    bcftools call --ploidy 1 -m -v -o ${sample_id}_bcf.vcf ${sample_id}.mpileup
    vcfutils.pl varFilter ${sample_id}_bcf.vcf > ${sample_id}_samtools.vcf

    # freebayes (calls mnps and snps)
    #freebayes -f  ${ref} -F 0.05 -E 3 --report-all-haplotype-alleles --haplotype-length -1 ${sample_id}_picard_readgroup.bam > ${sample_id}_Freebayes.vcf
    freebayes -f  ${ref} -F 0.05 -E 3 --report-all-haplotype-alleles ${sample_id}_picard_readgroup.bam > ${sample_id}_Freebayes.vcf

    # GATK
    gatk HaplotypeCaller --native-pair-hmm-threads 8 -R ${ref} -I ${sample_id}_picard_readgroup.bam  --min-base-quality-score 0 -O ${sample_id}_gatk.vcf


    # vardict
    vardict -G ${ref}  -f 0.05 -N ${sample_id} -b ${sample_id}_picard_readgroup.bam -c 1 -S 2 -E 3 -g 4 'mars_genes.bed'| var2vcf_valid.pl -N ${sample_id} -E -f 0.05 >  ${sample_id}_vardict.vcf


    """

}
