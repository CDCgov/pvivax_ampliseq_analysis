process vcf_to_DF {
   label 'nfNest'
    publishDir "${params.out}/vcf_to_DF/${sample_id}", mode : "copy"


    input:
      tuple val(sample_id), path("${sample_id}_samtools_vartype.vcf"), path("${sample_id}_freeBayes_vartype.vcf"), path("${sample_id}_HaplotypeCaller_vartype.vcf"), path("${sample_id}_vardict_vartype.vcf")

    output:
     tuple val(sample_id), path("${sample_id}_samtools.csv"), path("${sample_id}_freeBayes.csv"), path("${sample_id}_HaplotypeCaller.csv"),path("${sample_id}_vardict.csv"), emit: csv_annotate

    script:
        """

        ${params.pyscripts}/mars_scripts/vcf_merge.py -n ${sample_id}_samtools_vartype.vcf
        ${params.pyscripts}/mars_scripts/vcf_merge.py -n ${sample_id}_freeBayes_vartype.vcf
        ${params.pyscripts}/mars_scripts/vcf_merge.py -n ${sample_id}_HaplotypeCaller_vartype.vcf
        ${params.pyscripts}/mars_scripts/vcf_merge.py -n ${sample_id}_vardict_vartype.vcf
        """
}

process csv_merge {
   label 'nfNest'
    publishDir "${params.out}/CSV_merge/${sample_id}", mode : "copy"


    input:
      tuple  val(sample_id), path("${sample_id}_samtools.csv"), path("${sample_id}_freeBayes.csv"), path("${sample_id}_HaplotypeCaller.csv"), path("${sample_id}_vardict.csv")


    output:
      tuple val(sample_id), path("${sample_id}_merge.csv"), emit: CSV_merge
      tuple val(sample_id), path("${sample_id}_introns.csv"), emit: Introns_file

    script:

        """


        ${params.pyscripts}/mars_scripts/csv_merge.py -v1 "${sample_id}_samtools.csv" -v2  '${sample_id}_freeBayes.csv' -v3 '${sample_id}_HaplotypeCaller.csv' -v4 '${sample_id}_vardict.csv'


        """
}
