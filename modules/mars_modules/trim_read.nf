// Adapter and quality trimming
process Trim_reads {
    label 'nfNest'
    errorStrategy 'ignore'

    tag { "Trim_reads${sample_id}"}

    publishDir "${params.out}/Trimmed_Fastq/${sample_id}", pattern: "*.fastq", mode : "copy"
    publishDir "${params.out}/Trimmed_Fastq/${sample_id}/Stats", pattern: "*.txt", mode : "copy"

    input:
    tuple path (adapter), val(sample_id), path(reads)



    output:
    tuple val(sample_id), path("*.bbduk.fastq"), emit: Trimmed_fastq
    // file "${sample_id}.stats.txt"
    tuple val(sample_id), path("${sample_id}.stats.txt"), emit: Trimmed_stats

    path "*.clean_merged.fastq", emit: readyReads
    path "*GENOTYPE*"


    script:


     """

     bbduk.sh  -Xmx1g  ktrimright=t k=27 hdist=1 edist=0 ref=${adapter} \\
     qtrim=rl trimq=30  minlength=50 trimbyoverlap=t minoverlap=24 qin=33 in=${reads[0]} in2=${reads[1]} \\
     out=${sample_id}.R1.bbduk.fastq out2=${sample_id}.R2.bbduk.fastq stats=${sample_id}.stats.txt


    newName="\$(echo ${sample_id} | cut -c 1-${params.nameLength} | sed 's/-/_/g')"



    cat ${sample_id}.R1.bbduk.fastq ${sample_id}.R2.bbduk.fastq > \$newName.clean_merged.fastq

    touch \$newName.GENOTYPE_noJunction
    touch \$newName.GENOTYPE_withJunction
    touch \$newName.GENOTYPE_final

    """


}
