process PreFastqC {
     label 'nfNest'

    tag { "PreFastqC ${reads}"}


    publishDir "${params.out}/PrefastqC/ZIP", pattern: "*.zip",  mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path ("*_fastqc.zip")

    script:
    """
    fastqc ${reads}

    """

}

process pre_multiqc{
     label 'nfNest'


     publishDir "${params.out}/multiqc1/", pattern: "*.html",  mode:'copy'


     input:
     path ("*")

     output:
     path("multiqc_report.html")

     script:
     """

      multiqc .

     """

}

process PostFastqC{
     label 'nfNest'

    //tag { "PostFastqC ${reads}"}

    publishDir "${params.out}/PostfastqC/",  mode:'copy'


    input:
    tuple val(sample_id), path(reads)

    output:
    path ("*_fastqc.zip")

    script:
    """

    fastqc ${reads}

    """

}

process multiqc{
     label 'nfNest'


     publishDir "${params.out}/multiqc/", pattern: "*.html",  mode:'copy'


     input:
     path ("*")

     output:
     path("multiqc_report.html")

     script:
     """

      multiqc .

     """

}
