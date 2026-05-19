process MULTIQC {
  
    publishDir path: "${params.outdir}/multiQC", mode: 'copy'
    
    input:
    path(nanoplot_dirs)

    output:
    path "*multiqc_report.html"

    script:
    """
    multiqc . -o multiqc_report 
    mv multiqc_report/multiqc_report.html ${params.pool_ID}_ONT_multiqc_report.html
    """
}

