process DEMUX_SUMMARIZE {
    publishDir path: "${params.outdir}/demux_summary", mode: 'copy', pattern: '*.csv'
       
    input:
     path(fq_files)

    output:
     path("*.csv")

    script:
    """
    summarize_barcodes.sh "${params.pool_ID}" "${params.error_rate}"
    """
}

