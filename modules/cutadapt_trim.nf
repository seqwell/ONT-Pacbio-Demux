
process CUTADAPT_TRIM {
    tag "${pair_id}"
    publishDir path: "${params.outdir}/demuxed_fastq_flat/", mode: 'copy', pattern: "*seqWell*"
    publishDir path: "${params.outdir}/demuxed_fastq/", mode: 'copy', saveAs: { filename ->
        if (filename.contains("seqWell")) "${pair_id}/${filename}"
        else if (filename == "unknown.fastq.gz") filename
        else null
    }
    publishDir path: "${params.outdir}/demuxed_fastq", mode: 'copy', pattern: "unknown.fastq.gz"
    publishDir path: "${params.outdir}/other/ME_tagged_fastq", mode: 'copy', pattern: '*tag*'

    input:
    tuple val(pair_id), path(fq)

    output:
    path("*.ME.tagged.fastq.gz")
    path("{*seqWell,*.ME.tagged}.fastq.gz"), emit: fq

    script:
    """
    # Step 1: Adapter trim
    cutadapt \\
        -g AGATGTGTATAAGAGACAG \\
        -g AATGATACGGCGACCACCGAGATCTACAC \\
        --minimum-length ${params.length_filter} \\
        -e 0.106 \\
        -O 12 \\
        --report=minimal \\
        -o ${pair_id}.step1.fastq.gz \\
        ${fq} > ${pair_id}.cutadapt_report.step1.txt

    # Step 2: ME tag detection
    cutadapt \\
        -b AGATGTGTATAAGAGACAG \\
        -b CTGTCTCTTATACACATCT \\
        -e 0 \\
        -O 19 \\
        --action=none \\
        --report=minimal \\
        --untrimmed-output ${pair_id}.seqWell.fastq.gz \\
        -o ${pair_id}.ME.tagged.fastq.gz \\
        ${pair_id}.step1.fastq.gz > ${pair_id}.cutadapt_report.txt

    rm *step*.fastq.gz
    """
}
