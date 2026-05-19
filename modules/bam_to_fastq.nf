/*
 * BAM_TO_FASTQ
 *
 * For BAM-input runs only.
 * Converts a per-barcode demuxed BAM (from COMBINE_BARCODES_BAM) into a
 * gzipped FASTQ so that CUTADAPT_TRIM can run on it.
 * The resulting FASTQ is then used both for cutadapt filtering AND for
 * determining which reads to keep in the final BAM (via FILTER_BAM_BY_FASTQ).
 */
process BAM_TO_FASTQ {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz"), emit: fastq

    container 'quay.io/biocontainers/samtools:1.19--h50ea8bc_0'

    script:
    """
    samtools fastq -T '*' ${bam} | gzip > ${sample_id}.fastq.gz
    """
}
