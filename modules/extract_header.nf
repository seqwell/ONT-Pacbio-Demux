/*
 * EXTRACT_HEADER
 *
 * For FASTQ-input runs only.
 * Extracts header fields from the input FASTQ and writes them to a TSV
 * (uuid_tags.tsv) so that REHEADER_READS can restore the original header
 * tags after dorado demux strips them.
 *
 * Not used for BAM input — BAM aux tags are preserved end-to-end in the
 * demuxed BAM, and no header extraction step is required.
 */
process EXTRACT_HEADER {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(input_file)

    output:
    path "${sample_id}.uuid_tags.tsv"

    script:
    """
    extract_fastq_tags.py ${input_file} ${sample_id}.uuid_tags.tsv
    """
}
