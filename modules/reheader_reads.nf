process REHEADER_READS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(fastq)
    path tags_tsv

    output:
    tuple val(sample_id), path("${sample_id}.reheadered.fastq.gz"), emit: fq

    script:
    """
    reheader_reads.py \\
        ${fastq} \\
        ${sample_id}.reheadered.fastq.gz \\
        ${tags_tsv}
    """
}