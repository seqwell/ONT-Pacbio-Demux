process FASTQ_TO_BAM {
    tag "$sample_id"
    publishDir path: "${params.outdir}/demuxed_bam_flat/",       mode: 'copy', pattern: "*seqWell*.bam"
    publishDir path: "${params.outdir}/demuxed_bam/${sample_id}/",mode: 'copy', pattern: "*seqWell*.bam"
    publishDir path: "${params.outdir}/other/ME_tagged_bam/",   mode: 'copy', pattern: "*tagged.bam"
    publishDir path: "${params.outdir}/demuxed_bam/",            mode: 'copy', pattern: "unknown.bam"

    input:
    tuple val(sample_id), path(fastq)
    path tags_tsv

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam

    script:
    """
    fastq_to_bam.py ${fastq} ${sample_id}.bam ${tags_tsv}
    """
}
