process CHECK_ME_BBDUK {

    tag "${sample_id}"

    publishDir "${params.outdir}/check_me/${sample_id}", mode: 'copy'


    input:
    tuple val(sample_id), path(fastq)
 

    output:
    tuple val(sample_id), path("${sample_id}.matched.fastq.gz"),   emit: matched
    tuple val(sample_id), path("${sample_id}.unmatched.fastq.gz"), emit: unmatched
    path  "${sample_id}.stats.2_error.txt",                         emit: stats
    path  "${sample_id}.log",                                       emit: log

    script:
    """
    bbduk.sh -Xmx${task.memory.toGiga()}g \\
        in=${fastq} \\
        ref=${params.me_sequence} \\
        out=${sample_id}.unmatched.fastq \\
        outm=${sample_id}.matched.fastq \\
        k=19 \\
        hdist=2 \\
        maskmiddle=f \\
        overwrite=t \\
        stats=${sample_id}.stats.2_error.txt \\
        2>> ${sample_id}.log

    gzip ${sample_id}.matched.fastq
    gzip ${sample_id}.unmatched.fastq
    """
}
