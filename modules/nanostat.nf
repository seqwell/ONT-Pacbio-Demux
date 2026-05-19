
process NANOSTAT {
    tag "$sample_id"
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(fastq)

    output:
    path("${sample_id}_nanostat.txt")

    script:
    """
    NanoStat --fastq ${fastq} --name ${sample_id}_nanostat.txt --threads ${task.cpus}
    """
}
