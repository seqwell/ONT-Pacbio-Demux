process FILTER_BAM_BY_FASTQ {
    tag "$sample_id"
    publishDir path: "${params.outdir}/demuxed_bam_flat/",        mode: 'copy', pattern: "*seqWell*.bam"
    publishDir path: "${params.outdir}/demuxed_bam/${sample_id}/", mode: 'copy', pattern: "*seqWell*.bam"
    publishDir path: "${params.outdir}/other/ME_tagged_bam/",     mode: 'copy', pattern: "*tagged.bam"
    publishDir path: "${params.outdir}/demuxed_bam/",             mode: 'copy', pattern: "unknown.bam"

    input:
    tuple val(sample_id), path(bam), path(untagged_fastq)
    // untagged_fastq is e.g. barcode130.untagged.fastq.gz

    output:
    tuple val(sample_id), path("*tagged.bam"),  emit: tagged_bam,  optional: true
    tuple val(sample_id), path("*seqWell.bam"), emit: seqwell_bam

    script:
    // barcode130.untagged.fastq.gz → base = barcode130.untagged
    // seqWell BAM: barcode130.seqWell.bam  (reads IN untagged FASTQ)
    // tagged BAM:  barcode130.ME.tagged.bam (reads NOT IN untagged FASTQ)
    def base        = untagged_fastq.name.replaceAll(/\.fastq\.gz$/, '')   // barcode192.seqWell
    def prefix      = base.replaceAll(/\.seqWell$/, '')                    // barcode192
    def seqwell_bam = "${prefix}.seqWell.bam"                              // barcode192.seqWell.bam
    def tagged_bam  = "${prefix}.ME.tagged.bam"                            // barcode192.ME.tagged.bam
    """
    # Step 1: Extract read IDs from the untagged FASTQ
    if [[ "${untagged_fastq}" == *.gz ]]; then
        zcat "${untagged_fastq}" | awk 'NR%4==1 {sub(/^@/,""); split(\$0,a," "); print a[1]}' > untagged_reads.txt
    else
        awk 'NR%4==1 {sub(/^@/,""); split(\$0,a," "); print a[1]}' "${untagged_fastq}" > untagged_reads.txt
    fi

    # Step 2: Reads IN untagged FASTQ → seqWell BAM (normal/untagged reads)
    samtools view -N untagged_reads.txt -o "${seqwell_bam}" "${bam}"
    samtools index "${seqwell_bam}"

    # Step 3: Reads NOT IN untagged FASTQ → tagged BAM (ME-tagged reads)
    # -U writes reads that do NOT pass the filter
    samtools view -N untagged_reads.txt -U "${tagged_bam}" -o /dev/null "${bam}"
    samtools index "${tagged_bam}"
    """
}