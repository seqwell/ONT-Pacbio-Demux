process COMBINE_BARCODES {
    tag "combine_barcodes"

    // publishDir "${params.outdir}/merged_fastq", mode: 'copy'

    input:
    path dummy  // collected fastq_dir paths

    output:
    path("output_fastq/*"), emit: fastq

    script:
    """
    #!/bin/bash
    set -euo pipefail
    shopt -s nullglob
    mkdir -p output_fastq

    dirs=( *_Demux_fastq )
    (( \${#dirs[@]} )) || { echo "No Demux_fastq directories found" >&2; exit 1; }

    for i in \$(seq -w 001 384); do
        files=()
        for d in "\${dirs[@]}"; do
            for f in "\$d"/*_sw_barcode\${i}.fastq; do files+=("\$f"); done
        done
        if (( \${#files[@]} )); then
            cat "\${files[@]}" | gzip > "output_fastq/barcode\${i}.fastq.gz"
        fi
    done

    files=()
    for d in "\${dirs[@]}"; do
        for f in "\$d"/*_unclassified.fastq; do files+=("\$f"); done
    done
    (( \${#files[@]} )) && cat "\${files[@]}" | gzip > "output_fastq/unclassified.fastq.gz" || true
    """
}

process COMBINE_BARCODES_BAM {
    tag "combine_barcodes_bam"

    input:
    path dummy  // collected bam_dir paths

    output:
    path("output_bam/*"), emit: bam

    container 'quay.io/biocontainers/samtools:1.19--h50ea8bc_0'

    script:
    """
    #!/bin/bash
    set -euo pipefail
    shopt -s nullglob
    mkdir -p output_bam

    dirs=( *_Demux_bam )
    (( \${#dirs[@]} )) || { echo "No Demux_bam directories found" >&2; exit 1; }

    for i in \$(seq -w 001 384); do
        files=()
        for d in "\${dirs[@]}"; do
            for f in "\$d"/*_sw_barcode\${i}.bam; do files+=("\$f"); done
        done
        if (( \${#files[@]} )); then
            samtools merge -f "output_bam/barcode\${i}.bam" "\${files[@]}"
            samtools index "output_bam/barcode\${i}.bam"
        fi
    done

    files=()
    for d in "\${dirs[@]}"; do
        for f in "\$d"/*_unclassified.bam; do files+=("\$f"); done
    done
    if (( \${#files[@]} )); then
        samtools merge -f "output_bam/unclassified.bam" "\${files[@]}"
        samtools index "output_bam/unclassified.bam"
    fi
    """
}
