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

