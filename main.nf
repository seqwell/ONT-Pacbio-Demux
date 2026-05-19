nextflow.enable.dsl = 2


include { DORADO_DEMUX         } from './modules/dorado_demux.nf'
include { COMBINE_BARCODES     } from './modules/combine_barcodes.nf'
include { COMBINE_BARCODES_BAM } from './modules/combine_barcodes.nf'
include { BAM_TO_FASTQ         } from './modules/bam_to_fastq.nf'
include { EXTRACT_HEADER       } from './modules/extract_header.nf'
include { REHEADER_READS       } from './modules/reheader_reads.nf'
include { CUTADAPT_TRIM        } from './modules/cutadapt_trim.nf'
include { FILTER_BAM_BY_FASTQ  } from './modules/filter_bam_by_fastq.nf'
include { DEMUX_SUMMARIZE      } from './modules/demux_summarize.nf'
include { NANOSTAT             } from './modules/nanostat.nf'
include { MULTIQC              } from './modules/multiQC.nf'
include { READ_LENGTH          } from './modules/read_length.nf'


workflow {

    // ---------------------------------------------------------------
    // Validate data_type param
    // ---------------------------------------------------------------
    if (!params.data_type) {
        error "Please specify --data_type [bam|fastq]"
    }
    if (!['bam','fastq'].contains(params.data_type)) {
        error "Invalid --data_type '${params.data_type}'. Must be 'bam' or 'fastq'."
    }

    use_bam = (params.data_type == 'bam')

    if (use_bam) {
        log.info "data_type=bam  →  BAM input: dorado demux (BAM) → BAM-to-FASTQ → cutadapt → filter demux-BAM by cutadapt-FASTQ read IDs"
        input_ch = Channel.fromPath(params.input + "/*.bam")
                     .map { it -> tuple(it.baseName, it) }
    } else {
        log.info "data_type=fastq  →  FASTQ input: extract headers → dorado demux → reheader → cutadapt → FASTQ output only (no BAM created)"
        input_ch = Channel.fromPath(params.input + "/*.fastq.gz")
                     .map { it -> tuple(it.baseName.replace(".fastq", ""), it) }
    }

    barcode_fasta    = file(params.barcodes)
    arrangement_toml = file(params.arrangement_toml)


    // ---------------------------------------------------------------
    // Step 1: Demux
    //   bam input   → dorado emits BAM dir only (--no emit-fastq)
    //   fastq input → dorado emits FASTQ dir only (--emit-fastq)
    // ---------------------------------------------------------------
    DORADO_DEMUX(input_ch, barcode_fasta, arrangement_toml)


    // ================================================================
    //  BAM INPUT BRANCH
    // ================================================================
    if (use_bam) {

        // -----------------------------------------------------------
        // Step 2 (BAM): Combine per-sample BAM demux dirs → per-barcode BAMs
        // -----------------------------------------------------------
        COMBINE_BARCODES_BAM(DORADO_DEMUX.out.bam_dir.collect())

        bam_barcode_ch = COMBINE_BARCODES_BAM.out.bam
                           .flatMap { it instanceof List ? it : [it] }
                           .filter { it.name.endsWith('.bam') }
                           .map { bam -> tuple(bam.baseName, bam) }

        // -----------------------------------------------------------
        // Step 3 (BAM): Convert each per-barcode demuxed BAM → FASTQ
        //   This FASTQ is used as input to cutadapt.
        //   No header extraction is needed for BAM input.
        // -----------------------------------------------------------
        BAM_TO_FASTQ(bam_barcode_ch)

        // -----------------------------------------------------------
        // Step 4 (BAM): Cutadapt trim on the converted FASTQ
        // -----------------------------------------------------------
        CUTADAPT_TRIM(BAM_TO_FASTQ.out.fastq)

        // trimmed_ch keyed by bare barcode ID (strip suffix after first dot)
        // e.g. "barcode001.seqWell.fastq.gz" → key "barcode001"
        // This must match bam_barcode_ch which is also keyed by "barcode001"
        base_trimmed_ch = CUTADAPT_TRIM.out.fq
                          .flatMap { it instanceof List ? it : [it] }
                          .filter { fq -> !fq.name.contains('tagged') && fq.size() > 20 }

        trimmed_ch = base_trimmed_ch
                      .map { fq -> tuple(fq.baseName.replace('.seqWell.fastq',''), fq) }

        trimmed_ch_for_nanostat = base_trimmed_ch
                      .map { fq -> tuple(fq.baseName.replace('.fastq',''), fq) }



        // -----------------------------------------------------------
        // Step 5 (BAM): Split the demux-BAM to match cutadapt-FASTQ results
        //   Join on bare barcode ID so the keys match on both sides:
        //     bam_barcode_ch key: "barcode001"
        //     trimmed_ch key:     "barcode001"
        //   Then subset the BAM to only reads whose names appear in the FASTQ,
        //   propagating cutadapt adapter/length/ME-tag filtering onto the BAM.
        // -----------------------------------------------------------
  
        bam_fastq_ch = bam_barcode_ch.join(trimmed_ch, by: 0)
                         // produces: tuple(barcode_id, bam, fastq)
   

        FILTER_BAM_BY_FASTQ(bam_fastq_ch)

        demuxed_ch = trimmed_ch_for_nanostat

    // ================================================================
    //  FASTQ INPUT BRANCH
    // ================================================================
    } else {

        // -----------------------------------------------------------
        // Step 2 (FASTQ): Extract FASTQ headers BEFORE demuxing
        //   awk extracts FASTQ header fields → uuid_tags.tsv
        //   Consumed by REHEADER_READS after demux to restore original headers.
        // -----------------------------------------------------------
        EXTRACT_HEADER(input_ch)
        merged_tags = EXTRACT_HEADER.out.collect()

        // -----------------------------------------------------------
        // Step 3 (FASTQ): Combine per-sample FASTQ demux dirs → per-barcode FASTQs
        // -----------------------------------------------------------
        COMBINE_BARCODES(DORADO_DEMUX.out.fastq_dir.collect())

        fastq_barcode_ch = COMBINE_BARCODES.out.fastq
                             .flatMap { it instanceof List ? it : [it] }
                             .map { fq -> tuple(fq.baseName.replace(".fastq", ""), fq) }

        // -----------------------------------------------------------
        // Step 4 (FASTQ): Restore original header tags into demuxed FASTQs
        // -----------------------------------------------------------
        REHEADER_READS(fastq_barcode_ch, merged_tags)

        // -----------------------------------------------------------
        // Step 5 (FASTQ): Cutadapt trim — on reheadered FASTQ
        //   No BAM files are created for FASTQ input.
        // -----------------------------------------------------------
        CUTADAPT_TRIM(REHEADER_READS.out.fq)

        trimmed_ch = CUTADAPT_TRIM.out.fq
                       .flatMap { it instanceof List ? it : [it] }
                       .map { fq -> tuple(fq.baseName.replace(".fastq", ""), fq) }

        demuxed_ch = trimmed_ch
    }


    // ---------------------------------------------------------------
    // Step 6: Downstream QC — runs on trimmed FASTQ (both branches)
    // ---------------------------------------------------------------
    DEMUX_SUMMARIZE(CUTADAPT_TRIM.out.fq.collect())

    valid_ids_ch = DEMUX_SUMMARIZE.out
                    .splitCsv()
                    .filter { row ->
                        row[0].startsWith("barcode") || row[0] == "unknown"
                    }
                    .map { row -> tuple(row[0].trim(), true) }

    filtered_ch = demuxed_ch
                    .join(valid_ids_ch, by: 0)
                    .map { sample_id, fq, flag -> tuple(sample_id, fq) }
                    

    READ_LENGTH(filtered_ch)
    NANOSTAT(filtered_ch)
    MULTIQC(NANOSTAT.out.collect())
}
