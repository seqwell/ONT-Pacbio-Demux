include { CHECK_ME_BBDUK       } from './modules/check_me_bbduk.nf'
include { DORADO_DEMUX         } from './modules/dorado_demux.nf'
include { COMBINE_BARCODES     } from './modules/combine_barcodes.nf'
include { COMBINE_BARCODES_BAM } from './modules/combine_barcodes_bam.nf'
include { BAM_TO_FASTQ as  BAM_TO_FASTQ_for_bbduk        } from './modules/bam_to_fastq.nf'
include { BAM_TO_FASTQ         } from './modules/bam_to_fastq.nf'
include { EXTRACT_HEADER       } from './modules/extract_header.nf'
include { REHEADER_READS       } from './modules/reheader_reads.nf'
include { CUTADAPT_TRIM        } from './modules/cutadapt_trim.nf'
include { FILTER_BAM_BY_FASTQ  } from './modules/filter_bam_by_fastq.nf'
include { FILTER_BAM_BY_FASTQ  as  FILTER_BAM_BY_FASTQ_for_bbduk } from './modules/filter_bam_by_fastq.nf'
include { DEMUX_SUMMARIZE      } from './modules/demux_summarize.nf'
include { NANOSTAT             } from './modules/nanostat.nf'
include { MULTIQC              } from './modules/multiQC.nf'
include { READ_LENGTH          } from './modules/read_length.nf'
include { NANOFILT             } from './modules/nanofilt.nf'

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

    barcode_fasta    = file(params.barcodes)
    arrangement_toml = file(params.arrangement_toml)

    // ================================================================
    //  BAM INPUT BRANCH — Step 0: convert raw BAM → FASTQ → ME filter → filter BAM
    // ================================================================
    if (use_bam) {
        log.info "data_type=bam  →  BAM input: BAM-to-FASTQ → ME filter → filter raw BAM → dorado demux (BAM) → BAM-to-FASTQ → cutadapt → filter demux-BAM by cutadapt-FASTQ read IDs"

        raw_bam_ch = Channel.fromPath(params.input + "/*.bam")
                       .map { it -> tuple(it.baseName, it) }

        // -----------------------------------------------------------
        // Step 0a (BAM): Convert raw input BAM → FASTQ for ME check
        // -----------------------------------------------------------
        BAM_TO_FASTQ_for_bbduk(raw_bam_ch)

        // -----------------------------------------------------------
        // Step 0b (BAM): Check ME adapter presence on raw FASTQ
        // -----------------------------------------------------------
        CHECK_ME_BBDUK(BAM_TO_FASTQ_for_bbduk.out.fastq)
        
        NANOFILT(CHECK_ME_BBDUK.out.matched)

        // -----------------------------------------------------------
        // Step 0c (BAM): Filter raw BAM to only ME-matched nanofilt filtered read IDs
        //   Reuses FILTER_BAM_BY_FASTQ: tuple(sample_id, bam, matched_fastq)
        // -----------------------------------------------------------
      
        me_filter_bam_ch = raw_bam_ch.join(NANOFILT.out, by: 0)

        FILTER_BAM_BY_FASTQ_for_bbduk(me_filter_bam_ch)

        // ME-filtered BAM is the input to the rest of the BAM pipeline
        input_ch = FILTER_BAM_BY_FASTQ_for_bbduk.out.seqwell_bam

    // ================================================================
    //  FASTQ INPUT BRANCH — Step 0: ME filter → NANOFILT
    // ================================================================
    } else {
        log.info "data_type=fastq  →  FASTQ input: ME filter → nanofilt → extract headers → dorado demux → reheader → cutadapt → FASTQ output only (no BAM created)"

        raw_fastq_ch = Channel.fromPath(params.input + "/*.fastq.gz")
                         .map { it -> tuple(it.baseName.replace(".fastq", ""), it) }

        // -----------------------------------------------------------
        // Step 0a (FASTQ): Check ME adapter presence on raw FASTQ
        // -----------------------------------------------------------
        CHECK_ME_BBDUK(raw_fastq_ch)

        // -----------------------------------------------------------
        // Step 0b (FASTQ): Quality filter ME-matched reads
        // -----------------------------------------------------------
        input_ch = NANOFILT(CHECK_ME_BBDUK.out.matched)
    }


    // ---------------------------------------------------------------
    // Step 1: Demux
    //   bam input   → dorado emits BAM dir only (--no-emit-fastq)
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
                           .flatMap { bam -> bam in List ? bam : [bam] }
                           .filter { bam -> bam.name.endsWith('.bam') }
                           .map { bam -> tuple(bam.baseName, bam) }
        bam_barcode_ch.view()

        // -----------------------------------------------------------
        // Step 3 (BAM): Convert each per-barcode demuxed BAM → FASTQ
        // -----------------------------------------------------------
        BAM_TO_FASTQ(bam_barcode_ch)

        // -----------------------------------------------------------
        // Step 4 (BAM): Cutadapt trim on the converted FASTQ
        // -----------------------------------------------------------
        CUTADAPT_TRIM(BAM_TO_FASTQ.out.fastq)

        base_trimmed_ch = CUTADAPT_TRIM.out.fq
                          .flatMap { it -> it in List ? it : [it] }
                          .filter { fq -> !fq.name.contains('tagged') && fq.size() > 20 }

        trimmed_ch = base_trimmed_ch
                      .map { fq -> tuple(fq.baseName.replace('.seqWell.fastq',''), fq) }

        trimmed_ch_for_nanostat = base_trimmed_ch
                      .map { fq -> tuple(fq.baseName.replace('.fastq',''), fq) }

        // -----------------------------------------------------------
        // Step 5 (BAM): Filter demux-BAM by cutadapt-FASTQ read IDs
        // -----------------------------------------------------------
        bam_fastq_ch = bam_barcode_ch.join(trimmed_ch, by: 0)

        FILTER_BAM_BY_FASTQ(bam_fastq_ch)

        demuxed_ch = trimmed_ch_for_nanostat

    // ================================================================
    //  FASTQ INPUT BRANCH
    // ================================================================
    } else {

        // -----------------------------------------------------------
        // Step 2 (FASTQ): Extract FASTQ headers BEFORE demuxing
        // -----------------------------------------------------------
        EXTRACT_HEADER(input_ch)
        merged_tags = EXTRACT_HEADER.out.collect()

        // -----------------------------------------------------------
        // Step 3 (FASTQ): Combine per-sample FASTQ demux dirs → per-barcode FASTQs
        // -----------------------------------------------------------
        COMBINE_BARCODES(DORADO_DEMUX.out.fastq_dir.collect())

        fastq_barcode_ch = COMBINE_BARCODES.out.fastq
                             .flatMap { it -> it in List ? it : [it] }
                             .map { fq -> tuple(fq.baseName.replace(".fastq", ""), fq) }

        // -----------------------------------------------------------
        // Step 4 (FASTQ): Restore original header tags into demuxed FASTQs
        // -----------------------------------------------------------
        REHEADER_READS(fastq_barcode_ch, merged_tags)

        // -----------------------------------------------------------
        // Step 5 (FASTQ): Cutadapt trim — on reheadered FASTQ
        // -----------------------------------------------------------
        CUTADAPT_TRIM(REHEADER_READS.out.fq)

        trimmed_ch = CUTADAPT_TRIM.out.fq
                       .flatMap { it ->  it in List ? it : [it] }
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
                    .map { it -> tuple(it[0], it[1]) }

    READ_LENGTH(filtered_ch)
    NANOSTAT(filtered_ch)
    MULTIQC(NANOSTAT.out.collect())
}
