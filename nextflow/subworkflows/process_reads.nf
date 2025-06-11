/*
 *  Process raw sequencing reads, emitting reads ready for ASV inference
 */


//// modules to import
include { FASTQC                                    } from '../modules/fastqc'
include { NANOPLOT                                  } from '../modules/nanoplot'
include { MISEQ_QC                                  } from '../modules/miseq_qc'
include { SPLIT_LOCI                                } from '../modules/split_loci'
include { PRIMER_TRIM                               } from '../modules/primer_trim'
include { READ_FILTER                               } from '../modules/read_filter' 
include { FILTER_QUALPLOTS as FILTER_QUALPLOTS_PRE  } from '../modules/filter_qualplots'
include { FILTER_QUALPLOTS as FILTER_QUALPLOTS_POST } from '../modules/filter_qualplots'

workflow PROCESS_READS {

    take:
    ch_sample_reads
    ch_sample_primers_reads
    ch_primer_params
    seq_type 
    paired
    ch_read_groups


    main:

    ch_seq_type = channel.value ( seq_type )
    ch_paired = channel.value ( paired )

    if ( params.miseq_dir ) {
        ch_miseq_dir = Channel.fromPath(params.miseq_dir).first()
    } else {
        ch_miseq_dir = Channel.empty()
    }

    //// create empty channels
    ch_read_tracker_samples =   // read-tracking for sample-level processes; card: path(.csv)
        Channel.empty()  
    ch_processed_reads = 
        Channel.empty()
    ch_processed_fwd =
        Channel.empty()
    ch_processed_rev =
        Channel.empty()
    ch_processed_single =
        Channel.empty()


    //// run SEQ_QC per flow cell if data is internal MiSeq
    if ( params.miseq_internal == true ) {

        if ( !ch_miseq_dir ) {
            error "'--miseq_dir' must be specified when using '--miseq_internal'"
        }

        MISEQ_QC ( 
            ch_read_groups,
            ch_miseq_dir
        ) 
    }

    //// run fastqc on input reads
    FASTQC (
        ch_sample_reads,
        ch_seq_type,
        ch_paired
        )
        
    //// run nanoplot on nanopore reads
    if ( seq_type == "nanopore" ) {
        NANOPLOT (
            ch_sample_reads,
            ch_seq_type,
            ch_paired
            )
    }

    //// combine map of SPLIT_LOCI primer parameters to read channel
    ch_primer_params
        .map { primers, primer_params ->
            def process_params = 
                primer_params.subMap('for_primer_seq','rev_primer_seq','locus') + 
                [ 'seq_type': params.seq_type, 'paired': params.paired, 'primer_error_rate': params.primer_error_rate ]
            [ primers, process_params ] }
        .set { ch_split_loci_params }
    
    ch_sample_primers_reads
        .combine ( ch_split_loci_params, by: 0 )
        .set { ch_split_loci_input }

    //// split sample reads by primers (based on primer seq.)
    SPLIT_LOCI ( 
        ch_split_loci_input
    ) 

    //// combine map of PRIMER_TRIM primer parameters to read channel
    ch_primer_params
        .map { primers, primer_params ->
            def process_params = 
                primer_params.subMap('for_primer_seq','rev_primer_seq','locus') + 
                [ 'seq_type': params.seq_type, 'paired': params.paired, 'primer_error_rate': params.primer_error_rate, 'primer_n_trim': params.primer_n_trim ]
            [ primers, process_params ] }
        .set { ch_primer_trim_params }
    
    SPLIT_LOCI.out.reads
        .combine ( ch_primer_trim_params, by: 0 )
        .set { ch_primer_trim_input }

    //// trim primer sequences from the start and end of reads
    PRIMER_TRIM ( 
        ch_primer_trim_input 
    )

    //// filter reads using dada2 and input parameters
    READ_FILTER ( 
        PRIMER_TRIM.out.reads,
        ch_seq_type,
        ch_paired
        )

    //// create plots of read quality pre- and post-filtering, per flowcell (optional)
    FILTER_QUALPLOTS_PRE ( 
        PRIMER_TRIM.out.reads,
        ch_seq_type,
        ch_paired,
        "pre"
        )

    FILTER_QUALPLOTS_POST ( 
        READ_FILTER.out.reads,
        ch_seq_type,
        ch_paired,
        "post"
        )

    //// concat read_tracker outputs
    ch_read_tracker_samples
        .concat( SPLIT_LOCI.out.input_counts )
        .concat( SPLIT_LOCI.out.read_tracking )
        .concat( PRIMER_TRIM.out.read_tracking )
        .concat( READ_FILTER.out.read_tracking )
        .set { ch_read_tracker_samples }

    //// split filtered reads into channels per flowcell, primers and direction
    if ( params.paired == true ) {
        //// forward read channel
        READ_FILTER.out.reads
            .map { meta, reads -> 
                    [ "forward", meta.pcr_primers, meta.fcid, meta, reads[0] ] }
            // .groupTuple( by: [0,1,2] )
            .set { ch_processed_fwd }

        //// reverse read channel
        READ_FILTER.out.reads
            .map { meta, reads -> 
                    [ "reverse", meta.pcr_primers, meta.fcid, meta, reads[1] ] }
            // .groupTuple ( by: [0,1,2] )
            .set { ch_processed_rev }

    } else if ( params.paired == false ) {
        //// single-end read channel
        READ_FILTER.out.reads
            .map { meta, reads -> 
                    [ "single", meta.pcr_primers, meta.fcid, meta, reads ] }
            // .groupTuple ( by: [0,1,2] )
            .set { ch_processed_single }
    
    } else {
        error ("Disallowed 'params.paired' value")
    }

    ch_processed_reads
        .concat (ch_processed_fwd)
        .concat (ch_processed_rev)
        .concat (ch_processed_single)
        .set { ch_processed_reads }


    emit:

    ch_processed_reads
    ch_read_tracker_samples


}