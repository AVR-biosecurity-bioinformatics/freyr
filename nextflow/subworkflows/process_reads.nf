/*
 *  Process raw sequencing reads, emitting reads ready for ASV inference
 */


//// modules to import
include { SEQ_QC                                    } from '../modules/seq_qc'
include { SPLIT_LOCI                                } from '../modules/split_loci'
include { PRIMER_TRIM                               } from '../modules/primer_trim'
include { READ_FILTER                               } from '../modules/read_filter' 
include { FILTER_QUALPLOTS as FILTER_QUALPLOTS_PRE  } from '../modules/filter_qualplots'
include { FILTER_QUALPLOTS as FILTER_QUALPLOTS_POST } from '../modules/filter_qualplots'

workflow PROCESS_READS {

    take:
    ch_sample_locus_reads
    seq_type 
    pair_type
    miseq_internal
    ch_fcid


    main:

    //// create empty channels
    ch_read_tracker_samples =   // read-tracking for sample-level processes; card: path(.csv)
        Channel.empty()  


    //// run SEQ_QC per flow cell if data is internal MiSeq
    if ( miseq_internal == "true" ) {
        SEQ_QC ( ch_fcid ) 
        }

    //// split sample reads by locus (based on primer seq.)
    SPLIT_LOCI ( 
        ch_sample_locus_reads,
        seq_type,
        pair_type
        ) 

    //// trim primer sequences from the start and end of reads
    PRIMER_TRIM ( 
        SPLIT_LOCI.out.reads 
        )

    //// filter reads using dada2 and input parameters
    READ_FILTER ( 
        PRIMER_TRIM.out.reads 
        )

    //// create plots of read quality pre- and post-filtering, per flowcell (optional)
    FILTER_QUALPLOTS_PRE ( 
        PRIMER_TRIM.out.reads 
        )

    FILTER_QUALPLOTS_POST ( 
        READ_FILTER.out.reads 
        )

    //// concat read_tracker outputs
    ch_read_tracker_samples
        .concat( SPLIT_LOCI.out.input_counts )
        .concat( SPLIT_LOCI.out.read_tracking )
        .concat( PRIMER_TRIM.out.read_tracking )
        .concat( READ_FILTER.out.read_tracking )
        .set { ch_read_tracker_samples }

    ///// split filtered reads into lists of reads per flowcell, primers and direction
    //// forward read channel
    READ_FILTER.out.reads
        .map { meta, reads -> 
                [ "forward", meta.pcr_primers, meta.fcid, reads[0] ] }
        .groupTuple( by: [0,1,2] )
        .set { ch_error_input_fwd }

    //// reverse read channel
    READ_FILTER.out.reads
        .map { meta, reads -> 
                [ "reverse", meta.pcr_primers, meta.fcid, reads[1] ] }
        .groupTuple ( by: [0,1,2] )
        .set { ch_error_input_rev }


    processed_reads = READ_FILTER.out.reads


    emit:

    processed_reads

    /*
    Need to emit:
    - read tracking channel, concat per module (split, primer, filter)
    - channel with inputs for ERROR_MODEL, but combined with first list element "forward", "reverse" or "single"
        depending on `pair_type`; this can be branched later
    */




}