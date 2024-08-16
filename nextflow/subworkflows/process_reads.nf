/*
 *  Process raw sequencing reads, emitting reads ready for ASV inference
 */


//// modules to import
include { MISEQ_QC                                  } from '../modules/miseq_qc'
include { SPLIT_LOCI                                } from '../modules/split_loci'
include { PRIMER_TRIM                               } from '../modules/primer_trim'
include { READ_FILTER                               } from '../modules/read_filter' 
include { FILTER_QUALPLOTS as FILTER_QUALPLOTS_PRE  } from '../modules/filter_qualplots'
include { FILTER_QUALPLOTS as FILTER_QUALPLOTS_POST } from '../modules/filter_qualplots'

workflow PROCESS_READS {

    take:
    ch_sample_locus_reads
    seq_type 
    paired
    ch_fcid


    main:

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
        MISEQ_QC ( ch_fcid ) 
        }

    //// split sample reads by locus (based on primer seq.)
    SPLIT_LOCI ( 
        ch_sample_locus_reads,
        seq_type,
        paired
        ) 

    //// trim primer sequences from the start and end of reads
    PRIMER_TRIM ( 
        SPLIT_LOCI.out.reads,
        seq_type,
        paired 
        )

    //// filter reads using dada2 and input parameters
    READ_FILTER ( 
        PRIMER_TRIM.out.reads,
        seq_type,
        paired
        )

    //// create plots of read quality pre- and post-filtering, per flowcell (optional)
    FILTER_QUALPLOTS_PRE ( 
        PRIMER_TRIM.out.reads,
        seq_type,
        paired,
        "pre"
        )

    FILTER_QUALPLOTS_POST ( 
        READ_FILTER.out.reads,
        seq_type,
        paired,
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