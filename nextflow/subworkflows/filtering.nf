/*
 *  Infer ASVs from processed reads
 */


//// modules to import
include { FILTER_SEQTAB                             } from '../modules/filter_seqtab'

workflow FILTERING {

    take:
    ch_seqtab

    main:

    //// filter sequence table
    FILTER_SEQTAB ( 
        ch_seqtab
    )

    ch_seqtab_filtered = 
        FILTER_SEQTAB.out.seqtab

    emit:

    ch_seqtab_filtered
    ch_read_tracker_grouped = FILTER_SEQTAB.out.read_tracking


}