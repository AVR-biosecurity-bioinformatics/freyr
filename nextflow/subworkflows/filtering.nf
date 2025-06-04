/*
 *  Infer ASVs from processed reads
 */


//// modules to import
include { FILTER_CHIMERA                            } from '../modules/filter_chimera'
include { FILTER_FRAME                              } from '../modules/filter_frame'
include { FILTER_LENGTH                             } from '../modules/filter_length'
include { FILTER_PHMM                               } from '../modules/filter_phmm'
include { FILTER_SEQTAB                             } from '../modules/filter_seqtab'
include { MERGE_FILTERS                             } from '../modules/merge_filters'

workflow FILTERING {

    take:
    ch_seqtab
    ch_fcid

    main:

    //// group all flowcells (per pcr primer pair) into a single element for filtering
    ch_seqtab
        .map { pcr_primers, fcid, meta, seqtab_tibble, fasta -> 
            [ pcr_primers, meta, seqtab_tibble, fasta ] }
        .groupTuple ( by: 0 )
        .set { ch_filter_input }

    //// filter chimeras
    FILTER_CHIMERA (
        ch_filter_input,
        params.chimera_sample_frac
    )

    //// filter by length
    FILTER_LENGTH (
        ch_filter_input
    )

    //// filter by PHMM
    FILTER_PHMM (
        ch_filter_input
    )

    //// filter by frameshifts/stop codons
    FILTER_FRAME (
        ch_filter_input
    )

    //// modify filter input channel to just keep seqtabs and fasta
    ch_filter_input
        .map { pcr_primers, meta, seqtab_tibble_list, fasta_list -> 
            [ pcr_primers, seqtab_tibble_list, fasta_list ] }
        .set { ch_grouped_seqtabs }

    //// group filtering outputs (filters, setabs and fastas) by pcr_primers
    FILTER_CHIMERA.out.tibble
        .mix ( FILTER_LENGTH.out.tibble )
        .mix ( FILTER_PHMM.out.tibble )
        .mix ( FILTER_FRAME.out.tibble )
        .groupTuple ( by: 0 )
        .join ( ch_grouped_seqtabs, by: 0 )
        .set { ch_merge_filters_input }
    
    //// merge filters together, createfilter  read tracking tibble, and create filter plots
    MERGE_FILTERS (
        ch_merge_filters_input,
        ch_fcid.collect()
    )

    // //// filter sequence table (old version)
    // FILTER_SEQTAB ( 
    //     ch_seqtab
    // )


    emit:

    ch_seqtab_filtered = MERGE_FILTERS.out.filtered
    ch_read_tracker_grouped = MERGE_FILTERS.out.read_tracking


}