/*
 *  Infer ASVs from processed reads
 */


//// modules to import
include { FILTER_CHIMERA                            } from '../modules/filter_chimera'
include { FILTER_LENGTH                             } from '../modules/filter_length'
include { FILTER_SEQTAB                             } from '../modules/filter_seqtab'

workflow FILTERING {

    take:
    ch_seqtab

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

    // //// filter by PHMM
    // FILTER_PHMM (

    // )

    // //// filter by frameshifts/stop codons
    // FILTER_FRAME (

    // )

    //// filter sequence table (old version)
    FILTER_SEQTAB ( 
        ch_seqtab
    )

    ch_seqtab_filtered = 
        FILTER_SEQTAB.out.seqtab

    emit:

    ch_seqtab_filtered
    ch_read_tracker_grouped = FILTER_SEQTAB.out.read_tracking


}