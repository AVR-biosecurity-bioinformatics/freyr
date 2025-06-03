/*
 *  Infer ASVs from processed reads
 */


//// modules to import
include { FILTER_CHIMERA                            } from '../modules/filter_chimera'
include { FILTER_SEQTAB                             } from '../modules/filter_seqtab'

workflow FILTERING {

    take:
    ch_seqtab

    main:

    //// group all flowcells (per pcr primer pair) into a single element for chimera filtering
    ch_seqtab
        .map { pcr_primers, fcid, meta, seqtab_tibble, fasta -> 
            [ pcr_primers, seqtab_tibble, fasta ] }
        .groupTuple ( by: 0 )
        .set { ch_chimera_input }

    //// filter chimeras
    FILTER_CHIMERA (
        ch_chimera_input,
        params.chimera_sample_frac
    )



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