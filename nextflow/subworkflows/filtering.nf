/*
 *  Infer ASVs from processed reads
 */


//// modules to import
include { FILTER_CHIMERA                            } from '../modules/filter_chimera'
include { FILTER_FRAME                              } from '../modules/filter_frame'
include { FILTER_LENGTH                             } from '../modules/filter_length'
include { FILTER_PHMM                               } from '../modules/filter_phmm'
include { MERGE_FILTERS                             } from '../modules/merge_filters'

workflow FILTERING {

    take:
    ch_seqtab
    ch_read_groups
    ch_primer_params
    ch_samplesheet_split

    main:

    //// group all flowcells (per pcr primer pair) into a single element for filtering
    ch_seqtab
        .map { pcr_primers, read_groups, seqtab_tibble, fasta -> 
            [ pcr_primers, seqtab_tibble, fasta ] }
        .groupTuple ( by: 0 )
        .set { ch_filter_input }


    //// combine FILTER_CHIMERA parameters to input channel
    ch_primer_params
        .map { primers, primer_params ->
            [ primers, [ 'chimera_sample_frac': params.chimera_sample_frac ] ] }
        .set { ch_filter_chimera_params } 

    ch_filter_input
        .combine ( ch_filter_chimera_params, by: 0 )
        .set { ch_filter_chimera_input }

    //// filter chimeras
    FILTER_CHIMERA (
        ch_filter_chimera_input
    )


    //// combine FILTER_LENGTH parameters to input channel
    ch_primer_params
        .map { primers, primer_params ->
            def process_params = 
            primer_params.subMap('asv_min_length','asv_max_length')
            [ primers, process_params ] }
        .set { ch_filter_length_params } 

    ch_filter_input
        .combine ( ch_filter_length_params, by: 0 )
        .set { ch_filter_length_input }

    //// filter by length
    FILTER_LENGTH (
        ch_filter_length_input
    )


    //// combine FILTER_PHMM parameters to input channel
    ch_primer_params
        .map { primers, primer_params ->
            def process_params = 
            primer_params.subMap('for_primer_seq','rev_primer_seq', 'coding', 'phmm')
            [ primers, process_params ] }
        .set { ch_filter_phmm_params } 

    ch_filter_input
        .combine ( ch_filter_phmm_params, by: 0 )
        .set { ch_filter_phmm_input }

    //// filter by PHMM
    FILTER_PHMM (
        ch_filter_phmm_input
    )


    //// combine FILTER_FRAME parameters to input channel
    ch_primer_params
        .map { primers, primer_params ->
            def process_params = 
            primer_params.subMap('coding', 'genetic_code')
            [ primers, process_params ] }
        .set { ch_filter_frame_params } 

    ch_filter_input
        .combine ( ch_filter_frame_params, by: 0 )
        .set { ch_filter_frame_input }

    //// filter by frameshifts/stop codons
    FILTER_FRAME (
        ch_filter_frame_input
    )


    //// group filtering outputs (filters, setabs and fastas) by pcr_primers
    FILTER_CHIMERA.out.tibble
        .mix ( FILTER_LENGTH.out.tibble )
        .mix ( FILTER_PHMM.out.tibble )
        .mix ( FILTER_FRAME.out.tibble )
        .groupTuple ( by: 0 )
        .join ( ch_filter_input, by: 0 )
        .set { ch_merge_filters_input }
    
    //// merge filters together, create filter read tracking tibble, and create filter plots
    MERGE_FILTERS (
        ch_merge_filters_input,
        ch_samplesheet_split
    )


    emit:

    ch_seqtab_filtered = MERGE_FILTERS.out.filtered
    ch_read_tracker_grouped = MERGE_FILTERS.out.read_tracking


}