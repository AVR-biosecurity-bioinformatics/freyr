/*
 *  Create result summaries, using phyloseq
 */


//// modules to import
include { ACCUMULATION_CURVE                        } from '../modules/accumulation_curve'
include { PHYLOSEQ_UNFILTERED                       } from '../modules/phyloseq_unfiltered'
include { PHYLOSEQ_FILTER                           } from '../modules/phyloseq_filter'
include { PHYLOSEQ_MERGE                            } from '../modules/phyloseq_merge'
include { READ_TRACKING                             } from '../modules/read_tracking'


workflow RESULT_SUMMARIES {

    take:

    ch_seqtab_filtered  // pcr_primers, seqtab (including sequence), filters per sequence, .fasta combined across flowcells 
    ch_mergetax_output
    ch_loci_samdf
    ch_loci_params
    ch_read_tracker_samples
    ch_read_tracker_grouped


    main:

    // combine taxtables, seqtables and parameters by pcr_primers
    ch_mergetax_output
        .join ( ch_seqtab_filtered, by: 0 )
        .join ( ch_loci_samdf, by: 0 )
        .join ( ch_loci_params, by: 0 )
        .set { ch_phyloseq_input } // pcr_primers, merged_tax, seqtab, filters, fasta, loci_samdf, map[loci_params]

    //// create phyloseq objects per locus; output unfiltered summary tables and accumulation curve plot
    PHYLOSEQ_UNFILTERED ( ch_phyloseq_input )

    //// create accumulation curve plots
    ACCUMULATION_CURVE (
        PHYLOSEQ_UNFILTERED.out.ps
    )

    //// apply taxonomic and minimum abundance filtering per locus (from loci_params), then combine to output filtered summary tables
    PHYLOSEQ_FILTER ( PHYLOSEQ_UNFILTERED.out.ps )

    //// combine phyloseq outputs to merge across loci
    PHYLOSEQ_UNFILTERED.out.ps // val(pcr_primers), path("ps_unfiltered_*.rds"), val(loci_params)
        .map { pcr_primers, ps, filters_tibble, loci_params ->
            [ ps ] }
        .collect()
        .set { ch_ps_unfiltered }

    PHYLOSEQ_FILTER.out.ps // val(pcr_primers), path("ps_filtered_*.rds"), val(loci_params)
        .map { pcr_primers, ps, loci_params ->
            [ ps ] }
        .collect()
        .set { ch_ps_filtered }

    PHYLOSEQ_MERGE ( 
        ch_ps_unfiltered, 
        ch_ps_filtered,
        PHYLOSEQ_UNFILTERED.out.asv_fasta.collect(),
        PHYLOSEQ_FILTER.out.asv_fasta.collect()
        )
    
    ch_read_tracker_grouped = 
        ch_read_tracker_grouped
        .concat( PHYLOSEQ_MERGE.out.read_tracking.flatten() ) 
        .collect()

    // track reads and sequences across the pipeline
    READ_TRACKING ( 
        ch_read_tracker_samples.collect(), 
        ch_read_tracker_grouped 
        )

    ch_read_tracker = READ_TRACKING.out.csv


    emit:

    ch_read_tracker /// TODO: replace with real subworkflow outputs



}