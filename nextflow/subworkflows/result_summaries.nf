/*
 *  Create result summaries, using phyloseq
 */


//// modules to import
include { ACCUMULATION_CURVE                        } from '../modules/accumulation_curve'
include { PHYLOSEQ_CLUSTERED                        } from '../modules/phyloseq_clustered'
include { PHYLOSEQ_FILTER                           } from '../modules/phyloseq_filter'
include { PHYLOSEQ_MERGE                            } from '../modules/phyloseq_merge'
include { PHYLOSEQ_UNFILTERED                       } from '../modules/phyloseq_unfiltered'
include { READ_TRACKING                             } from '../modules/read_tracking'


workflow RESULT_SUMMARIES {

    take:

    ch_seqtab_filtered  // primers, seqtab (including sequence), filters per sequence, .fasta combined across flowcells 
    ch_mergetax_output
    ch_samplesheet_split
    ch_sample_metadata
    ch_primer_params
    ch_read_tracker_samples
    ch_read_tracker_grouped


    main:

    //// combine PHYLOSEQ_UNFILTERED parameters to input channel
    ch_primer_params
        .map { primers, primer_params ->
            def process_params = 
                primer_params.subMap('cluster_threshold')
            [ primers, process_params ] }
        .set { ch_phyloseq_unfiltered_params } 

    ch_mergetax_output
        .join ( ch_seqtab_filtered, by: 0 )
        .combine ( ch_samplesheet_split )
        .combine ( ch_sample_metadata )
        .combine ( ch_phyloseq_unfiltered_params, by: 0 )
        .set { ch_phyloseq_unfiltered_input } // primers, merged_tax, seqtab, filters, fasta, samplesheet_split, sample_metadata, process_params

    //// create phyloseq objects per locus; output unfiltered summary tables and accumulation curve plot
    PHYLOSEQ_UNFILTERED ( 
        ch_phyloseq_unfiltered_input 
    )

    //// combine ACCUMULATION_CURVE parameters to input channel
    ch_primer_params
        .map { primers, primer_params ->
            def process_params = 
                primer_params.subMap('min_sample_reads')
            [ primers, process_params ] }
        .set { ch_accumulation_curve_params } 

    PHYLOSEQ_UNFILTERED.out.ps
        .combine ( ch_accumulation_curve_params, by: 0 )
        .set { ch_accumulation_curve_input }

    //// create accumulation curve plots
    ACCUMULATION_CURVE (
        ch_accumulation_curve_input
    )


    //// combine PHYLOSEQ_FILTER parameters to input channel
    ch_primer_params
        .map { primers, primer_params ->
            def process_params = 
                primer_params.subMap('target_kingdom','target_phylum','target_class','target_order','target_family','target_genus','target_species','min_sample_reads','min_taxa_reads','min_taxa_ra')
            [ primers, process_params ] }
        .set { ch_phyloseq_filter_params } 

    PHYLOSEQ_UNFILTERED.out.ps
        .combine ( ch_phyloseq_filter_params, by: 0 )
        .set { ch_phyloseq_filter_input }

    //// apply taxonomic and minimum abundance filtering per locus (from loci_params), then combine to output filtered summary tables
    PHYLOSEQ_FILTER ( 
        ch_phyloseq_filter_input
    )


    //// merge each filtered ASV cluster into a single representative sequence
    // PHYLOSEQ_CLUSTERED (

    // )

    //// combine phyloseq outputs to merge across loci
    PHYLOSEQ_UNFILTERED.out.ps // val(primers), path("ps_unfiltered_*.rds")
        .map { primers, ps, filters_tibble -> [ ps ] }
        .collect()
        .set { ch_ps_unfiltered }

    PHYLOSEQ_FILTER.out.ps // val(primers), path("ps_filtered_*.rds")
        .map { primers, ps -> [ ps ] }
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
        ch_read_tracker_grouped,
        ch_samplesheet_split
    )

    ch_read_tracker = READ_TRACKING.out.csv


    emit:

    ch_read_tracker /// TODO: replace with real subworkflow outputs



}