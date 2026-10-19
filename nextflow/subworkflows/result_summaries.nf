/*
 *  Create result summaries, using phyloseq
 */


//// modules to import
include { ACCUMULATION_CURVE                        } from '../modules/accumulation_curve/accumulation_curve'
include { PHYLOSEQ_CLUSTERED                        } from '../modules/phyloseq_clustered/phyloseq_clustered'
include { PHYLOSEQ_FILTER                           } from '../modules/phyloseq_filter/phyloseq_filter'
include { PHYLOSEQ_UNFILTERED                       } from '../modules/phyloseq_unfiltered/phyloseq_unfiltered'
include { PHYLOSEQ_MERGE as MERGE_UNFILTERED        } from '../modules/phyloseq_merge/phyloseq_merge'
include { PHYLOSEQ_MERGE as MERGE_FILTERED          } from '../modules/phyloseq_merge/phyloseq_merge'
include { READ_TRACKING                             } from '../modules/read_tracking/read_tracking'


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

    PHYLOSEQ_UNFILTERED.out.components
        .map { primers, seqtab, _taxtab, samdf, _fasta, _filters -> [ primers, seqtab, samdf ]}
        .combine(ch_accumulation_curve_params, by: 0)
        .set { ch_accumulation_curve_input }

    //// create accumulation curve plots
    ACCUMULATION_CURVE (
        ch_accumulation_curve_input
    )


    //// combine PHYLOSEQ_FILTER parameters to input channel
    ch_primer_params
        .map { primers, primer_params ->
            def process_params = 
                primer_params.subMap('cluster_threshold','target_kingdom','target_phylum','target_class','target_order','target_family','target_genus','target_species','min_sample_reads','min_taxa_reads','min_taxa_ra')
            [ primers, process_params ] }
        .set { ch_phyloseq_filter_params } 

    PHYLOSEQ_UNFILTERED.out.components
        .combine(ch_phyloseq_filter_params, by: 0)
        .set { ch_phyloseq_filter_input }


    //// apply taxonomic and minimum abundance filtering per locus (from loci_params), then combine to output filtered summary tables
    PHYLOSEQ_FILTER ( 
        ch_phyloseq_filter_input
    )

    if ( params.merge_clusters ){
        
        PHYLOSEQ_FILTER.out.csvs
            .combine ( PHYLOSEQ_FILTER.out.ps, by: 0 )
            .combine ( PHYLOSEQ_FILTER.out.clusters, by: 0)
            .set { ch_phyloseq_clustered_input } // primers, seqtab, taxtab, samdf, raw, summary, ps, clusters

        //// merge each filtered ASV cluster into a single representative sequence
        PHYLOSEQ_CLUSTERED (
            ch_phyloseq_clustered_input,
            params.merge_clusters
        )
        
    }

    //// Collect unfiltered components across loci.
    PHYLOSEQ_UNFILTERED.out.components
        .map { primers, seqtab, taxtab, samdf, fasta, _filters ->
            tuple(seqtab, taxtab, samdf, fasta)
        }
        .collect(flat: false)
        .map { components ->
            tuple(
                "unfiltered",
                components.collect { component -> component[0] },
                components.collect { component -> component[1] },
                components.collect { component -> component[2] },
                components.collect { component -> component[3] }
            )
        }
        .set { ch_merge_unfiltered_input }

    //// Collect filtered components across loci.
    PHYLOSEQ_FILTER.out.components
        .map { primers, seqtab, taxtab, samdf, fasta, _filters ->
            tuple(seqtab, taxtab, samdf, fasta)
        }
        .collect(flat: false)
        .map { components ->
            tuple(
                "filtered",
                components.collect { component -> component[0] },
                components.collect { component -> component[1] },
                components.collect { component -> component[2] },
                components.collect { component -> component[3] }
            )
        }
        .set { ch_merge_filtered_input }


    //// Merge filtered and unfiltered objects in separate process tasks.
    MERGE_UNFILTERED(
        ch_merge_unfiltered_input
    )

    MERGE_FILTERED(
        ch_merge_filtered_input
    )

    //// Add merged filtered and unfiltered read-tracking outputs.
    ch_read_tracker_grouped_final = ch_read_tracker_grouped
        .concat(
            MERGE_UNFILTERED.out.read_tracking,
            MERGE_FILTERED.out.read_tracking
        )
        .collect()

    //// Track reads and sequences across the pipeline.
    READ_TRACKING(
        ch_read_tracker_samples.collectFile(
            name: "read_tracker_samples.csv"
        ),
        ch_read_tracker_grouped_final,
        ch_samplesheet_split
    )

    ch_read_tracker = READ_TRACKING.out.csv


    emit:

    ch_read_tracker /// TODO: replace with real subworkflow outputs



}