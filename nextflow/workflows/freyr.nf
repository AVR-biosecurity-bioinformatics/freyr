/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//// import subworkflows
include { PROCESS_READS                             } from '../subworkflows/process_reads/'
include { DADA2                                     } from '../subworkflows/dada2'
include { FILTERING                                 } from '../subworkflows/filtering'
include { TAXONOMY                                  } from '../subworkflows/taxonomy'
include { RESULT_SUMMARIES                          } from '../subworkflows/result_summaries'

//// import modules
include { DOWNSAMPLE_READS                          } from '../modules/downsample_reads/downsample_reads'
include { PARSE_INPUTS                              } from '../modules/parse_inputs/parse_inputs'
include { TRAIN_IDTAXA                              } from '../modules/train_idtaxa/train_idtaxa'

//// utility processes for development and debugging
include { STOP                                      } from '../modules/stop/stop'
workflow FREYR {

    main:

    //// check input option combinations are allowed
    if ( params.seq_type == "illumina" && params.paired ) {
        ch_reads_type = "illumina_paired"
    } else if ( params.seq_type == "nanopore" && !params.paired ) {
        ch_reads_type = params.seq_type
    } else if ( params.illumina && !params.paired ) {
        error "Illumina single-end reads are currently not supported: check --seq_type and --paired" 
    } else if ( params.seq_type == "nanopore" && params.paired ) {
        error "Nanopore reads cannot be paired-end: check --seq_type and --paired" 
    } else if ( params.seq_type == "pacbio" ) {
        error "PacBio reads are currently not supported: check --seq_type" 
    } else { 
        error "Disallowed combination of --seq_type and --paired"
    }

    //// parse path channels
    ch_samplesheet_file = channel.fromPath( params.samplesheet, checkIfExists: true, type: 'file' )

    if ( params.primer_params ){
        ch_primer_params_file = channel.fromPath( params.primer_params, checkIfExists: true, type: 'file' )
        ch_primer_params_type = "user"
    } else {
        ch_primer_params_file = channel.fromPath( "assets/primer_params_default.csv", checkIfExists: true, type: 'file' )
        ch_primer_params_type = "default"
    }

    //// parse subsample parameter
    ch_subsample = params.subsample ?: "false"

    //// parse extension parameter 
    ch_extension = params.extension ?: "null"

    //// Create empty channels
    ch_read_tracker_grouped =   // read-tracking for grouped processes
        channel.empty()      

    //// create pp_* parameter map as input for PARSE_INPUTS
    channel.from(params)
        .map { params_map ->
            def pp_params = 
                params_map.subMap(
                    'pp_primers', 'pp_locus', 'pp_for_primer_seq', 'pp_rev_primer_seq', 
                    'pp_max_primer_mismatch', 'pp_read_min_length', 'pp_read_max_length', 'pp_read_max_ee', 
                    'pp_read_trunc_length', 'pp_read_trim_left', 'pp_read_trim_right', 'pp_asv_min_length', 
                    'pp_asv_max_length', 'pp_concat_unmerged', 'pp_genetic_code', 'pp_coding', 
                    'pp_phmm', 'pp_idtaxa_db', 'pp_ref_fasta', 'pp_idtaxa_confidence', 
                    'pp_run_blast', 'pp_blast_min_identity', 'pp_blast_min_coverage', 'pp_cluster_threshold', 
                    'pp_target_kingdom', 'pp_target_phylum', 'pp_target_class', 'pp_target_order', 
                    'pp_target_family', 'pp_target_genus', 'pp_target_species', 'pp_min_sample_reads', 
                    'pp_min_taxa_reads', 'pp_min_taxa_ra'
                )
            pp_params }
        .set { ch_pp_params } 

    //// read-in samplesheet and primer_params .csv files, validate their contents, and produce inputs for rest of pipeline
    PARSE_INPUTS ( 
        ch_samplesheet_file, 
        ch_primer_params_file,
        ch_primer_params_type,
        params.seq_type,
        params.paired,
        ch_subsample,
        ch_extension,
        ch_pp_params,
        workflow.launchDir.toString()
    )

    ch_samplesheet_split = PARSE_INPUTS.out.samplesheet_split.first()

    ch_sample_metadata = PARSE_INPUTS.out.sample_metadata.first()

    //// create channel of just read files for PROCESS_READS
    if ( params.paired ){
        PARSE_INPUTS.out.samplesheet_unsplit
            .flatten ()
            .splitCsv ( header: true )
            .map { header -> 
                [ header.sample, [ file(header.fwd, checkIfExists: true), file(header.rev, checkIfExists: true) ] ] 
                }
            .set { ch_sample_reads }
    } else if ( !params.paired ) {
        PARSE_INPUTS.out.samplesheet_unsplit
            .flatten ()
            .splitCsv ( header: true )
            .map { header -> 
                [ header.sample, file(header.single, checkIfExists: true) ] 
                }
            .set { ch_sample_reads }
    } else {
        error " '--paired' must be 'true' or 'false'. "
    }
    
    //// create channel of read files with sample, read_group, primers and sample_primers
    if ( params.paired ){
        ch_samplesheet_split
            .splitCsv ( header: true )
            .map { header -> 
                [ header.primers, header.read_group, header.sample, header.sample_primers, [ file(header.fwd, checkIfExists: true), file(header.rev, checkIfExists: true) ] ] 
                }
            .set { ch_sample_primers_reads }
    } else if ( !params.paired ) {
        ch_samplesheet_split
            .splitCsv ( header: true )
            .map { header -> 
                [ header.primers, header.read_group, header.sample, header.sample_primers, file(header.single, checkIfExists: true) ] 
                }
            .set { ch_sample_primers_reads }
    } else {
        error " '--paired' must be 'true' or 'false'. "
    }

    //// create channel of primer parameters
    PARSE_INPUTS.out.primer_params_parsed
        .flatten()
        .splitCsv ( header: true )
        .map { header -> 
            [ header.primers, header ]
            }
        .set { ch_primer_params }
    
    //// train IDTAXA model from reference database .fasta
    if ( params.train_idtaxa ) {
        
        //// create input channel for TRAIN_IDTAXA
        ch_primer_params
            .map { primers, primer_params ->
                [ primers, primer_params.ref_fasta ]  }
            .set { ch_train_idtaxa_input }

        //// train model
        TRAIN_IDTAXA (
            ch_train_idtaxa_input,
            params.train_idtaxa
        )

        //// update primer_params with new models
        ch_primer_params
            .join ( TRAIN_IDTAXA.out.model, by: 0 ) // join by primers
            .map { primers, primer_params, new_idtaxa ->
                [ primers, primer_params + [ idtaxa_db: new_idtaxa ] ] }
            .set { ch_primer_params }
    
    } 

    //// downsample reads if params.downsample is defined
    if ( params.downsample_reads ) {

        DOWNSAMPLE_READS (
            ch_sample_primers_reads,
            params.seq_type,
            params.paired,
            params.downsample_reads
        )

        ch_process_reads_input = DOWNSAMPLE_READS.out.reads

    } else {
        ch_process_reads_input = ch_sample_primers_reads
    }

    //// subworkflow: process sequencing reads
    PROCESS_READS (
        ch_sample_reads,
        ch_process_reads_input,
        ch_primer_params,
        params.seq_type,
        params.paired
    )


    //// subworkflow: run DADA2 to infer ASVs
    DADA2 (
        PROCESS_READS.out.ch_processed_reads,
        ch_primer_params
    )

    //// subworkflow: apply soft filters to ASVs
    FILTERING (
        DADA2.out.ch_seqtab,
        ch_primer_params,
        ch_samplesheet_split
    )

    //// subworkflow: assign taxonomy
    TAXONOMY (
        DADA2.out.ch_seqtab,
        ch_primer_params
    )

    //// combine read tracking grouped channels
    ch_read_tracker_grouped = DADA2.out.ch_read_tracker_grouped
        .concat(FILTERING.out.ch_read_tracker_grouped)

   
    //// subworkflow: create result summaries
    RESULT_SUMMARIES (
        FILTERING.out.ch_seqtab_filtered,
        TAXONOMY.out.ch_mergetax_output,
        ch_samplesheet_split,
        ch_sample_metadata,
        ch_primer_params,
        PROCESS_READS.out.ch_read_tracker_samples,
        ch_read_tracker_grouped
    )
    

    ///// VISUALISATION
    /*
    - stacked bar charts of ASV abundance across samples -- filtered vs unfiltered
    - plots showing proportion of sequences that passed or failed the final filters, per sample/flowcell/locus etc
    - taxonomy heatmap tree
    */

    // TODO: adopt workflow outputs

}