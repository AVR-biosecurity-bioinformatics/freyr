#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AVR-biosecurity-bioinformatics/freyr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/AVR-biosecurity-bioinformatics/freyr
----------------------------------------------------------------------------------------
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include functions from nf-schema
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema' 

// 
// from https://github.com/nextflow-io/nextflow/issues/1129

def startupMessage() {
    log.info pipelineHeader()
    log.info "~~~ freyr: Metabarcoding analysis for biosecurity and biosurveillance ~~~"
    log.info " "
}

def pipelineHeader(){
    return """
    
    :::::::::: :::::::::  :::::::::: :::   ::: :::::::::  
    :+:        :+:    :+: :+:        :+:   :+: :+:    :+: 
    +:+        +:+    +:+ +:+         +:+ +:+  +:+    +:+ 
    :#::+::#   +#++:++#:  +#++:++#     +#++:   +#++:++#:  
    +#+        +#+    +#+ +#+           +#+    +#+    +#+ 
    #+#        #+#    #+# #+#           #+#    #+#    #+# 
    ###        ###    ### ##########    ###    ###    ### 
    
    """.stripIndent()
}

startupMessage()

workflow.onComplete {
    if ( workflow.success ) {
      log.info "[$workflow.complete] >> Pipeline finished SUCCESSFULLY after $workflow.duration"
    } else {
      log.info "[$workflow.complete] >> Pipeline finished with ERRORS after $workflow.duration"
    }
    /*
    TODO: use other metadata (https://www.nextflow.io/docs/latest/metadata.html) to display at the end of a run.
    */

}


// Print help message, supply typical command line usage for the pipeline
if (params.help) {
//    log.info startupMessage()
   log.info paramsHelp("nextflow run AVR-biosecurity-bioinformatics/freyr --samplesheet samplesheet.csv --loci_params loci_params.csv") // TODO: add typical commands for pipeline
   exit 0
}

// Validate input parameters using schema
validateParameters( parameters_schema: 'nextflow_schema.json' )


// Print summary of supplied parameters (that differ from defaults)
log.info paramsSummaryLog(workflow)
// make it clear that samples are being subsampled
if (params.subsample) {
    log.info "*** NOTE: Input samples are being randomly subsampled to $params.subsample per primer x flowcell combination (params.subsample = $params.subsample) ***"
}

if( !nextflow.version.matches('=23.04.5') ) {
    println " "
    println "*** ERROR ~ This pipeline currently requires Nextflow version 23.04.5 -- You are running version ${nextflow.version}. ***"
    error "*** You can use version 23.04.5 by appending 'NXF_VER=23.04.5' to the front of the 'nextflow run' command. ***"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
// ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
// ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
// ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

// TODO: What channels do I need for config? This could be used to define custom visualisation or output summary formats
// TODO: Set maximum number of cores/cpus at the total number of read pairs/samples

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT AND VARIABLES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//// Inputs





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//// import subworkflows
include { PROCESS_READS                             } from './nextflow/subworkflows/process_reads'
include { DADA2                                     } from './nextflow/subworkflows/dada2'
include { TAXONOMY                                  } from './nextflow/subworkflows/taxonomy'
include { RESULT_SUMMARIES                          } from './nextflow/subworkflows/result_summaries'

//// import modules
include { PARSE_INPUTS                              } from './nextflow/modules/parse_inputs'
include { MISEQ_QC                                  } from './nextflow/modules/miseq_qc'
include { SPLIT_LOCI                                } from './nextflow/modules/split_loci'
include { PRIMER_TRIM                               } from './nextflow/modules/primer_trim'
include { READ_FILTER                               } from './nextflow/modules/read_filter' 
include { FILTER_QUALPLOTS as FILTER_QUALPLOTS_PRE  } from './nextflow/modules/filter_qualplots'
include { FILTER_QUALPLOTS as FILTER_QUALPLOTS_POST } from './nextflow/modules/filter_qualplots'
include { ERROR_MODEL as ERROR_MODEL_F              } from './nextflow/modules/error_model'
include { ERROR_MODEL as ERROR_MODEL_R              } from './nextflow/modules/error_model'
include { DENOISE as DENOISE1_F                     } from './nextflow/modules/denoise'
include { DENOISE as DENOISE1_R                     } from './nextflow/modules/denoise'
include { PRIORS as PRIORS_F                        } from './nextflow/modules/priors'
include { PRIORS as PRIORS_R                        } from './nextflow/modules/priors'
include { DENOISE as DENOISE2_F                     } from './nextflow/modules/denoise'
include { DENOISE as DENOISE2_R                     } from './nextflow/modules/denoise'
include { MAKE_SEQTAB_PAIRED                        } from './nextflow/modules/make_seqtab_paired'
include { FILTER_SEQTAB                             } from './nextflow/modules/filter_seqtab'
include { TAX_IDTAXA                                } from './nextflow/modules/tax_idtaxa'
include { TAX_BLAST                                 } from './nextflow/modules/tax_blast'
include { JOINT_TAX                                 } from './nextflow/modules/joint_tax'
include { MERGE_TAX                                 } from './nextflow/modules/merge_tax'
include { ASSIGNMENT_PLOT                           } from './nextflow/modules/assignment_plot'
include { TAX_SUMMARY                               } from './nextflow/modules/tax_summary'
include { TAX_SUMMARY_MERGE                         } from './nextflow/modules/tax_summary_merge'
include { PHYLOSEQ_UNFILTERED                       } from './nextflow/modules/phyloseq_unfiltered'
include { PHYLOSEQ_FILTER                           } from './nextflow/modules/phyloseq_filter'
include { PHYLOSEQ_MERGE                            } from './nextflow/modules/phyloseq_merge'
include { READ_TRACKING                             } from './nextflow/modules/read_tracking'

// utility processes for development and debugging
include { STOP                                      } from './nextflow/modules/stop'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

// include { FASTQC                            } from '../modules/nf-core/fastqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE AND RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FREYR {


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

    //// read-in samplesheet and loci_params .csv files, validate their contents, and produce inputs for rest of pipeline
    PARSE_INPUTS ( params.samplesheet, params.loci_params )

    // ch_versions = Channel.empty()
    
    //// Create empty channels
    ch_read_tracker_grouped =   // read-tracking for grouped processes
        Channel.empty()      
    

    //// parse samplesheets that contain locus-specific parameters
    PARSE_INPUTS.out.samplesheet_locus
        .flatten ()
        .splitCsv ( header: true )
        .map { row -> 
            def meta = row.subMap(
                'sample_id','sample_name','extraction_rep','amp_rep',
                'client_name','experiment_name','sample_type','collection_method',
                'collection_location','latitude','longitude','environment','collection_date',
                'operator_name','description','assay','extraction_method',
                'amp_method','target_gene','pcr_primers','for_primer_seq',
                'rev_primer_seq','index_plate','index_well','i7_index_id',
                'i7_index','i5_index_id','i5_index','seq_platform',
                'fcid','for_read_length','rev_read_length','seq_run_id',
                'seq_id','seq_date','analysis_method','notes','max_primer_mismatch','read_min_length','read_max_length',
                'read_max_ee','read_trunc_length','read_trim_left','read_trim_right',
                'asv_min_length','asv_max_length','high_sensitivity','concat_unmerged','genetic_code','coding',
                'phmm','idtaxa_db','ref_fasta','idtaxa_confidence',
                'run_blast','blast_min_identity','blast_min_coverage','target_kingdom',
                'target_phylum','target_class','target_order','target_family',
                'target_genus','target_species','min_sample_reads','min_taxa_reads',
                'min_taxa_ra','threads'
                )
            [ meta, [
                file(row.fwd, checkIfExists: true),
                file(row.rev, checkIfExists: true)
            ] ]  
            }
        .set { ch_sample_locus_reads }


    //// create channel that links locus-specific samplesheets to pcr_primer key, in the format 'pcr_primers, csv_file'
    PARSE_INPUTS.out.samplesheet_locus
        .flatten()
        .map { csv -> 
            def csv_name = csv.getFileName().toString()
            ( pcr_primers, rest ) = csv_name.split("__")
            [ pcr_primers, csv ]
            }
        // .dump (tag: 'ch_loci_samdf')
        .set { ch_loci_samdf }  

    //// get names and count of the multiplexed loci used
    PARSE_INPUTS.out.samplesheet_locus
        .flatten ()
        .splitCsv ( header: true )
        .map { row -> row.target_gene }
        .unique ()
        .toList ()
        .set { ch_loci_names } // value channel; list

    ch_loci_names
        .flatten ()
        .count ()
        .set { ch_loci_number } // value channel; integer


    //// get names of flow cells ('fcid') as channel
    //// TODO: Move these to outputs of PARSE_INPUTS, to tidy up the pipeline logic
    // extract fcid from metadata
    ch_sample_locus_reads 
        .map { meta, reads ->
            def fcid = meta.fcid
            return tuple(fcid, reads) } 
        .groupTuple() 
        .set { ch_sample_reads_fcid }
    
    // extract flow cell IDs as channel
    ch_sample_reads_fcid 
        .map { group -> group[0] }
        .set { ch_fcid }

    // create channel linking pcr_primers and databases (from params)
    ch_sample_locus_reads 
        .map { meta, reads -> 
                [ meta.pcr_primers, meta.target_gene, meta.idtaxa_db, meta.ref_fasta ] }
        .unique()
        .set { ch_loci_info }

    //// create channel of loci parameters
    PARSE_INPUTS.out.loci_params_parsed
        .splitCsv ( header: true )
        .map { row -> 
                [ row.pcr_primers, row ] }
        .set { ch_loci_params } // cardinality: pcr_primers, map(all params, incl. pcr_primers)


    //// subworkflow: process sequencing reads
    PROCESS_READS (
        ch_sample_locus_reads,
        params.seq_type,
        params.paired,
        ch_fcid
        )


    //// subworkflow: run DADA2 to infer ASVs
    DADA2 (
        PROCESS_READS.out.ch_processed_reads,
        ch_read_tracker_grouped
    )


    //// subworkflow: assign taxonomy
    TAXONOMY (
        DADA2.out.ch_seqtab_meta,
        ch_loci_params
    )
    
    // //// inputs for PHYLOSEQ_UNFILTERED
    // ch_taxtables_locus = 
    //     MERGE_TAX.out.merged_tax // pcr_primers, path("*_merged_tax.rds")

    // ch_seqtables_locus = 
    //     ch_seqtab
    //     .map { pcr_primers, fcid, seqtab -> [ pcr_primers, seqtab ] } // remove fcid field
    //     .groupTuple ( by: 0 ) // group seqtabs into lists per locus
    
    // // combine taxtables, seqtables and parameters
    // ch_taxtables_locus
    //     // .dump(tag: '1')
    //     .combine ( ch_seqtables_locus, by: 0 )
    //     // .dump(tag: '2')
    //     .combine ( ch_loci_samdf, by: 0 )
    //     // .dump(tag: '3')
    //     .combine ( ch_loci_params, by: 0 )
    //     // .dump(tag: '4')
    //     .set { ch_phyloseq_input }

    
    // //// create phyloseq objects per locus; output unfiltered summary tables and accumulation curve plot
    // PHYLOSEQ_UNFILTERED ( ch_phyloseq_input )

    // //// apply taxonomic and minimum abundance filtering per locus (from loci_params), then combine to output filtered summary tables
    // PHYLOSEQ_FILTER ( PHYLOSEQ_UNFILTERED.out.ps )

    // //// combine phyloseq outputs to merge across loci
    // PHYLOSEQ_UNFILTERED.out.ps // val(pcr_primers), path("ps_unfiltered_*.rds"), val(loci_params)
    //     .map { pcr_primers, ps, loci_params ->
    //         [ ps ] }
    //     .collect()
    //     .set { ch_ps_unfiltered }

    // PHYLOSEQ_FILTER.out.ps // val(pcr_primers), path("ps_filtered_*.rds"), val(loci_params)
    //     .map { pcr_primers, ps, loci_params ->
    //         [ ps ] }
    //     .collect()
    //     .set { ch_ps_filtered }

    // PHYLOSEQ_MERGE ( 
    //     ch_ps_unfiltered, 
    //     ch_ps_filtered 
    //     )
    // ch_read_tracker_grouped = 
    //     ch_read_tracker_grouped
    //     .concat( PHYLOSEQ_MERGE.out.read_tracking.flatten() ) 
    //     .collect()

    //// track reads and sequences across the pipeline
    // READ_TRACKING ( 
    //     PROCESS_READS.out.ch_read_tracker_samples.collect(), 
    //     ch_read_tracker_grouped 
    //     )


    ///// VISUALISATION
    /*
    - stacked bar charts of OTU abundance across samples -- filtered vs unfiltered
    - plots showing proportion of sequences that passed or failed the final filters, per sample/flowcell/locus etc
    - taxonomy heatmap tree
    */

    /// perhaps a process called "REPORTS()" that generates HTML reports of the pipeline output


    //// move final output files to the output report directory

    // make directories for unfiltered and filtered outputs

    //// process called "OUTPUT()" that collates all the relevant output files and copies them into particular directories
}

workflow {
    FREYR()
}