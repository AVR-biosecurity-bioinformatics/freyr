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

// check Nextflow version matches 23.05.0-edge, due to incompatibility of higher versions with Shifter
if( !nextflow.version.matches('=23.05.0-edge') ) {
    println " "
    println "*** ERROR ~ This pipeline currently requires Nextflow version 23.05.0-edge -- You are running version ${nextflow.version}. ***"
    error "*** Either: append 'NXF_VER=23.05.0-edge' to the front of the 'nextflow run' command, or run 'export NXF_VER=23.05.0-edge'. ***"
}

// Validate input parameters using schema
validateParameters( parameters_schema: 'nextflow_schema.json' )

// Print summary of supplied parameters (that differ from defaults)
log.info paramsSummaryLog(workflow)
// make it clear that samples are being subsampled
if (params.subsample) {
    log.info "***"
    log.info "NOTE: Input samples are being randomly subsampled to $params.subsample per primer x flowcell combination (params.subsample = $params.subsample)"
    log.info "***"
}

if (params.downsample_reads) {
    log.info "***"
    log.info "NOTE: Input samples are being randomly downsampled to $params.downsample_reads reads (params.subsample = $params.downsample_reads)"
    log.info "***"
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
include { DOWNSAMPLE_READS                          } from './nextflow/modules/downsample_reads'
include { TRAIN_IDTAXA                              } from './nextflow/modules/train_idtaxa'

//// utility processes for development and debugging
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

    //// parse path channels
    ch_samplesheet_file = channel.fromPath( params.samplesheet, checkIfExists: true, type: 'file' )
    ch_loci_params_file = channel.fromPath( params.loci_params, checkIfExists: true, type: 'file' )

    //// read-in samplesheet and loci_params .csv files, validate their contents, and produce inputs for rest of pipeline
    PARSE_INPUTS ( 
        ch_samplesheet_file, 
        ch_loci_params_file,
        params.seq_type,
        params.paired
        )

    // ch_versions = Channel.empty()
    
    //// Create empty channels
    ch_read_tracker_grouped =   // read-tracking for grouped processes
        Channel.empty()      
    ch_idtaxa_db_new = 
        Channel.empty()


    //// parse samplesheets that contain locus-specific parameters
    if ( params.paired == true ) {
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
            [ meta, [ file(row.fwd, checkIfExists: true), file(row.rev, checkIfExists: true) ] ]  
            }
        .set { ch_sample_locus_reads }
    } else if ( params.paired == false ) {
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
            [ meta, file(row.single, checkIfExists: true) ]  
            }
        .set { ch_sample_locus_reads }
    } else {
        error " 'params.paired' must be 'true' or 'false'. "
    }

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


    //// train IDTAXA model from reference database .fasta
    if ( params.train_idtaxa ) {
        
        //// create input channel for TRAIN_IDTAXA
        ch_loci_params
            .map { pcr_primers, loci_params ->
                [ pcr_primers, loci_params.ref_fasta ]  }
            .set { ch_train_idtaxa_input }

        //// train model
        TRAIN_IDTAXA (
            ch_train_idtaxa_input,
            params.train_idtaxa
        )

        ch_idtaxa_db_new
            .concat ( TRAIN_IDTAXA.out.model )
            .set { ch_idtaxa_db_new }
    } else {
        ch_idtaxa_db_new = channel.value( "no_new_model" )
    }

    //// downsample reads if params.downsample is defined
    if ( params.downsample_reads ) {
        DOWNSAMPLE_READS (
            ch_sample_locus_reads,
            params.seq_type,
            params.paired,
            params.downsample_reads
        )

        ch_sample_locus_reads = DOWNSAMPLE_READS.out.reads
    }

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
        DADA2.out.ch_seqtab,
        ch_loci_params,
        ch_idtaxa_db_new
        )


    //// subworkflow: create result summaries
    RESULT_SUMMARIES (
        DADA2.out.ch_seqtab,
        TAXONOMY.out.ch_mergetax_output,
        TAXONOMY.out.ch_seqtab,
        ch_loci_samdf,
        ch_loci_params,
        PROCESS_READS.out.ch_read_tracker_samples,
        DADA2.out.read_tracker_grouped
        )
    

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