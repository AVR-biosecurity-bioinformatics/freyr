/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'
// this allows 'fromSamplesheet' command to pull data from samplesheet file

// def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
// def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
// def summary_params = paramsSummaryMap(workflow)

// // Print parameter summary log to screen
// log.info logo + paramsSummaryLog(workflow) + citation

// WorkflowAmpliseq.initialise(params, log)

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT AND VARIABLES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Input

// directory in projectDir where data files are stored; "test_data" for test data
// def data_loc = "test_data" // not needed as parameters are parsed in R script

// report sources

// Set non-params Variables


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PARAMETER_SETUP                           } from '../modules/parameter_setup'
include { SEQ_QC                                    } from '../modules/seq_qc'
include { SPLIT_LOCI                                } from '../modules/split_loci'
include { PRIMER_TRIM                               } from '../modules/primer_trim'
include { READ_FILTER                               } from '../modules/read_filter' 
include { FILTER_QUALPLOTS as FILTER_QUALPLOTS_PRE  } from '../modules/filter_qualplots'
include { FILTER_QUALPLOTS as FILTER_QUALPLOTS_POST } from '../modules/filter_qualplots'
include { ERROR_MODEL as ERROR_MODEL_F              } from '../modules/error_model'
include { ERROR_MODEL as ERROR_MODEL_R              } from '../modules/error_model'
include { DENOISE as DENOISE_F                      } from '../modules/denoise'
include { DENOISE as DENOISE_R                      } from '../modules/denoise'
include { DADA_PRIORS as DADA_PRIORS_F              } from '../modules/dada_priors'
include { DADA_PRIORS as DADA_PRIORS_R              } from '../modules/dada_priors'
include { DENOISE as DENOISE2_F                     } from '../modules/denoise'
include { DENOISE as DENOISE2_R                     } from '../modules/denoise'
// include { DADA                                      } from '../modules/dada'
// include { FILTER_SEQTAB                             } from '../modules/filter_seqtab'
// include { MERGE_SEQTAB                              } from '../modules/merge_seqtab'
// include { TAX_IDTAXA                                } from '../modules/tax_idtaxa'
// include { TAX_BLAST                                 } from '../modules/tax_blast'
// include { JOINT_TAX                                 } from '../modules/joint_tax'
// include { MERGE_TAX                                 } from '../modules/merge_tax'
// include { ASSIGNMENT_PLOT                           } from '../modules/assignment_plot'
// include { TAX_SUMMARY                               } from '../modules/tax_summary'
// include { PHYLOSEQ_CREATE                           } from '../modules/phyloseq_create'
// include { PHYLOSEQ_SUMMARY                          } from '../modules/phyloseq_summary'
// include { ACCUMULATION_CURVE                        } from '../modules/accumulation_curve'
// include { PHYLOSEQ_FILTER                           } from '../modules/phyloseq_filter'
// include { READ_TRACKING                             } from '../modules/read_tracking'




//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

// include { PARSE_INPUT                   } from '../subworkflows/local/parse_input'


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
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPERLINE {

    ch_versions = Channel.empty()

    //
    // Create input channels
    //


    //// input samplesheet and loci parameters
    PARAMETER_SETUP ( )
    
    //// TODO: remove this bit if it's not used in final pipeline
    //// get read paths and metadata for each sample from the sample sheet .csv
    //// from: https://training.nextflow.io/advanced/grouping/#grouping-using-submap 
    PARAMETER_SETUP.out.samdf
        | splitCsv ( header: true )
        | map { row -> 
            meta = row.subMap(
                'sample_id','sample_name','extraction_rep','amp_rep',
                'client_name','experiment_name','sample_type','collection_method',
                'collection_location','lat_lon','environment','collection_date',
                'operator_name','description','assay','extraction_method',
                'amp_method','target_gene','pcr_primers','for_primer_seq',
                'rev_primer_seq','index_plate','index_well','i7_index_id',
                'i7_index','i5_index_id','i5_index','seq_platform',
                'fcid','for_read_length','rev_read_length','seq_run_id',
                'seq_id','seq_date','analysis_method','notes','base'
                )
                [ meta, [
                    file(row.fwd, checkIfExists: true),
                    file(row.rev, checkIfExists: true)
                ]]  
            }
        | set { ch_sample_reads }

    // ch_sample_reads | view // check output

    //// parse samplesheets that contain locus-specific parameters
    PARAMETER_SETUP.out.samdf_locus
        | flatten ()
        | splitCsv ( header: true )
        | map { row -> 
            meta = row.subMap(
                'sample_id','sample_name','extraction_rep','amp_rep',
                'client_name','experiment_name','sample_type','collection_method',
                'collection_location','lat_lon','environment','collection_date',
                'operator_name','description','assay','extraction_method',
                'amp_method','target_gene','pcr_primers','for_primer_seq',
                'rev_primer_seq','index_plate','index_well','i7_index_id',
                'i7_index','i5_index_id','i5_index','seq_platform',
                'fcid','for_read_length','rev_read_length','seq_run_id',
                'seq_id','seq_date','analysis_method','notes',
                'base','max_primer_mismatch','read_min_length','read_max_length',
                'read_max_ee','read_trunc_length','read_trim_left','read_trim_right',
                'asv_min_length','asv_max_length','genetic_code','coding',
                'phmm','idtaxa_db','ref_fasta','idtaxa_confidence',
                'run_blast','blast_min_identity','blast_min_coverage','target_kingdom',
                'target_phylum','target_class','target_order','target_family',
                'target_genus','target_species','min_sample_reads','min_taxa_reads',
                'min_taxa_ra','threads'
                )
                [ meta, [
                    file(row.fwd, checkIfExists: true),
                    file(row.rev, checkIfExists: true)
                ]]  
            }
        | set { ch_sample_locus_reads }

        // ch_sample_locus_reads | view // check output

    //// get names and count of the multiplexed loci used
    PARAMETER_SETUP.out.samdf_locus
    | flatten ()
    | splitCsv ( header: true )
    | map { row -> row.target_gene }
    | unique ()
    | toList ()
    | set { ch_loci_names } // value channel; list
    
    ch_loci_names
    | flatten ()
    | count ()
    | set { ch_loci_number } // value channel; integer


    //// get names of flow cells ('fcid') as channel
    // extract fcid from metadata
    ch_sample_reads 
    | map { meta, reads ->
        def fcid = meta.fcid
        return tuple(fcid, reads) } 
    | groupTuple() 
    | set { ch_sample_reads_fcid }
    // extract flow cell IDs as channel
    ch_sample_reads_fcid 
    | map { group -> group[0] }
    | set { ch_fcid }


    // run SEQ_QC per flow cell 
    // SEQ_QC ( ch_fcid ) // optional step for testing

    /// TODO: develop method to count reads as they move through the pipeline
    //      ie. input file, after splitting, after primer trimming, after qual trim etc.

    //// split sample reads by locus (based on primer seq.)
    SPLIT_LOCI ( ch_sample_locus_reads ) 

    // SPLIT_LOCI.out.reads.view()

    //// trim primer sequences from the start and end of reads
    PRIMER_TRIM ( SPLIT_LOCI.out.reads )

    // PRIMER_TRIM.out.reads.view()

    //// filter reads using dada2 and input parameters
    READ_FILTER ( PRIMER_TRIM.out.reads )

    // READ_FILTER.out.reads.view()

    //// create plots of read quality pre- and post-filtering, per flowcell (optional)
    // FILTER_QUALPLOTS_PRE ( PRIMER_TRIM.out.reads )

    // FILTER_QUALPLOTS_POST ( READ_FILTER.out.reads )

    /// TODO: Use FILTER_QUALPLOTS_COMBINE to combine plots by fcid and type into one PDF


    ///// split filtered reads into lists of reads per flowcell, primers and direction
    //// forward read channel
    READ_FILTER.out.reads
    | map { meta, reads -> 
            [ "forward", meta.fcid, meta.pcr_primers, reads[0] ] }
    | groupTuple( by: [0,1,2] )
    | set { ch_error_input_fwd }

    //// reverse read channel
    READ_FILTER.out.reads
    | map { meta, reads -> 
            [ "reverse", meta.fcid, meta.pcr_primers, reads[1] ] }
    | groupTuple ( by: [0,1,2] )
    | set { ch_error_input_rev }


    //// error model on forward reads
    ERROR_MODEL_F ( ch_error_input_fwd )

    //// error model on reverse reads
    ERROR_MODEL_R ( ch_error_input_rev )

    //// input channel for denoising forward reads
    READ_FILTER.out.reads
    | map { meta, reads -> 
            [ "forward", meta.fcid, meta.pcr_primers, meta, reads[0] ] }
    | combine ( ERROR_MODEL_F.out.errormodel, by: [0,1,2] ) // combine with error model
    | map { direction, fcid, pcr_primers, meta, reads, errormodel -> // add empty file path for priors
            [ direction, fcid, pcr_primers, meta, reads, errormodel, "$projectDir/assets/NO_FILE" ] }
    | set { ch_denoise_input_forward }

    // ch_denoise_input_forward | view ()

    //// input channel for denoising reverse reads
    READ_FILTER.out.reads
    | map { meta, reads -> 
            [ "reverse", meta.fcid, meta.pcr_primers, meta, reads[1] ] }
    | combine ( ERROR_MODEL_R.out.errormodel, by: [0,1,2] ) // combine with error model
    | map { direction, fcid, pcr_primers, meta, reads, errormodel -> // add empty file path for priors
            [ direction, fcid, pcr_primers, meta, reads, errormodel, "$projectDir/assets/NO_FILE" ] }
    | set { ch_denoise_input_reverse }


    //// first pass of denoising (always done)
    // channels for denoising
    ch_firstpass = Channel.value ( "first" )
    ch_secondpass = Channel.value ( "second" )

    //// denoise forward reads per flowcell, primer and sample
    DENOISE_F ( ch_denoise_input_forward, ch_firstpass )
    
    //// denoise forward reads per flowcell, primer and sample
    DENOISE_R ( ch_denoise_input_reverse, ch_firstpass )

    // high sensitivity mode condition
    if ( params.high_sensitivity ) { // run prior extraction and second pass denoising
        //// group priors for each read file
        /// forward reads
        DENOISE_F.out.seq
        | map { direction, fcid, pcr_primers, meta, reads, priors ->
                [ direction, fcid, pcr_primers, priors ] }
        | groupTuple ( by: [0,1,2] )
        | set { ch_priors_f_pre }

        /// reverse reads
        DENOISE_R.out.seq
        | map { direction, fcid, pcr_primers, meta, reads, priors ->
                [ direction, fcid, pcr_primers, priors ] }
        | groupTuple ( by: [0,1,2] )
        | set { ch_priors_r_pre }

        //// get priors for forward reads
        DADA_PRIORS_F ( ch_priors_f_pre )
        
        /// combine with forward read data channel
        ch_denoise_input_forward
        | map { direction, fcid, pcr_primers, meta, reads, errormodel, priors -> // remove null priors
                [ direction, fcid, pcr_primers, meta, reads, errormodel ] }
        | combine ( DADA_PRIORS_F.out.priors, by: [0,1,2] )
        | set { ch_denoise2_input_forward }

        //// get priors for reverse reads
        DADA_PRIORS_R ( ch_priors_r_pre )

        /// combine with reverse read data channel
        ch_denoise_input_reverse
        | map { direction, fcid, pcr_primers, meta, reads, errormodel, priors -> // remove null priors
                [ direction, fcid, pcr_primers, meta, reads, errormodel ] }
        | combine ( DADA_PRIORS_R.out.priors, by: [0,1,2] )
        | set { ch_denoise2_input_reverse }

        //// run pseudo-pooled 2nd denoising with priors on forward reads
        DENOISE2_F ( ch_denoise2_input_forward, ch_secondpass )

        //// run pseudo-pooled 2nd denoising with priors on reverse reads
        DENOISE2_R ( ch_denoise2_input_reverse, ch_secondpass )

        //// set output as input for merging and seqtab construction
        /// join F and R denoise2 outputs
        // prepare forward reads
        DENOISE2_F.out.seq
        | map { direction, fcid, pcr_primers, meta, reads, seqF ->
                [ meta.sample_id, fcid, pcr_primers, meta, seqF ] }
        | set { ch_seq_forward }

        // prepare reverse reads
        DENOISE2_R.out.seq
        | map { direction, fcid, pcr_primers, meta, reads, seqR ->
                [ meta.sample_id, fcid, pcr_primers, meta, seqR ] }
        | set { ch_seq_reverse }

        // join
        ch_seq_forward
        | combine ( ch_seq_reverse, by: [0,1,2,3] ) 
        | view ()

    } else { // don't run second denoising

        //// set output as input for merging and seqtab construction

    }





    //// denoise forward reads per flowcell, primer and sample
    // DENOISE_R ( ch_denoise_input_reverse, ch_high_sen_mode, ch_nopriors )

    /// Condition for processing reads in "high sensitivity mode", ie. with priors and pseudo-pooled 2nd denoising
    
    
    

    
    

    
    
    
    
    // }





    
    
}
