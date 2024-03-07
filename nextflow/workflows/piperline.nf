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
def data_loc = "test_data"

// report sources

// Set non-params Variables


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { RENAME_RAW_DATA_FILES         } from '../modules/local/rename_raw_data_files'

include { PARAMETER_SETUP                   }         from '../modules/parameter_setup'
// include { SEQ_QC                            }         from '../modules/seq_qc'
// include { SWITCHING_QC                      }         from '../modules/switching_qc'
// include { PRIMER_TRIM                       }         from '../modules/primer_trim'
// include { SAMDF2                            }         from '../modules/samdf2'
// include { READ_FILTER                       }         from '../modules/read_filter' 
// include { SAMDF3                            }         from '../modules/samdf3'
// include { PREFILT_QUALPLOTS                 }         from '../modules/prefilt_qualplots'
// include { POSTFILT_QUALPLOTS                }         from '../modules/postfilt_qualplots'
// include { ERROR_MODEL as ERROR_MODEL_F      }         from '../modules/error_model'
// include { ERROR_MODEL as ERROR_MODEL_R      }         from '../modules/error_model'
// include { DENOISE as DENOISE_F              }         from '../modules/denoise'
// include { DENOISE as DENOISE_R              }         from '../modules/denoise'
// include { DENOISE as DENOISE2_F             }         from '../modules/denoise'
// include { DENOISE as DENOISE2_R             }         from '../modules/denoise'
// include { DADA                              }         from '../modules/dada'
// include { FILTER_SEQTAB                     }         from '../modules/filter_seqtab'
// include { MERGE_SEQTAB                      }         from '../modules/merge_seqtab'
// include { TAX_IDTAXA                        }         from '../modules/tax_idtaxa'
// include { TAX_BLAST                         }         from '../modules/tax_blast'
// include { JOINT_TAX                         }         from '../modules/joint_tax'
// include { MERGE_TAX                         }         from '../modules/merge_tax'
// include { ASSIGNMENT_PLOT                   }         from '../modules/assignment_plot'
// include { TAX_SUMMARY                       }         from '../modules/tax_summary'
// include { PHYLOSEQ_CREATE                   }         from '../modules/phyloseq_create'
// include { PHYLOSEQ_SUMMARY                  }         from '../modules/phyloseq_summary'
// include { ACCUMULATION_CURVE                }         from '../modules/accumulation_curve'
// include { PHYLOSEQ_FILTER                   }         from '../modules/phyloseq_filter'
// include { READ_TRACKING                     }         from '../modules/read_tracking'




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

    // def file_samdf = "${projectDir}/sample_data/Sample_info.csv"
    // def file_params = "${projectDir}/sample_data/loci_params.csv"

    //file_samdf = file("${projectDir}/sample_data/Sample_info.csv",  checkIfExists: true)
    //file_params = file("${projectDir}/sample_data/loci_params.csv",  checkIfExists: true)
    // file_samdf.view { "path: $it" }

    // // for debugging if paths exist
    // ch_samdf.view { "path: $it" }
    // ch_loci_params.view { "path: $it" }

    // input samplesheet and loci parameters
    PARAMETER_SETUP ( data_loc )
    
    // get read paths and metadata for each sample from the sample sheet .csv
    // from: https://training.nextflow.io/advanced/grouping/#grouping-using-submap 
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
        | set { ch_reads_meta }

    ch_reads_meta | view

    // this takes 
    // SEQ_QC ()


}
