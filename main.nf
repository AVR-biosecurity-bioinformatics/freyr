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
    Imports & Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 

//// workflow import
include { FREYR  } from './nextflow/workflows/freyr'


// function imports
include {
    paramsHelp
    validateParameters
    paramsSummaryLog
} from 'plugin/nf-schema'

// Function definitions
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    startupMessage()

    // Print help message, supply typical command line usage for the pipeline
    if (params.help) {
        def help_message = paramsHelp() as String
        log.info help_message
        exit 0, ''
    }

    // validate input params
    validateParameters()
    
    // Print summary of supplied parameters (that differ from defaults)
    log.info(
        paramsSummaryLog() as String
    )

    workflow.onComplete = {
        if (workflow.success) {
            log.info "[${workflow.complete}] >> Pipeline finished SUCCESSFULLY after ${workflow.duration}"
        } else {
            log.info "[${workflow.complete}] >> Pipeline finished with ERRORS after ${workflow.duration}"
        }
        /*
        TODO: use other metadata (https://www.nextflow.io/docs/latest/metadata.html) to display at the end of a run.
        */
    }

    // make it clear that samples are being subsampled
    if (params.subsample) {
        log.info "***"
        log.info "NOTE: Input samples are being pseudorandomly subsampled to $params.subsample per primer x read group combination (--subsample = $params.subsample)"
        log.info "***"
    }

    if (params.downsample_reads) {
        log.info "***"
        log.info "NOTE: Input samples are being pseudorandomly downsampled to $params.downsample_reads reads (--downsample_reads = $params.downsample_reads)"
        log.info "***"
    }

    // Run freyr workflow
    FREYR()
}
