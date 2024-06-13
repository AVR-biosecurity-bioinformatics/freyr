process SEQ_QC {
    def module_name = "seq_qc"
    tag "$fcid"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    val(fcid) 

    output:
    path "${fcid}_flowcell_qc.pdf"
    path "${fcid}_index_switching.pdf" 

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when:

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fcid = "$fcid"
    
    ## global variables
    projectDir = "$projectDir"
    params_dict = "$params"
    
    ### source functions and themes, load packages, and import Nextflow params
    ### from "bin/process_start.R"
    sys.source("${projectDir}/bin/process_start.R", envir = .GlobalEnv)

    ### run module code
    sys.source(
        "${projectDir}/bin/$module_script", # run script
        envir = .GlobalEnv # this allows import of existing objects like projectDir
    )
    
    ### save .RData for debugging
    if ("${params.rdata}" == "true") {
        save.image()
    } else {
        NULL
    }

    """
}