process MISEQ_QC {
    def module_name = "miseq_qc"
    tag "$fcid"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    val(fcid) 
    val(miseq_dir)

    output:
    path "*_flowcell_qc.pdf"
    path "*_index_switching.pdf" 
    path "*.csv"

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when:

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fcid = "${fcid}"
    miseq_dir = "${miseq_dir}"
    
    ## global variables
    projectDir = "$projectDir"
    params_dict = "$params"
    
    tryCatch({
    ### source functions and themes, load packages, and import Nextflow params
    ### from "bin/process_start.R"
    sys.source("${projectDir}/bin/process_start.R", envir = .GlobalEnv)

    ### run module code
    sys.source(
        "${projectDir}/bin/$module_script", # run script
        envir = .GlobalEnv # this allows import of existing objects like projectDir
    )
    }, finally = {
    ### save R environment for debugging
    if ("${params.rdata}" == "true") { save.image(file = "${task.process}.rda") } 
    })

    """
}