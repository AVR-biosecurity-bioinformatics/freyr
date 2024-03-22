process SEQ_QC {
    def module_name = "seq_qc"
    tag: flowcell_id
    // label: 

    input:
    // val data_loc
    val flowcell_id 

    output:
    path "${flowcell_id}_flowcell_qc.pdf"
    path "${flowcell_id}_index_switching.pdf" 

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when:

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    flowcell_id = "$flowcell_id"
    
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
    
    """
}