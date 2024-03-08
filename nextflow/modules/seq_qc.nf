process SEQ_QC {
    def module_name = "seq_qc"
    // tag:
    // label: 

    input:
    val data_loc
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
    # defining Nextflow environment variables as R variables
    projectDir = "$projectDir"
    data_loc = "$data_loc"
    flowcell_id = "$flowcell_id"

    ### source functions and themes, and load packages from "bin/process_start.R"
    sys.source("${projectDir}/bin/process_start.R", envir = .GlobalEnv)

    ### run module code
    sys.source(
        "${projectDir}/bin/$module_script", # run script
        envir = .GlobalEnv # this allows import of existing objects like projectDir
        )
    """
}