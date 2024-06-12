process PRIMER_TRIM_R {
    def module_name = "primer_trim_r"
    // tag: 
    // label: 
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(meta), path(reads)

    val data_loc
    
    output: 

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when:

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    # defining Nextflow environment variables as R variables
    projectDir = "$projectDir"
    data_loc = "$params.data_folder"
    sample_id = "$meta.sample_id"
    for_primer_seq = "$meta.for_primer_seq"
    rev_primer_seq = "$meta.rev_primer_seq"
    pcr_primers = "$meta.pcr_primers"
    fcid = "$meta.fcid"

    ### source functions and themes, and load packages from "bin/process_start.R"
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