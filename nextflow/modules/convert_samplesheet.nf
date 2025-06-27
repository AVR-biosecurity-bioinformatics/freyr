process PARSE_INPUTS {
    def module_name = "parse_inputs"
    tag "Whole dataset"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(samplesheet)
    path(primer_params)
    val(seq_type)
    val(paired)
    val(subsample)

    output: 
    path("sample_metadata.csv"),                    emit: sample_metadata
    path("samplesheet_unsplit.csv"),                emit: samplesheet_unsplit
    path("samplesheet_split.csv"),                  emit: samplesheet_split
    path("primer_params_parsed.csv"),               emit: primer_params_parsed

    publishDir "${projectDir}/output/modules/${module_name}",  mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    samplesheet =           "${samplesheet}"
    primer_params =         "${primer_params}"
    seq_type =              "${seq_type}"
    paired =                "${paired}"
    subsample =             "${subsample}"

    ## global variables
    launchDir = "$launchDir"
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