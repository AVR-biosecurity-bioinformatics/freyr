process FILTER_QUALPLOTS {
    def module_name = "filter_qualplots"
    tag "$sample_primers"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(reads), val(process_params)
    val(file_suffix)

    output:   
    path("*_qualplots.pdf")                 , emit: plots, optional: true

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"

    def reads_paths = reads.join(";") // concatenate reads into a single string
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    reads_paths =       "${reads_paths}"
    sample_id =         "${sample_primers}"
    fcid =              "${read_group}"
    target_gene =       "${process_params.locus}"
    pcr_primers =       "${primers}"
    seq_type =          "${process_params.seq_type}"
    paired =            "${process_params.paired}"
    file_suffix =       "${file_suffix}"
    
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