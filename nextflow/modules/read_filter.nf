process READ_FILTER {
    def module_name = "read_filter"
    tag "$sample_primers"
    label "medium"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(reads), val(process_params)

    output:   
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path("*_filter_R{0,1,2}.fastq.gz"),      emit: reads, optional: true
    path("*_readsout.csv"),                                                                                         emit: read_tracking

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
    read_min_length =   "${process_params.read_min_length}"
    read_max_length =   "${process_params.read_max_length}"
    read_max_ee =       "${process_params.read_max_ee}"
    read_trunc_length = "${process_params.read_trunc_length}"
    read_trim_left =    "${process_params.read_trim_left}"
    read_trim_right =   "${process_params.read_trim_right}"
    sample_primers =    "${sample_primers}"
    locus =             "${process_params.locus}"
    primers =           "${primers}"
    read_group =        "${read_group}"
    seq_type =          "${process_params.seq_type}"
    paired =            "${process_params.paired}"
    
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