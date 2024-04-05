process PARAMETER_SETUP {
    def module_name = "parameter_setup"
    tag "Whole dataset"
    // label: 

    input:


    output:  
    path("params.csv") ,                emit: loci_params
    path("samdf_*.csv") ,               emit: samdf_locus
    // path("fastq_paths_all.csv"),        emit: fastq_paths_all
    // path("fastq_paths_samples.csv"),    emit: fastq_paths_samples
    // path("fastq_paths_ud.csv"),         emit: fastq_paths_ud
    // path("params_primer.csv"),          emit: params_primer
    // path("params_readfilter.csv"),      emit: params_readfilter
    // path("params_dada.csv"),            emit: params_dada
    // path("params_asvfilter.csv"),       emit: params_asvfilter
    // path("params_database.csv"),        emit: params_database
    // path("params_ps.csv"),              emit: params_ps


    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables

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