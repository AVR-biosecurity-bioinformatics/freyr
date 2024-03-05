process SEQ_QC {
    def module_name = "seq_qc"
    // tag:
    // label: 

    input:
    path sample_info
    path loci_params

    output:

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'


    // when:

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    #set working directory to projectDir so everything works
    setwd("${projectDir}")

    # load R libraries
    source("running_scripts/_targets_packages.R") # make sure this is the current location
    # also check that I can't just load a few packages like purrr, dplyr and tibble to save loading all the packages in the targets pipeline
    source("${functions_r}") # this should source the "functions.R" script
    source("${themes_r}") # this should source the "themes.R" script


    """
}