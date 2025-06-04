process MERGE_FILTERS {
    def module_name = "merge_filters"
    tag "$pcr_primers"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(pcr_primers), path(filter_tibble_list), path(seqtab_tibble_list), path(fasta_list)
    val(fcid_list)

    output:
    tuple val(pcr_primers), path("*_seqtab_combined.csv"), path("*_filters.csv"), path("*_seqs.fasta"),      emit: seqtab
    path("*_ASV_cleanup.csv"),                                                                                  emit: cleanup
    path("*_asv_abundance.pdf"),                                                                                emit: abundance_plot
    path("*_asv_count.pdf"),                                                                                    emit: count_plot
    path("*_readsout.csv"),                                                                                     emit: read_tracking

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
        
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    pcr_primers =           "${pcr_primers}"
    filter_tibble_list =    "${filter_tibble_list}"
    seqtab_tibble_list =    "${seqtab_tibble_list}"
    fasta_list =            "${fasta_list}"
    fcid_list =             "${fcid_list}"
    
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