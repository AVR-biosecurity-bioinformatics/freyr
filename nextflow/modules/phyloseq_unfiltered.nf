process PHYLOSEQ_UNFILTERED {
    def module_name = "phyloseq_unfiltered"
    tag "$primers"
    label "phyloseq"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), path(taxtab), path(seqtab), path(filters), path(fasta), path(samplesheet_split), path(sample_metadata), val(process_params)

    output:
    tuple val(primers), path("seqtab_unfiltered_*.csv"), path("taxtab_unfiltered_*.csv"), path("samdf_unfiltered_*.csv"), path("raw_unfiltered_*.csv"), path("summary_unfiltered_*.csv"), emit: csvs
    tuple val(primers), path("ps_unfiltered_*.rds"), path("filters_*.csv"),                        emit: ps 
    path("asvs_unfiltered_*.fasta"),                                                               emit: asv_fasta

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    primers =                   "${primers}"
    taxtab_file =               "${taxtab}"
    seqtab_file =               "${seqtab}"
    filters_file =              "${filters}"
    fasta_file =                "${fasta}"
    samplesheet_split_file =    "${samplesheet_split}"
    sample_metadata_file =      "${sample_metadata}"
    cluster_threshold =         "${process_params.cluster_threshold}"
     
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