process PHYLOSEQ_CLUSTERED {
    def process_name = "phyloseq_clustered"
    tag "$primers"
    label "phyloseq"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), path(seqtab), path(taxtab), path(samdf), path(raw), path(summary), path(ps), path(clusters)
    val(merge_clusters)

    output:
    path("*.csv"),                                                          emit: csvs
    tuple val(primers), path("ps_clustered_*.rds"),                         emit: ps 
    path("asvs_clustered*.fasta"),                                          emit: asv_fasta

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy'

    // when: 

    script:
    """
    
    ${process_name}.R \
        --process_name "$process_name" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --seqtab_file "$seqtab" \
        --taxtab_file "$taxtab" \
        --samdf_file "$samdf" \
        --raw_file "$raw" \
        --summary_file "$summary" \
        --ps_file "$ps" \
        --clusters_file "$clusters" \
        --merge_clusters "$merge_clusters" 
        
    """
    
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    primers =               "${primers}"
    seqtab_file =           "${seqtab}"
    taxtab_file =           "${taxtab}"
    samdf_file =            "${samdf}"
    raw_file =              "${raw}"
    summary_file =          "${summary}"
    ps_file =               "${ps}"
    clusters_file =         "${clusters}"
    merge_clusters =        "${merge_clusters}"
     
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