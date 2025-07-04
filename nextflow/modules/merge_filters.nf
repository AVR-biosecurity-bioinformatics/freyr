process MERGE_FILTERS {
    def process_name = "merge_filters"
    tag "$primers"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), path(filter_tibble_list), path(seqtab_tibble_list), path(fasta_list)
    path(samplesheet_split)

    output:
    tuple val(primers), path("*_seqtab_combined.csv"), path("*_filters.csv"), path("*_combseqs.fasta"),         emit: filtered
    path("*_ASV_cleanup.csv"),                                                                                  emit: cleanup
    path("*_asv_abundance.pdf"),                                                                                emit: abundance_plot
    path("*_asv_count.pdf"),                                                                                    emit: count_plot
    path("*_readsout.csv"),                                                                                     emit: read_tracking

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    // when: 

    script:
    """
    
    ${process_name}.R \
        --process_name "$process_name" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --filter_tibble_list "$filter_tibble_list" \
        --seqtab_tibble_list "$seqtab_tibble_list" \
        --fasta_list "$fasta_list" \
        --samplesheet_split "$samplesheet_split" 
        
    """

}