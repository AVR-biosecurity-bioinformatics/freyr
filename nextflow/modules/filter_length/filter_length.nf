process FILTER_LENGTH {
    tag "$primers"
    label "medium"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/filter_length", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(primers), path(seqtab_tibble_list), path(fasta_list), val(process_params)

    output:
    tuple val(primers), path("*_length_filter.csv"),            emit: tibble


    script:
    """
    
    filter_length.R \
        --process_name "filter_length" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --seqtab_tibble_list "$seqtab_tibble_list" \
        --fasta_list "$fasta_list" \
        --asv_min_length "$process_params.asv_min_length" \
        --asv_max_length "$process_params.asv_max_length"
    
    """

}