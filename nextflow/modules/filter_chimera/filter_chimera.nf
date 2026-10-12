process FILTER_CHIMERA {
    tag "$primers"
    label "very_high"
    container "jackscanlan/piperline-multi:0.0.1"
    publishDir "${launchDir}/output/modules/filter_chimera", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(primers), path(seqtab_tibble_list), path(fasta_list), val(process_params)

    output:
    tuple val(primers), path("*_chimera_filter.csv"),        emit: tibble

    script:
    """
    
    filter_chimera.R \
        --process_name "filter_chimera" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --seqtab_tibble_list "$seqtab_tibble_list" \
        --fasta_list "$fasta_list" \
        --minSampleFraction "$process_params.chimera_sample_frac"
    
    """

}