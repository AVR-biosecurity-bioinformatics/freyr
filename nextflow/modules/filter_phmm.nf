process FILTER_PHMM {
    def process_name = "filter_phmm"
    tag "$primers"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), path(seqtab_tibble_list), path(fasta_list), val(process_params)

    output:
    tuple val(primers), path("*_phmm_filter.csv"),      emit: tibble

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
        --seqtab_tibble_list "$seqtab_tibble_list" \
        --fasta_list "$fasta_list" \
        --phmm "$process_params.phmm" \
        --for_primer_seq "$process_params.for_primer_seq" \
        --rev_primer_seq "$process_params.rev_primer_seq" \
        --coding "$process_params.coding" 
    
    """

}