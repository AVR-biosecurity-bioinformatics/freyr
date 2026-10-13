process FILTER_PHMM {
    tag "$primers"
    label "small"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/filter_phmm", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(primers), path(fasta), val(process_params)

    output:
    tuple val(primers), path("*_phmm_filter.csv"),      emit: tibble

    script:
    """
    
    ${moduleDir}/filter_phmm.R \
        --process_name "filter_phmm" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --fasta "$fasta" \
        --phmm "$process_params.phmm" \
        --for_primer_seq "$process_params.for_primer_seq" \
        --rev_primer_seq "$process_params.rev_primer_seq" \
        --coding "$process_params.coding" 
    
    """

}