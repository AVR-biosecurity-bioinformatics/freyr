process FILTER_PHMM {
    tag "$primers"
    label "small"
    label "env_r_phmm"
    publishDir "${params.outdir}/modules/filter_phmm", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(primers), path(fasta), val(process_params)

    output:
    tuple val(primers), path("*_phmm_filter.csv"),      emit: tibble

    script:
    """
    
    Rscript --vanilla "${moduleDir}/filter_phmm.R" \
        --process_name "filter_phmm" \
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