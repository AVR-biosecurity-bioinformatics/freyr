process FILTER_FRAME {
    tag "$primers"
    label "medium"
    label "env_r_phyloseq"
    publishDir "${params.outdir}/modules/filter_frame", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(primers), path(seqtab_tibble_list), path(fasta_list), val(process_params)

    output:
    tuple val(primers), path("*_frame_filter.csv"),       emit: tibble

    script:
    """
    
    Rscript --vanilla "${moduleDir}/filter_frame.R" \
        --process_name "filter_frame" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --seqtab_tibble_list "$seqtab_tibble_list" \
        --fasta_list "$fasta_list" \
        --coding "$process_params.coding" \
        --genetic_code "$process_params.genetic_code"
    
    """

}