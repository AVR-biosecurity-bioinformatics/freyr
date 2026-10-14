process MERGE_FILTERS {
    tag "$primers"
    label "medium"
    label "env_r_core"
    publishDir "${params.outdir}/modules/merge_filters", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${params.outdir}/results/sequence_filtering", pattern: '{asv_cleanup,asv_abundance,asv_count}_*.{csv,pdf}', mode: 'copy'

    input:
    tuple val(primers), path(filter_tibble_list, name: 'input_filter*.csv'), path(seqtab_tibble_list), path(fasta_list, name: 'input_seqs*.fasta')
    path(samplesheet_split)

    output:
    tuple val(primers), path("*_seqtab_combined.csv"), path("*_filters.csv"), path("*_combseqs.fasta"),         emit: filtered
    path("{asv_cleanup,asv_abundance,asv_count}_*.{csv,pdf}"),                                                  emit: cleanup
    path("*_readsout.csv"),                                                                                     emit: read_tracking

    script:
    """
    
    Rscript --vanilla "${moduleDir}/merge_filters.R" \
        --process_name "merge_filters" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --filter_tibble_list "$filter_tibble_list" \
        --seqtab_tibble_list "$seqtab_tibble_list" \
        --fasta_list "$fasta_list" \
        --samplesheet_split "$samplesheet_split" 
        
    """

}