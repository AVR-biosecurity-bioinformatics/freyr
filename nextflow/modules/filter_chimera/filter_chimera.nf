process FILTER_CHIMERA {
    tag "$primers"
    label "high"
    label "env_r_dada2"
    publishDir { "${params.outdir}/modules/filter_chimera" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()

    input:
    tuple val(primers), path(seqtab_tibble_list), path(fasta_list), val(process_params)

    output:
    tuple val(primers), path("*_chimera_filter.csv"),        emit: tibble

    script:
    """
    
    Rscript --vanilla "${moduleDir}/filter_chimera.R" \
        --process_name="filter_chimera" \
        --cpus="$task.cpus" \
        --rdata="$params.rdata" \
        --primers="$primers" \
        --seqtab_tibble_list="$seqtab_tibble_list" \
        --fasta_list="$fasta_list" \
        --minSampleFraction="$process_params.chimera_sample_frac"
    
    """

}