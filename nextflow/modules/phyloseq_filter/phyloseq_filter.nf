process PHYLOSEQ_FILTER {
    tag "$primers"
    label "medium"
    label "env_r_phyloseq"
    publishDir { "${params.outdir}/modules/phyloseq_filter" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()
    publishDir { "${params.outdir}/results/main_outputs/filtered" }, mode: 'copy'

    input:
    tuple val(primers), path(ps), path(filters_tibble), val(process_params)

    output:
    tuple val(primers), path("seqtab_filtered_*.csv"), path("taxtab_filtered_*.csv"), path("samdf_filtered_*.csv"), path("raw_filtered_*.csv"), path("summary_filtered_*.csv"), emit: csvs
    tuple val(primers), path("ps_filtered_*.rds"),                          emit: ps 
    path("*.fasta"),                                                        emit: asv_fasta
    tuple val(primers), path("clusters_*.csv"),                             emit: clusters

    script:
    """
     
    Rscript --vanilla "${moduleDir}/phyloseq_filter.R" \
        --process_name="phyloseq_filter" \
        --cpus="$task.cpus" \
        --rdata="$params.rdata" \
        --primers="$primers" \
        --ps="$ps" \
        --filters_tibble="$filters_tibble" \
        --cluster_threshold="$process_params.cluster_threshold" \
        --target_kingdom="$process_params.target_kingdom" \
        --target_phylum="$process_params.target_phylum" \
        --target_class="$process_params.target_class" \
        --target_order="$process_params.target_order" \
        --target_family="$process_params.target_family" \
        --target_genus="$process_params.target_genus" \
        --target_species="$process_params.target_species" \
        --min_sample_reads="$process_params.min_sample_reads" \
        --min_taxa_reads="$process_params.min_taxa_reads" \
        --min_taxa_ra="$process_params.min_taxa_ra" 
        
    """

}