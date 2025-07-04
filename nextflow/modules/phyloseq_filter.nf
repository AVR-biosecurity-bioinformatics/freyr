process PHYLOSEQ_FILTER {
    def process_name = "phyloseq_filter"
    tag "$primers"
    label "phyloseq"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), path(ps), path(filters_tibble), val(process_params)

    output:
    tuple val(primers), path("seqtab_filtered_*.csv"), path("taxtab_filtered_*.csv"), path("samdf_filtered_*.csv"), path("raw_filtered_*.csv"), path("summary_filtered_*.csv"), emit: csvs
    tuple val(primers), path("ps_filtered_*.rds"),                          emit: ps 
    path("*.fasta"),                                                        emit: asv_fasta
    tuple val(primers), path("clusters_*.csv"), emit: clusters

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
        --ps "$ps" \
        --filters_tibble "$filters_tibble" \
        --cluster_threshold "$process_params.cluster_threshold" \
        --target_kingdom "$process_params.target_kingdom" \
        --target_phylum "$process_params.target_phylum" \
        --target_class "$process_params.target_class" \
        --target_order "$process_params.target_order" \
        --target_family "$process_params.target_family" \
        --target_genus "$process_params.target_genus" \
        --target_species "$process_params.target_species" \
        --min_sample_reads "$process_params.min_sample_reads" \
        --min_taxa_reads "$process_params.min_taxa_reads" \
        --min_taxa_ra "$process_params.min_taxa_ra" 
        
    """

}