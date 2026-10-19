process PHYLOSEQ_CLUSTERED {
    tag "$primers"
    label "medium"
    label "env_r_phyloseq"
    publishDir { "${params.outdir}/modules/phyloseq_clustered" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()
    publishDir { "${params.outdir}/results/main_outputs/clustered" }, mode: 'copy'

    input:
    tuple val(primers), path(seqtab), path(taxtab), path(samdf), path(raw), path(summary), path(ps), path(clusters)
    val(merge_clusters)

    output:
    path("*.csv"),                                                          emit: csvs
    tuple val(primers), path("ps_clustered_*.rds"),                         emit: ps 
    path("asvs_clustered*.fasta"),                                          emit: asv_fasta

    script:
    """
    
    Rscript --vanilla "${moduleDir}/phyloseq_clustered.R" \
        --process_name="phyloseq_clustered" \
        --cpus="$task.cpus" \
        --rdata="$params.rdata" \
        --primers="$primers" \
        --seqtab_file="$seqtab" \
        --taxtab_file="$taxtab" \
        --samdf_file="$samdf" \
        --raw_file="$raw" \
        --summary_file="$summary" \
        --ps_file="$ps" \
        --clusters_file="$clusters" \
        --merge_clusters="$merge_clusters" 
        
    """


}