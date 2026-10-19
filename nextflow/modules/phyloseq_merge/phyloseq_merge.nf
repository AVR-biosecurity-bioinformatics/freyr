process PHYLOSEQ_MERGE {
    tag "Whole dataset"
    label "high"
    label "env_r_phyloseq"
    publishDir { "${params.outdir}/modules/phyloseq_merge" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()
    publishDir { "${params.outdir}/results/main_outputs/merged" }, pattern: '*{_filtered,_unfiltered}.{csv,rds,fasta}', mode: 'copy'

    input:
    tuple val(dataset),
        path(seqtab_files),
        path(taxtab_files),
        path(samdf_files),
        path(fasta_files)

    output:
    path("*{_filtered,_unfiltered}.{csv,rds,fasta}")
    path("*_readsout.csv"), emit: read_tracking

    script:
    """
        
    Rscript --vanilla "${moduleDir}/phyloseq_merge.R" \
        --process_name="phyloseq_merge" \
        --cpus="$task.cpus" \
        --rdata="$params.rdata" \
        --dataset="$dataset" \
        --seqtab_files="${seqtab_files.join(' ')}" \
        --taxtab_files="${taxtab_files.join(' ')}" \
        --samdf_files="${samdf_files.join(' ')}" \
        --fasta_files="${fasta_files.join(' ')}"

    """
    
}