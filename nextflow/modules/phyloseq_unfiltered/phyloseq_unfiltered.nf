process PHYLOSEQ_UNFILTERED {
    tag "$primers"
    label "env_r_phyloseq"
    publishDir { "${params.outdir}/modules/phyloseq_unfiltered" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()
    publishDir { "${params.outdir}/results/main_outputs/unfiltered" }, mode: 'copy'

    input:
    tuple val(primers), path(taxtab), path(seqtab), path(filters), path(fasta), path(samplesheet_split), path(sample_metadata), val(process_params)

    output:
    tuple val(primers),
      path("seqtab_unfiltered_*.csv"),
      path("taxtab_unfiltered_*.csv"),
      path("samdf_unfiltered_*.csv"),
      path("asvs_unfiltered_*.fasta"),
      path("filters_*.csv"),
      emit: components    

    tuple val(primers), path("seqtab_unfiltered_*.csv"), path("taxtab_unfiltered_*.csv"), path("samdf_unfiltered_*.csv"), path("raw_unfiltered_*.csv"), path("summary_unfiltered_*.csv"), emit: csvs
    tuple val(primers), path("ps_unfiltered_*.rds"), path("filters_*.csv"),                        emit: ps 
    path("asvs_unfiltered_*.fasta"),                                                               emit: asv_fasta

    script:
    """
    
    Rscript --vanilla "${moduleDir}/phyloseq_unfiltered.R" \
        --process_name="phyloseq_unfiltered" \
        --cpus="$task.cpus" \
        --rdata="$params.rdata" \
        --primers="$primers" \
        --taxtab_file="$taxtab" \
        --seqtab_file="$seqtab" \
        --filters_file="$filters" \
        --fasta_file="$fasta" \
        --samplesheet_split_file="$samplesheet_split" \
        --sample_metadata_file="$sample_metadata" \
        --cluster_threshold="$process_params.cluster_threshold" 
        
    """

}