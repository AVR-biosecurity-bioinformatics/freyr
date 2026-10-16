process TRAIN_IDTAXA {
    tag "Whole pipeline"
    label "high"
    label "env_r_phyloseq"
    publishDir { "${params.outdir}/modules/train_idtaxa" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()
    publishDir { "${params.outdir}/results/taxonomic_assignment" }, mode: 'copy'

    input:
    tuple val(primers), path(ref_fasta)
    val(max_group_size)
    val(max_iterations)
    val(allow_group_removal)
    val(drop_malformed)

    output:   
    tuple val(primers), path("*_idtaxa_db.rds")             , emit: model

    script:
    """
    
    Rscript --vanilla "${moduleDir}/train_idtaxa.R" \
        --process_name "train_idtaxa" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --ref_fasta "$ref_fasta" \
        --max_group_size "$max_group_size" \
        --max_iterations "$max_iterations" \
        --allow_group_removal "$allow_group_removal" \
        --drop_malformed "$drop_malformed"
        
    """

}