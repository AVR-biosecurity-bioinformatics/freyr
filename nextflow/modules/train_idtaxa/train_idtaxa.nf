process TRAIN_IDTAXA {
    tag "Whole pipeline"
    label "high"
    label "env_r_phyloseq"
    publishDir "${params.outdir}/modules/process_name", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${params.outdir}/results/taxonomic_assignment", mode: 'copy'

    input:
    tuple val(primers), path(ref_fasta)
    val(fasta_type)

    output:   
    tuple val(primers), path("*_idtaxa_db.rds")             , emit: model

    script:
    """
    
    Rscript --vanilla "${moduleDir}/process_name.R" \
        --process_name "train_idtaxa" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --ref_fasta "$ref_fasta" \
        --fasta_type "$fasta_type" 
        
    """

}