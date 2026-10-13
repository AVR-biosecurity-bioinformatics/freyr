process TRAIN_IDTAXA {
    tag "Whole pipeline"
    label "high"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/process_name", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/taxonomic_assignment", mode: 'copy'

    input:
    tuple val(primers), path(ref_fasta)
    val(fasta_type)

    output:   
    tuple val(primers), path("*_idtaxa_db.rds")             , emit: model

    script:
    """
    
    ${moduleDir}/process_name.R \
        --process_name "train_idtaxa" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --ref_fasta "$ref_fasta" \
        --fasta_type "$fasta_type" 
        
    """

}