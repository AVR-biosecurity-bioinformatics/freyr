process TRAIN_IDTAXA {
    tag "Whole pipeline"
    label "high"
    container "jackscanlan/piperline-multi:0.0.1"
    publishDir "${launchDir}/output/modules/process_name", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/taxonomic_assignment", mode: 'copy'

    input:
    tuple val(primers), path(ref_fasta)
    val(fasta_type)

    output:   
    tuple val(primers), path("*_idtaxa_db.rds")             , emit: model

    script:
    """
    
    process_name.R \
        --process_name "train_idtaxa" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --ref_fasta "$ref_fasta" \
        --fasta_type "$fasta_type" 
        
    """

}