process TAX_IDTAXA {
    def process_name = "tax_idtaxa"
    tag "$primers; $read_group"
    label "medium"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), val(read_group), path(fasta), val(process_params)

    output:
    tuple val(primers), val(read_group), path("*_idtaxa_tax.csv"), emit: tax
    tuple val(primers), val(read_group), path("*_idtaxa_ids.rds"), emit: ids

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
        --read_group "$read_group" \
        --idtaxa_confidence "$process_params.idtaxa_confidence" \
        --idtaxa_db "$process_params.idtaxa_db" \
        --fasta "$fasta" 
        
    """

}