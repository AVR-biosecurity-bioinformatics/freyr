process JOINT_TAX {
    def process_name = "joint_tax"
    tag "$primers; $read_group"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), val(read_group), path(tax, name: "idtaxa*.csv"), path(blast, name: "blast*.csv")
    
    output:
    tuple val(primers), val(read_group), path("*_joint.csv"), emit: joint

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
        --idtaxa_list "$tax" \
        --blast_list "$blast" 
        
    """

}