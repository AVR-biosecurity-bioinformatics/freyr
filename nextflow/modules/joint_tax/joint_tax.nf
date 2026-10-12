process JOINT_TAX {
    tag "$primers; $read_group"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"
    publishDir "${launchDir}/output/modules/joint_tax", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(primers), val(read_group), path(tax, name: "idtaxa*.csv"), path(blast, name: "blast*.csv")
    
    output:
    tuple val(primers), val(read_group), path("*_joint.csv"), emit: joint

    script:
    """
    
    joint_tax.R \
        --process_name "joint_tax" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --read_group "$read_group" \
        --idtaxa_list "$tax" \
        --blast_list "$blast" 
        
    """

}