process PRIORS {
    def process_name = "priors"
    tag "$primers; $read_group"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(direction), val(primers), val(read_group), path(priors)

    output:   
    tuple val(direction), val(primers), val(read_group), path("*_priors{F,R,S}.rds"),              emit: priors

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
        --direction "$direction" \
        --read_group "$read_group" \
        --priors "$priors" 
        
    """
    
}