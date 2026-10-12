process PRIORS {
    tag "$primers; $read_group"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"
    publishDir "${launchDir}/output/modules/priors", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(direction), val(primers), val(read_group), path(priors)
    val(prior_min_samples)

    output:   
    tuple val(direction), val(primers), val(read_group), path("*_priors{F,R,S}.rds"),              emit: priors


    script:
    """
    
    priors.R \
        --process_name "priors" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --direction "$direction" \
        --read_group "$read_group" \
        --priors "$priors" \
        --prior_min_samples "$prior_min_samples" 
        
    """
    
}