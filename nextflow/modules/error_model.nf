process ERROR_MODEL {
    def process_name = "error_model"
    tag "$primers; $read_group"
    label "high"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(direction), val(primers), val(read_group), path(reads)

    output:   
    tuple val(direction), val(primers), val(read_group), path("*_errormodel{F,R,S}.rds"),   emit: errormodel
    path("*_errormodel.pdf"),                                           emit: plot

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    // when: 

    script:
    """

    ${process_name}.R \
        --process_name "$process_name" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --reads "$reads" \
        --primers "$primers" \
        --read_group "$read_group" \
        --direction "$direction"

    """
    
}