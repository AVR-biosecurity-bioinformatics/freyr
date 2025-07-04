process ACCUMULATION_CURVE {
    def process_name = "accumulation_curve"
    tag "$primers"
    label "phyloseq"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), path(ps_file), path(filters_tibble), val(process_params)

    output:
    path("accumulation_curve_*.pdf"),                                          emit: pdf

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
        --ps_file "$ps_file" \
        --min_sample_reads "$process_params.min_sample_reads"
        
    """

}