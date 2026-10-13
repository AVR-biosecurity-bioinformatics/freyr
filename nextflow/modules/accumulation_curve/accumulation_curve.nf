process ACCUMULATION_CURVE {
    tag "$primers"
    label "phyloseq"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/accumulation_curve", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/main_outputs/unfiltered", mode: 'copy'

    input:
    tuple val(primers), path(ps_file), path(filters_tibble), val(process_params)

    output:
    path("accumulation_curve_*.pdf"),                                          emit: pdf

    script:
    """
    ${moduleDir}/accumulation_curve.R \
        --process_name "accumulation_curve" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --ps_file "$ps_file" \
        --min_sample_reads "$process_params.min_sample_reads"
        
    """

}