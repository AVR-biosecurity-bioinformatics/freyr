process MERGE_TAX {
    tag "$primers"
    label "small"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/merge_tax", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(primers), path(taxtabs)

    output:
    tuple val(primers), path("*_merged_tax.csv"), emit: merged_tax

    script:
    """
    
    merge_tax.R \
        --process_name "merge_tax" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --taxtabs "$taxtabs" 
        
    """
    
}