process TAX_SUMMARY_MERGE {
    tag "Whole dataset"
    label "small"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/tax_summary_merge", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/taxonomic_assignment", mode: 'copy'

    input:
    path(tax_summary_list)

    output:
    path("taxonomic_assignment_summary_combined.csv")

    script:
    """
    
    ${moduleDir}/tax_summary_merge.R \
        --process_name "tax_summary_merge" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --tax_summary_list "$tax_summary_list"
        
    """

}