process TAX_SUMMARY_MERGE {
    tag "Whole dataset"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"
    publishDir "${launchDir}/output/modules/tax_summary_merge", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/taxonomic_assignment", mode: 'copy'

    input:
    path(tax_summary_list)

    output:
    path("taxonomic_assignment_summary_combined.csv")

    script:
    """
    
    tax_summary_merge.R \
        --process_name "tax_summary_merge" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --tax_summary_list "$tax_summary_list"
        
    """

}