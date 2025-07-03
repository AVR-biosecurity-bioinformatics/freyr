process TAX_SUMMARY_MERGE {
    def process_name = "tax_summary_merge"
    tag "Whole dataset"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(tax_summary_list)

    output:
    path("taxonomic_assignment_summary.csv")

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy'

    // when: 

    script:
    """
    
    ${process_name}.R \
        --process_name "$process_name" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --tax_summary_list "$tax_summary_list"
        
    """

}