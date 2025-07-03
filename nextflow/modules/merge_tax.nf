process MERGE_TAX {
    def process_name = "merge_tax"
    tag "$primers"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), path(taxtabs)

    output:
    tuple val(primers), path("*_merged_tax.csv"), emit: merged_tax

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy'

    // when: 

    script:
    """
    
    ${process_name}.R \
        --process_name "$process_name" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --taxtabs "$taxtabs" 
        
    """
    
}