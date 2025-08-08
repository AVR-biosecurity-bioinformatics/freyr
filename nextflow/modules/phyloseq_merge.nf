process PHYLOSEQ_MERGE {
    def process_name = "phyloseq_merge"
    tag "Whole dataset"
    label "phyloseq"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(ps_unfiltered)
    path(ps_filtered)
    path(unfiltered_fastas, name: "unfiltered_*.fasta")
    path(filtered_fastas, name: "filtered_*.fasta")

    output:
    path("*{_filtered,_unfiltered}.{csv,rds,fasta}")
    path("*_readsout.csv"), emit: read_tracking

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    publishDir "${launchDir}/output/results/main_outputs/merged", pattern: '*{_filtered,_unfiltered}.{csv,rds,fasta}', mode: 'copy'

    // when: 

    script:
    """
        
    ${process_name}.R \
        --process_name "$process_name" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --ps_unfiltered "$ps_unfiltered" \
        --ps_filtered "$ps_filtered" \
        --unfiltered_fastas "$unfiltered_fastas" \
        --filtered_fastas "$filtered_fastas" 
        
        
    """
    
}