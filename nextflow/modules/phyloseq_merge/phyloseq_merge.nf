process PHYLOSEQ_MERGE {
    tag "Whole dataset"
    label "phyloseq"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/phyloseq_merge", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/main_outputs/merged", pattern: '*{_filtered,_unfiltered}.{csv,rds,fasta}', mode: 'copy'

    input:
    path(ps_unfiltered)
    path(ps_filtered)
    path(unfiltered_fastas, name: "unfiltered_*.fasta")
    path(filtered_fastas, name: "filtered_*.fasta")

    output:
    path("*{_filtered,_unfiltered}.{csv,rds,fasta}")
    path("*_readsout.csv"), emit: read_tracking

    script:
    """
        
    phyloseq_merge.R \
        --process_name "phyloseq_merge" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --ps_unfiltered "$ps_unfiltered" \
        --ps_filtered "$ps_filtered" \
        --unfiltered_fastas "$unfiltered_fastas" \
        --filtered_fastas "$filtered_fastas" 
        
        
    """
    
}