process ASSIGNMENT_PLOT {
    tag "$primers; $read_group"
    label "small"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/assignment_plot", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/taxonomic_assignment", pattern: 'taxonomic_assignment_summary_*.pdf', mode: 'copy'

    input:
    tuple val(primers), val(read_group), path(fasta), path(blast, name: "blast*.rds"), path(tax), val(process_params)

    output:
    path("taxonomic_assignment_summary_*.pdf"),                 emit: plot
    tuple val(primers), val(read_group), path("*_joint.rds"),   emit: joint

    script:
    """
    
    ${moduleDir}/assignment_plot.R \
        --process_name "assignment_plot" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --read_group "$read_group" \
        --fasta "$fasta" \
        --tax "$tax" \
        --blast_list "$blast" \
        --idtaxa_db "$process_params.idtaxa_db" \
        --ref_fasta "$process_params.ref_fasta"
        
    """

}