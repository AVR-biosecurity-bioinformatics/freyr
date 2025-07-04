process ASSIGNMENT_PLOT {
    def process_name = "assignment_plot"
    tag "$primers; $read_group"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), val(read_group), path(fasta), path(blast, name: "blast*.rds"), path(tax), val(process_params)

    output:
    path("*_taxonomic_assignment_summary.pdf"),                 emit: plot
    tuple val(primers), val(read_group), path("*_joint.rds"),   emit: joint

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    publishDir "${launchDir}/output/results/taxonomic_assignment", pattern: '*_taxonomic_assignment_summary.pdf', mode: 'copy'


    // when: 

    script:
    """
    
    ${process_name}.R \
        --process_name "$process_name" \
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