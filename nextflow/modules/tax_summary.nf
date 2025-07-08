process TAX_SUMMARY {
    def process_name = "tax_summary"
    tag "$primers; $read_group"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), val(read_group), path(fasta), path(tax, name: "idtaxa*.csv"), path(ids, name: "idtaxa*.rds"), path(joint), val(process_params)

    output:
    tuple val(primers), val(read_group), path("taxonomic_assignment_summary_*.rds"), emit: rds
    tuple val(primers), val(read_group), path("taxonomic_assignment_summary_*.csv"), emit: csv

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    publishDir "${launchDir}/output/results/taxonomic_assignment", pattern: 'taxonomic_assignment_summary_*.csv', mode: 'copy'

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
        --tax_list "$tax" \
        --ids_list "$ids" \
        --joint_file "$joint" \
        --locus "$process_params.locus" \
        --idtaxa_db "$process_params.idtaxa_db" \
        --ref_fasta "$process_params.ref_fasta"
        
    """

}