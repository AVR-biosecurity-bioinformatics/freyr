process TAX_BLAST {
    def process_name = "tax_blast"
    tag "$primers; $read_group"
    label "medium"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), val(read_group), path(fasta), val(process_params)

    output:
    tuple val(primers), val(read_group), path("*_blast.csv"),                     emit: blast
    tuple val(primers), val(read_group), path("*_blast_spp_low.rds"),             emit: blast_assignment
    tuple val(primers), val(read_group), path("*_n_ranks.txt"),                   emit: n_ranks

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
        --read_group "$read_group" \
        --fasta "$fasta" \
        --ref_fasta "$process_params.ref_fasta" \
        --blast_min_identity "$process_params.blast_min_identity" \
        --blast_min_coverage "$process_params.blast_min_coverage" \
        --run_blast "$process_params.run_blast"
        
    """

}