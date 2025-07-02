process DENOISE {
    def process_name = "denoise"
    tag "$sample_primers"
    label "medium"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(direction), val(primers), val(read_group), val(sample), val(sample_primers), path(reads), path(errormodel), path(priors)
    val(n_pass)
    val(dada_band_size)
    val(dada_homopolymer)

    output:
    tuple val(direction), val(primers), val(read_group), val(sample), val(sample_primers), path(reads), path("*_dada{1,2}{F,R,S}.rds"), emit: seq

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy'

    // when: 

    script:
    // def priors_path = priors.name != "NO_FILE" ? priors.name : "NO_FILE"
    """

    ${process_name}.R \
        --process_name "$process_name" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --reads "$reads" \
        --sample_primers "$sample_primers" \
        --primers "$primers" \
        --read_group "$read_group" \
        --direction "$direction" \
        --errormodel "$errormodel" \
        --n_pass "$n_pass" \
        --priors "$priors" \
        --dada_band_size "$dada_band_size" \
        --dada_homopolymer "$dada_homopolymer"

    """
    
}