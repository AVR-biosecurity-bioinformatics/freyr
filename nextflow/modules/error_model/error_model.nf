process ERROR_MODEL {
    tag "$primers; $read_group"
    label "high"
    label "env_r_dada2"
    publishDir { "${params.outdir}/modules/error_model" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()
    publishDir { "${params.outdir}/results/qc/dada2_error_model" }, mode: 'copy'

    input:
    tuple val(direction), val(primers), val(read_group), path(reads)

    output:   
    tuple val(direction), val(primers), val(read_group), path("*_errormodel{F,R,S}.rds"),   emit: errormodel
    path("*_errormodel.pdf"),                                           emit: plot

    script:
    """

    Rscript --vanilla "${moduleDir}/error_model.R" \
        --process_name "error_model" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --reads "$reads" \
        --primers "$primers" \
        --read_group "$read_group" \
        --direction "$direction"

    """
    
}