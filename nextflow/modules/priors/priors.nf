process PRIORS {
    tag "$primers; $read_group"
    label "small"
    label "env_r_dada2"
    publishDir { "${params.outdir}/modules/priors" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()

    input:
    tuple val(direction), val(primers), val(read_group), path(priors)
    val(prior_min_samples)

    output:   
    tuple val(direction), val(primers), val(read_group), path("*_priors{F,R,S}.rds"),              emit: priors


    script:
    """
    
    Rscript --vanilla "${moduleDir}/priors.R" \
        --process_name="priors" \
        --cpus="$task.cpus" \
        --rdata="$params.rdata" \
        --primers="$primers" \
        --direction="$direction" \
        --read_group="$read_group" \
        --priors="$priors" \
        --prior_min_samples="$prior_min_samples" 
        
    """
    
}