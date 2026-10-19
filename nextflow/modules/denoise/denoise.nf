process DENOISE {
    tag "$sample_primers"
    label "medium"
    label "env_r_dada2"
    publishDir { "${params.outdir}/modules/denoise" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()

    input:
    tuple val(direction), val(primers), val(read_group), val(sample), val(sample_primers), path(reads), path(errormodel), path(priors)
    val(n_pass)
    val(dada_band_size)
    val(dada_homopolymer)

    output:
    tuple val(direction), val(primers), val(read_group), val(sample), val(sample_primers), path(reads), path("*_dada{1,2}{F,R,S}.rds"), emit: seq

    script:
    // def priors_path = priors.name != "NO_FILE" ? priors.name : "NO_FILE"
    """

    Rscript --vanilla "${moduleDir}/denoise.R" \
        --process_name="denoise" \
        --cpus="$task.cpus" \
        --rdata="$params.rdata" \
        --reads="$reads" \
        --sample_primers="$sample_primers" \
        --primers="$primers" \
        --read_group="$read_group" \
        --direction="$direction" \
        --errormodel="$errormodel" \
        --n_pass="$n_pass" \
        --priors="$priors" \
        --dada_band_size="$dada_band_size" \
        --dada_homopolymer="$dada_homopolymer"

    """
    
}