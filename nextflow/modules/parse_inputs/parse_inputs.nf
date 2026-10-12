process PARSE_INPUTS {
    publishDir "${launchDir}/output/modules/parse_inputs",  mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/parsed_inputs", mode: 'copy'

    tag "Whole dataset"
    label "small"
    conda "${moduleDir}/environment.yml"

    input:
    path(samplesheet)
    path(primer_params)
    val(pp_type)
    val(seq_type)
    val(paired)
    val(subsample)
    val(extension)
    val(pp_params)

    output: 
    path("sample_metadata.csv"),                    emit: sample_metadata
    path("samplesheet_unsplit.csv"),                emit: samplesheet_unsplit
    path("samplesheet_split.csv"),                  emit: samplesheet_split
    path("primer_params_parsed.csv"),               emit: primer_params_parsed

    script:
    """

    parse_inputs.R \
        --process_name "parse_inputs" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --launchDir "$launchDir" \
        --samplesheet "$samplesheet" \
        --primer_params "$primer_params" \
        --pp_type "$pp_type" \
        --seq_type "$seq_type" \
        --paired "$paired" \
        --subsample "$subsample" \
        --extension "$extension" \
        --pp_params "$pp_params"

    """
}