process PARSE_INPUTS {
    def process_name = "parse_inputs"
    tag "Whole dataset"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(samplesheet)
    path(primer_params)
    val(pp_type)
    val(seq_type)
    val(paired)
    val(subsample)
    val(extension)

    output: 
    path("sample_metadata.csv"),                    emit: sample_metadata
    path("samplesheet_unsplit.csv"),                emit: samplesheet_unsplit
    path("samplesheet_split.csv"),                  emit: samplesheet_split
    path("primer_params_parsed.csv"),               emit: primer_params_parsed

    publishDir "${launchDir}/output/modules/${process_name}",  mode: 'copy'

    // when: 

    script:
    """

    ${process_name}.R \
        --process_name "$process_name" \
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
        --params "$params"

    """
}