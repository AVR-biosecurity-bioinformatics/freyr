process PARSE_INPUTS {
    def process_name = "parse_inputs"
    def projDir = workflow.projectDir // redeclaring these two variables stops cache invalidation
    def launDir = workflow.launchDir

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
    val(pp_params)

    output: 
    path("sample_metadata.csv"),                    emit: sample_metadata
    path("samplesheet_unsplit.csv"),                emit: samplesheet_unsplit
    path("samplesheet_split.csv"),                  emit: samplesheet_split
    path("primer_params_parsed.csv"),               emit: primer_params_parsed

    publishDir "${launchDir}/output/modules/${process_name}",  mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    publishDir "${launchDir}/output/results/parsed_inputs", mode: 'copy'


    // when: 

    script:
    """

    ${process_name}.R \
        --process_name "$process_name" \
        --projectDir "$projDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --launchDir "$launDir" \
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