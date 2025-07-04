process PHYLOSEQ_UNFILTERED {
    def process_name = "phyloseq_unfiltered"
    tag "$primers"
    label "phyloseq"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), path(taxtab), path(seqtab), path(filters), path(fasta), path(samplesheet_split), path(sample_metadata), val(process_params)

    output:
    tuple val(primers), path("seqtab_unfiltered_*.csv"), path("taxtab_unfiltered_*.csv"), path("samdf_unfiltered_*.csv"), path("raw_unfiltered_*.csv"), path("summary_unfiltered_*.csv"), emit: csvs
    tuple val(primers), path("ps_unfiltered_*.rds"), path("filters_*.csv"),                        emit: ps 
    path("asvs_unfiltered_*.fasta"),                                                               emit: asv_fasta

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    // when: 

    script:
    """
    
    ${process_name}.R \
        --process_name "$process_name" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --taxtab_file "$taxtab" \
        --seqtab_file "$seqtab" \
        --filters_file "$filters" \
        --fasta_file "$fasta" \
        --samplesheet_split_file "$samplesheet_split" \
        --sample_metadata_file "$sample_metadata" \
        --cluster_threshold "$process_params.cluster_threshold" 
        
    """

}