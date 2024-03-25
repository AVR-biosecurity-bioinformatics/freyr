process SPLIT_LOCI {
    def module_name = "split_loci"
    // tag: 
    // label:  

    input:
    // input read pairs, primer seqs and locus name (to use in output file names)
    tuple val(meta), path(reads)

    output:   
    // path "output_file.txt"

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/bash

    ### run module code
    bash ${module_name}.sh
    
    """

}