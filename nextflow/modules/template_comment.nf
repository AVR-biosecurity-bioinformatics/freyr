process TEMPLATE {
    // Top is used for process directives like tag, label, time etc. 
    tag //$some_var
    /* Tags are used to associate each running process with a custom label, perhaps based on the file being used as input. 
    See here: https://www.nextflow.io/docs/latest/process.html#tag */
    
    // label 'process_low'
    /* Labels can be used to group processes based on certain qualities.
    This is commonly used to group by resource allocation and nf-core has standard ones. */

    input:
    path 


    output:
    /* using ", emit: name" at the end of an output line can be used to name this output later.
    For example, in ampliseq, RENAME_RAW_DATA_FILES uses ", emit: fastq" for the renamed files.
    In the main workflow, when RENAME_RAW_DATA_FILES is called, the fastq output can be used by another process
    as "RENAME_RAW_DATA_FILES.out.fastq". */

    /* The ampliseq modules/processes also output the version of the software used with additional
    bit of code. I could do this later but not worth trying for first implementation. 
    */

  


    when: 

    /* The when block allows for conditions to be set on the running of the process, for example
    if the input file name or other characteristic matches a pattern. 
    HOWEVER: this is usually better done in the workflow block (ie. in the code where the process
    is called) so the module is more portable. But there could be good uses for it.
    The ampliseq pipeline uses "task.ext.when == null || task.ext.when" in its when blocks, but
    I'm not sure what they mean yet. 
    */


    script:
    /* This is an R script, so needs to start with "#!/usr/bin/env Rscript".
    If I want to source the "functions.R" file every time, I can keep the path to that file as 
    a variable, perhaps "r_functions", and define that in the main pipeline. */ 
    
    /* This is based on the "prepare_inputs.R" script in the R pipeline. Note that "$" in R needs to be escaped like "\$" to avoid 
    it being parsed as a Nextflow variable. */
    
    """
    #!/usr/bin/env Rscript

    # source functions, themes and load packages
    source("${projectDir}/bin/functions.R")
    source("${projectDir}/bin/themes.R")
    source("${projectDir}/bin/_targets_packages.R")

    # run module code
    source("${projectDir}/bin/template.R")

    """

}