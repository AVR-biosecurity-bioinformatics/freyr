process STOP {
    
    // this process is used to halt the pipeline early, even if the last process finished successfully

    script:
    """
    #!/usr/bin/bash

    echo "***** Pipeline has been manually stopped with the STOP module *****"
    
    exit 1
    
    """
}