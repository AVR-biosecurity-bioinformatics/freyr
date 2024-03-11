# TODOs for the pipeline

### Allow for loci parameters to be input as a spreadsheet of some kind


> Make sure `jackscanlan/piperline` container has `bash:3.x` and `ps` installed, and make `/bin/bash` the container `ENTRYPOINT`. 

TODO: troubleshoot SEQ_QC; some issue with switching calcs

TODO: 



TODO: Implement `nf-validation` (ala [here](https://nextflow-io.github.io/nf-validation/latest/)) to handle samplesheets etc. 
- also read [plugin documentation](https://www.nextflow.io/docs/latest/plugins.html)


### LATER

TODO: Add reasonable process labels to each module. 

TODO: Add version-tracking for software used in each module ala. ampliseq, with a version channel. 

TODO: Find a method to get the output of a failed process, from within a sourced R script, to display as a log file or otherwise in the stderr or stdout on terminal
