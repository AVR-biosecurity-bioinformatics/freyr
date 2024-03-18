#!/usr/bin/env Rscript

#### TODO: Pull stats from top of index_switching.pdf and print to a final run report
### also flag index combos that have higher than average + (some stat)

#### TODO: use index_switch_calc.txt in jack_notes to run same process but in bash, which should be much faster

# run flow cell QC
step_seq_qc(flowcell_id)

# run index switching calculation
step_switching_calc(flowcell_id)

# copy output files to log folder
file.copy(paste0(flowcell_id,"_flowcell_qc.pdf"),paste0(projectDir,"/output/logs/",flowcell_id))
file.copy(paste0(flowcell_id,"_index_switching.pdf"),paste0(projectDir,"/output/logs/",flowcell_id))
