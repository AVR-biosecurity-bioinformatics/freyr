#!/usr/bin/env Rscript


## check and define variables


### run R code

# to run classic pipeline code, need:
# - merged seqtab (all sequences across all flowcells and loci)
# - merged taxtable (all sequences across all flowcells and loci, identified)
# - samplesheet (ie. samdf)

# should input all the files independently, merge within this code, then run filtering etc.
# tax tables (output of MERGE_TAX) are already per locus (merged across flowcells)
# sequence tables (output of FILTER_SEQTAB) are per locus x flowcell



## runs step_phyloseq()

## runs step_output_summary() and step_output_ps() on unfiltered data

## runs rareplot() and saves plot

## runs step_filter_phyloseq() on each locus independently

## runs merge_phyloseq_new() from functions.R

## runs step_output_summary() and step_output_ps() on filtered data





# stop(" *** stopped manually *** ") ##########################################
