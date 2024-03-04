> Add custom `.Rprofile` to GitHub repo (in `/jack_notes` dir?) that can be sourced at run time. This contains the location of the library that should have all the packages needed for analysis. 


Need to source copy of `.Rprofile` that contains new default library location and option to use `pak` in `renv`. 

    > source("./jack_notes/.Rprofile")

Add above code to the "R sourcing" code block at the start of every R script (could I wrap this in an R script?)


TODO: check for all files previously checked in the `parameters_setup.R` script, which no longer checks because Nextflow breaks the relative paths