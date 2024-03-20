##### this loads the listed R scripts and imports Nextflow variables from params

source(file.path(projectDir, "jack_notes/.Rprofile"))
source(file.path(projectDir, "bin/functions.R"))
source(file.path(projectDir, "bin/themes.R"))
source(file.path(projectDir, "bin/_targets_packages.R"))

### parse Nextflow params dictionary (aka. "params" in Nextflow) directly into R variables
if (!exists("params_dict")) {stop("'params_dict' not found; check it is defined in .nf")}

params_list <- params_dict %>% # convert Groovy list format into R nested list
    stringr::str_remove_all("\\[|\\]") %>% 
    stringr::str_split_1(", ") %>% 
    stringr::str_split(":")

for (i in 1:length(params_list)) { # loop through components of list
    row <- params_list[[i]] # select 'i'th row
    if (row[2] != "null") { # assign parameter to variable if value is not "null" (empty)
        assign(
            paste0("params.",row[1]),
            paste0(row[2])
        )
    }
}
### TODO: make a printable df of the new variables and their values?
