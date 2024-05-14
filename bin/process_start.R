### load listed R scripts 
if (!exists("projectDir")) { projectDir <- "." } # if not running in Nextflow, use current directory

source(file.path(projectDir, "jack_notes/.Rprofile"))
source(file.path(projectDir, "bin/functions.R"))
source(file.path(projectDir, "bin/themes.R"))
# source(file.path(projectDir, "bin/_targets_packages.R")) ### TODO: remove this once all R scripts directly load required packages

### parse Nextflow params dictionary (aka. "params" in Nextflow) directly into R variables
## NOTE: this works only as long as no parameter value contains ", " (but it can handle spaces)
if (!exists("params_dict")) {
    stop("'params_dict' not found; check it is defined in process .nf file!")
}

invisible(library(stringr)) # load stringr for parsing

params_list <- params_dict %>% # convert Groovy list format into R nested list
                stringr::str_remove_all("\\[|\\]") %>% 
                stringr::str_split_1(", ") %>% 
                stringr::str_split(":")

params_df = data.frame()

for (i in 1:length(params_list)) {
    row <- params_list[[i]]
    if (row[2] != "null") {
        assign(
            paste0("params.",row[1]),
            paste0(row[2])
        )
    }
output <- c("parameter" = paste0("params.",row[1]), "value" = row[2])
params_df <- rbind(params_df, output)
}

colnames(params_df) <- c("parameter","value")
# print(params_df) # access params if you need them
