#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

reads                       <- args$reads
primers                     <- args$primers
read_group                  <- args$read_group
direction                   <- args$direction
threads                     <- args$cpus

sys.source(paste0(args$projectDir,"/bin/functions.R"), envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "dada2",
    "ggplot2",
    "magrittr",
    "readr",
    "stringr",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### run R code 
if (!direction %in% c("forward","reverse","single")) { 
    stop(" Input reads direction needs to be 'forward', 'reverse' or 'single'! ")
}

if (direction == "forward") {
    direction_short <- "F"
} else if (direction == "reverse") {
    direction_short <- "R"
} else if (direction == "single") {
    direction_short <- "S"
} else {
    stop(" Input reads direction needs to be 'forward', 'reverse' or 'single'! ")
}

reads_list <- # convert input reads list from Groovy format to R format
    stringr::str_extract_all(
        reads, 
        pattern = "\\S+?\\.fastq\\.gz|\\S+?\\.fastq|\\S+?\\.fq\\.gz|\\S+?\\.fq" 
    ) %>% 
    unlist()

if(direction == "forward"){
    message(paste0("Modelling forward read error rates for primers: ", primers, " and flowcell: ", read_group))
} else if (direction == "reverse"){
    message(paste0("Modelling reverse read error rates for primers: ", primers, " and flowcell: ", read_group))
} else if (direction == "single"){
    message(paste0("Modelling single read error rates for primers: ", primers, " and flowcell: ", read_group))
} else {
    stop ("Read direction must be 'forward', 'reverse' or 'single'!")
}

## Learn error rates 
set.seed(1); err <- dada2::learnErrors(
    reads_list, 
    multithread = as.numeric(threads), 
    nbases = 1e8,
    randomize = FALSE, 
    # qualityType = "FastqQuality",
    qualityType = "Auto",
    verbose=TRUE
)
    

## save output as .rds file
saveRDS(err, paste0(read_group,"_",primers,"_errormodel",direction_short,".rds"))

## write out errors for diagnostics
readr::write_csv(as.data.frame(err$trans), paste0(read_group,"_",primers,"_err",direction_short,"_observed_transitions.csv"))
readr::write_csv(as.data.frame(err$err_out), paste0(read_group,"_",primers,"_err",direction_short,"_inferred_errors.csv"))

## output error plots to see how well the algorithm modelled the errors in the different runs
p1 <- dada2::plotErrors(err, nominalQ = TRUE) + ggtitle(paste0(primers, " ", read_group," ",direction," Reads"))
pdf(paste0(read_group,"_", primers, "_", direction_short,"_errormodel.pdf"), width = 11, height = 8 , paper="a4r")
    plot(p1)
try(dev.off(), silent=TRUE)

# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})