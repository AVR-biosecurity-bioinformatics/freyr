#!/usr/bin/env Rscript

# check variables defined


### run R code (from step_errormodel)
if (direction %in% c("forward","reverse")) { 
    message(sprintf("Read direction is %s!",direction))
} else (
    stop(" Input reads direction needs to be set! ")
)

reads_list <- # convert input reads list from Groovy format to R format
    stringr::str_extract_all(
        reads, 
        pattern = "\\S+?\\.fastq\\.gz|\\S+?\\.fastq|\\S+?\\.fq\\.gz|\\S+?\\.fq" 
        ) %>% 
    unlist()

# input_dir <- normalizePath(input_dir)
# output <- normalizePath(output)
# qc_dir <- normalizePath(qc_dir)

if(direction == "forward"){
    # filts <- list.files(input_dir, pattern= "*R1_001.*", full.names = TRUE)
    message(paste0("Modelling forward read error rates for primers: ", pcr_primers, " and flowcell: ", fcid))
} else if (direction == "reverse"){
    # filts <- list.files(input_dir, pattern="R2_001.*", full.names = TRUE)
    message(paste0("Modelling reverse read error rates for primers: ", pcr_primers, " and flowcell: ", fcid))
} else {
    stop ("Read direction must be 'forward' or 'reverse'!")
}
# # Subset fastqs to just the relevent pcr primers
# filts <- filts[str_detect(filts,paste0(pcr_primers, "(-|_|$)"))]
# message(paste0(length(filts), " fastq files to process for primers: ", pcr_primers, " and flowcell: ", fcid))

# Learn error rates from a subset of the samples and reads (rather than running self-consist with full dataset)
err <- dada2::learnErrors(
        reads_list, 
        multithread = FALSE, 
        nbases = 1e8,
        randomize = FALSE, 
        qualityType = "FastqQuality",
        verbose=TRUE
        )
    
### TODO: increase nbases or use default (1e8)

# print(err)

# stop(" *** stopped manually *** ") ##########################################


#write out errors for diagnostics
write_csv(as.data.frame(err$trans), paste0(fcid,"_",pcr_primers,"_err",direction,"_observed_transitions.csv"))
write_csv(as.data.frame(err$err_out), paste0(fcid,"_",pcr_primers,"_err",direction,"_inferred_errors.csv"))


# ##output error plots to see how well the algorithm modelled the errors in the different runs
# p1 <- dada2::plotErrors(err, nominalQ = TRUE) + ggtitle(paste0(pcr_primers, " ", fcid, " Forward Reads"))
# pdf(paste0(qc_dir,"/",fcid,"_", pcr_primers, "_", read,"_errormodel.pdf"), width = 11, height = 8 , paper="a4r")
#     plot(p1)
# try(dev.off(), silent=TRUE)