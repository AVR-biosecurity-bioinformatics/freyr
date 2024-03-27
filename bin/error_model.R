#!/usr/bin/env Rscript

# check variables defined

print(reads)
print(direction)
stringr::str_split(reads, pattern = " ") %>% as.data_frame() %>%
    write_csv(file = "fwd.csv", x = .)

stop(" *** stopped manually *** ") ##########################################

### run R code (from step_errormodel)
if (direction == "forward") { input_reads <- "fwd_reads"

input_dir <- normalizePath(input_dir)
output <- normalizePath(output)
qc_dir <- normalizePath(qc_dir)
if(read == "F"){
    filts <- list.files(input_dir, pattern= "*R1_001.*", full.names = TRUE)
    message(paste0("Modelling forward read error rates for primers: ", pcr_primers, " and flowcell: ", fcid))
} else if (read == "R"){
    filts <- list.files(input_dir, pattern="R2_001.*", full.names = TRUE)
    message(paste0("Modelling reverse read error rates for primers: ", pcr_primers, " and flowcell: ", fcid))
} else {
    stop ("read must be F or R!")
}
# Subset fastqs to just the relevent pcr primers
filts <- filts[str_detect(filts,paste0(pcr_primers, "(-|_|$)"))]
message(paste0(length(filts), " fastq files to process for primers: ", pcr_primers, " and flowcell: ", fcid))

# Learn error rates from a subset of the samples and reads (rather than running self-consist with full dataset)
err <-  dada2::learnErrors(filts, multithread = multithread, nbases = nbases,
                            randomize = randomize, qualityType = "FastqQuality", verbose=TRUE)
#write out errors for diagnostics
if(write_all){
    write_csv(as.data.frame(err$trans), paste0(qc_dir, "/", fcid, "_err",read,"_observed_transitions.csv"))
    write_csv(as.data.frame(err$err_out), paste0(qc_dir, "/", fcid, "_err",read,"_inferred_errors.csv"))
}

##output error plots to see how well the algorithm modelled the errors in the different runs
p1 <- dada2::plotErrors(err, nominalQ = TRUE) + ggtitle(paste0(pcr_primers, " ", fcid, " Forward Reads"))
pdf(paste0(qc_dir,"/",fcid,"_", pcr_primers, "_", read,"_errormodel.pdf"), width = 11, height = 8 , paper="a4r")
    plot(p1)
try(dev.off(), silent=TRUE)