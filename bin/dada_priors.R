#!/usr/bin/env Rscript

# check variables defined

### run R code


readRDS(denoise_fwd[stringr::str_detect(denoise_fwd,  paste0(.x,"_dada1F.rds"))])$sequence

# Only keep the ones that appear across more than one sample
priors <- unlist(process$priors)
priors <- names(table(priors))[table(priors) > 1]