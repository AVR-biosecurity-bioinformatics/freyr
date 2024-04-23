
# commands from: https://github.com/jackscanlan/piperline/blob/main/jack_notes/basc_run.md

library(renv)
library(pak)
library(targets)
library(tarchetypes)

source("_targets_packages.R")
source("R/functions_old.R") # uses the old functions.R file before I changed things for Nextflow
source("R/themes.R")

runs <- dir("data/") #Find all directories within data
SampleSheet <- list.files(paste0("data/", runs), pattern= "SampleSheet", full.names = TRUE)
runParameters <- list.files(paste0("data/", runs), pattern= "[Rr]unParameters.xml", full.names = TRUE)

# Create samplesheet containing samples and run parameters for all runs
samdf <- create_samplesheet(SampleSheet = SampleSheet, runParameters = runParameters, template = "V4") %>%
distinct()

# Check that sample_ids contain fcid, if not; attatch
samdf <- samdf %>%
mutate(sample_id = case_when(
    !str_detect(sample_id, fcid) ~ paste0(fcid,"_",sample_id),
    TRUE ~ sample_id
))

# Check that samples match samplesheet
fastqFs <- 
    purrr::map(list.dirs("data", recursive=FALSE),
                    list.files, pattern="_R1_", full.names = TRUE) %>%
    unlist() %>%
    str_remove(pattern = "^(.*)\\/") %>%
    str_remove(pattern = "(?:.(?!_S))+$")

# Filter undetermined reads from sample sheet
fastqFs <- fastqFs[!str_detect(fastqFs, "Undetermined")]

# Check for fastq files that are missing from samplesheet
if (length(setdiff(fastqFs, samdf$sample_id)) > 0) {warning("The fastq file/s: ", setdiff(fastqFs, samdf$sample_id), " are not in the sample sheet") }

# Check for sample_ids that dont have a corresponding fastq file
if (length(setdiff(samdf$sample_id, fastqFs)) > 0) {
warning(paste0("The fastq file: ",
                setdiff(samdf$sample_id, fastqFs),
                " is missing, dropping from samplesheet \n")) 
samdf <- samdf %>%
    filter(!sample_id %in% setdiff(samdf$sample_id, fastqFs))
}

# Write out sample tracking sheet
write_csv(samdf, "sample_data/Sample_info.csv")

# Add primers to sample sheet
samdf <- samdf %>%
    mutate(
        target_gene = "COI;EIF3L",
        pcr_primers = "fwhF2-fwhR2nDac;EIF3LminiF4-EIF3lminiR4",
        for_primer_seq = "GGDACWGGWTGAACWGTWTAYCCHCC;GATGCGYCGTTATGCYGATGC",
        rev_primer_seq = "GTRATWGCHCCIGCTAADACHGG;TTRAAYACTTCYARATCRCC"
        )

write_csv(samdf, "sample_data/Sample_info.csv")



# Params to add in step_add_parameters
params <- tibble(
    # Primer parameters
    pcr_primers = c("fwhF2-fwhR2nDac", "EIF3LminiF4-EIF3lminiR4"),
    target_gene=c("COI", "EIF3L"),
    max_primer_mismatch=1,

    # Read filtering
    read_min_length = 20,
    read_max_length = Inf,
    read_max_ee = 1,
    read_trunc_length = 150,
    read_trim_left = 0,
    read_trim_right = 0,
    
    # ASV filtering
    asv_min_length = c(195, 207),
    asv_max_length = c(215, 227),
    high_sensitivity = TRUE,
    concat_unmerged = FALSE,
    genetic_code = c("SGC4", "SGC0"),
    coding = c(TRUE, TRUE),
    phmm = c("reference/Bactrocera_COI.rds", "reference/Bactrocera_EIF3L.rds"),
    
    # Taxonomic assignment
    idtaxa_db = c("reference/COI_idtaxa.rds","reference/EIF3L_idtaxa.rds"),
    ref_fasta = c("reference/COI_hierarchial.fa.gz", "reference/EIF3L_hierarchial.fa.gz"),
    idtaxa_confidence = 60,
    run_blast=TRUE,
    blast_min_identity = 97,
    blast_min_coverage = 90,
    target_kingdom = "Metazoa",
    target_phylum = "Arthropoda",
    target_class = "Insecta",
    target_order = "Diptera",
    target_family = NA, 
    target_genus = NA,  
    target_species = NA,  
    
    # Sample & Taxon filtering
    min_sample_reads = c(1000, 1000),
    min_taxa_reads= NA, 
    min_taxa_ra = c(1e-4, 1e-4),
    
    # General pipeline parameters
    threads = 1
)

write_csv(params, "sample_data/loci_params.csv")


# run pipeline

tar_make(script = "_targets_rank7.R")
