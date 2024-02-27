
### This script just preps the input data for the pipeline like it used to, but doesn't run the targets pipeline like before. 

# load R libraries
source("_targets_packages.R")
source("R/functions.R")
source("R/themes.R")

# input data
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
    mutate(pcr_primers = "fwhF2-fwhR2n",
           for_primer_seq = "GGDACWGGWTGAACWGTWTAYCCHCC",
           rev_primer_seq = "GTRATWGCHCCDGCTARWACWGG"
           )

write_csv(samdf, "sample_data/Sample_info.csv")


# Params to add in step_add_parameters
params <- tibble(
    # Primer parameters
    pcr_primers = "fwhF2-fwhR2n",
    target_gene="COI",
    max_primer_mismatch=0,

    # Read filtering
    read_min_length = 20,
    read_max_length = Inf,
    read_max_ee = 1,
    read_trunc_length = 150,
    read_trim_left = 0, 
    read_trim_right = 0,
    
    # ASV filtering
    asv_min_length = 195, 
    asv_max_length = 215,
    high_sensitivity = TRUE,
    concat_unmerged = FALSE,
    genetic_code = "SGC4",
    coding = TRUE,
    phmm = "reference/folmer_fullength_model.rds",
    
    # Taxonomic assignment
    idtaxa_db = "reference/idtaxa_bftrimmed.rds",
    ref_fasta = "reference/insecta_hierarchial_bftrimmed.fa.gz",
    idtaxa_confidence = 60,
    run_blast=TRUE,
    blast_min_identity = 97,
    blast_min_coverage = 90,
    target_kingdom = "Metazoa",
    target_phylum = "Arthropoda",
    target_class = NA,
    target_order = NA,
    target_family = NA,
    target_genus = NA,
    target_species= NA,
    
    # Sample & Taxon filtering
    min_sample_reads = 1000,
    min_taxa_reads= NA,
    min_taxa_ra = 1e-4, #1e-4 is 0.01%
        
    # General pipeline parameters
    threads = 1
)

write_csv(params, "sample_data/loci_params.csv")
