# Working notes on how to run the R pipeline in BASC

This will run in the `dros_test` directory in `/group/pathogens/IAWS/Personal/JackS/piperline_tests`. 

    # Change into the main directory you wish to make the project in
    cd /group/pathogens/IAWS/Personal/JackS/piperline_tests

    # Clone the repository from my version
    git clone https://github.com/jackscanlan/piperline.git dros_test
    cd dros_test

Copy example demultiplexed `.fastq` and samplesheet from another metabarcoding project:

    cp -r /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3 ./data/

    # replace old samplesheet with new samplesheet (just for this example)
    rm ./data/JDYG3/SampleSheet.csv
    cp /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/SampleSheet_JDYG3.csv ./data/JDYG3/SampleSheet.csv

Pull reference database from example directoru; Zenodo download doesn't work: `/group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/reference`

    cp -r /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/reference .

Load interactive job, modules and R session:

    # Create new interactive SLURM session
    sinteractive --ntasks=1 --cpus-per-task=1 --mem-per-cpu=10GB --time=72:00:00

    # load modules
    module load R/4.2.0-foss-2021b
    module load pkgconfig/1.5.1-GCCcore-9.3.0-Python-3.8.2
    module load GDAL/3.3.0-foss-2021a
    module load BLAST+/2.11.0-gompi-2020a
    module load Pandoc/2.5
    module load ZeroMQ/4.3.2-GCCcore-9.3.0

    # make sure default R library exists
    mkdir -p ~/R_libs/default

    # create .Rprofile if it doesn't exist already
    echo '.libPaths("~/R_libs/default")' > ~/.Rprofile
    echo 'options(renv.config.pak.enabled = TRUE)' >> ~/.Rprofile
    # do above so renv uses pak to install packages (need to load pak first)

    # Start R session
    R

    # load user source file
    > source("~/.Rprofile")

I'm avoiding using the packages `crew`, `mirai` and `nanonext` as the latter has a dependency that is missing on BASC. As such, I'm using the `_targets_nocrew.R` targets file and the `renv_nocrew.lock` file to install packages. This file now also doesn't install `taxreturn`.

Install `renv` package (version `1.0.3`):

    > install.packages("renv", version = "1.0.3", repos = "http://cran.rstudio.com/")

Install `pak` package:

    > install.packages("pak", repos = "http://cran.rstudio.com/")

Load packages:

    > library(devtools) # trying this to install taxreturn manually
    > library(renv)
    > library(pak)

Install `taxreturn` manually using `devtools`:

    > install_github("alexpiper/taxreturn@e9dc03a", force = T)

Load taxreturn independently from personal library:

    > library(taxreturn)

Restore packages from `renv_nocrew.lock` file using `renv`:

    > renv::restore(lockfile = "./renv_nocrew.lock")
    # type 'y' to activate project, then "y" again to proceed installation; packages will install from download or cache (takes a while)

`pak` works really well and installs packages in parallel--much faster!

Create new `_targets_packages_nocrew.R` file that doesn't contain `crew`, `mirai` or `nanonext` and source it, installing the packages:

    > source("_targets_packages_nocrew.R")
    
    # Source ancillary functions
    > source("R/functions.R")
    > source("R/themes.R")

Create sample tracking sheet:


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

Add PCR primers to sample sheet:

    # Add primers to sample sheet
    samdf <- samdf %>%
        mutate(pcr_primers = "fwhF2-fwhR2n",
               for_primer_seq = "GGDACWGGWTGAACWGTWTAYCCHCC",
               rev_primer_seq = "GTRATWGCHCCDGCTARWACWGG"
               )

    write_csv(samdf, "sample_data/Sample_info.csv")

Create parameters file for single primer set:

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

Write pipeline with `_targets_nocrew.R` script:

    tar_make(script = "_targets_nocrew.R")

If running with `basc_run.R` script, run these commands first to convert DOS line breaks to UNIX format:

    dos2unix basc_run.R
    dos2unix basc.slurm

