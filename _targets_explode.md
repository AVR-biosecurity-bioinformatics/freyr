## Analysis of `_targets.R` file for conversion into Nextflow modules

### Notes
- Have to load require packages at the start of each process script
- Use `jackscanlan/piperline` container for each process (at the moment); see [this](https://www.nextflow.io/docs/latest/container.html)
    - Need to config for Shifter ala. here: https://www.nextflow.io/docs/latest/shifter.html 
    - Safer at this stage to define the specific container for each process, rather than workflow-wide, as will probably replace processes and modules later or add additional modules that don't run in the `piperline` container
- `tar_file()` is roughly equivalent to a Nextflow channel, I think
- `tar_target()` is equivalent to a Nextflow process


### Targets

> Note: Targets by themselves (eg. `samdf`) are actually `tar_target(samdf)` in the code.

#### "Parameter setup" targets

`samdf`: creates a default sample df, loads user samplesheet, checks the input, and asserts essential parameters. 

`params`: creates a default parameter df, loads user params, checks the input, 

`tar_file(params_primer_path)`: creates a params .csv with just primer info
- `params_primer`: reads .csv into df
- can be merged

`tar_file(params_readfilter_path)`: creates a params .csv with just read filter info
- `params_readfilter`: reads .csv into df
- can be merged 

`tar_file(params_dada_path)`: creates params .csv for DADA
- `params_dada`: reads .csv into df
- can be merged

`tar_file(params_asvfilter_path)`: creates params .csv for ASV filtering
- `params_asvfilter`: reads .csv into df
- can be merged

`tar_file(params_database_path)`: creates params .csv for reference database
- `params_database`: reads .csv into df
- can be merged

`tar_file(params_ps_path)`: creates params .csv for phyloseq object manip.
- `params_ps`: reads .csv into df

`create_dirs`: creates a series of data, reference, output and sample_data directories with the `step_validate_folders()` function

`tar_files(fastq_path)`: looks for sequencing reads and saves as an object

`temp_samdf1`: calls `step_check_files()` to check if the files in `fastq_path` match the `samdf` samplesheet

`tar_group_by(temp_samdf1_grouped, temp_samdf1, fcid)`: applies `temp_samdf1` but grouped by `fcid` (flow cell ID)
- this could probably be changed in Nextflow somehow to better use the implicit parallelisation

#### "Sequencing QC" targets

`seq_qc`: applies `step_seq_qc()` to `temp_samdf1_grouped`, which uses [`savR` package](https://www.bioconductor.org/packages/3.16/bioc/html/savR.html) to QC the sequencing run from the SAV files (InterOp and RunInfo.xml) 

`switching_qc`: applies `step_switching_calc()` to `temp_samdf1_grouped` to calculate index switching rates and then plots this info (heat map?) 

> `step_sample_qc()` function defined in `functions.R` but not used in pipeline. This would use FastQC on sample reads. 

> `step_multiqc` function defined in `functions.R` but not used in pipeline. This would use `ngsReports` package to write a QC report.

#### "Demultiplex and trim primers" targets

`primer_trim`: applies `step_primer_trim()` to files in `temp_samdf1`, trimming primer sequences from reads using the `trim_primers()` function
- in the future, could replace with a process using (non-R) software to trim primers

`primer_trim_path`: returns the file path for the primer-trimmed reads from `primer_trim`

`temp_samdf2`: makes another temp samdf from `temp_samdf1` using `step_demux_samdf()` to update sample and logging sheet to add newly demultiplexed files, and `step_check_files()` to check filenames and missing files and to also hash files

#### "Filter reads" targets

`read_filter`: filters and trims reads (from `primer_trim`) using `step_filter_reads()`, which uses `dada2::filterAndTrim()`

`read_filter_path`: returns the file path for the filtered reads from `read_filter`

`temp_samdf3`: makes another temp samdf from `temp_samdf2` using `step_check_files()` to check filenames and missing files and to also hash files

`prefilt_qualplots`: creates quality plots for pre-filter reads using `plot_read_quals()`

`write_prefilt_qualplots`: writes `prefilt_qualplots` plots to disk
- merge with above into one process

`postfilt_qualplots`: creates quality plots for post-filter reads using `plot_read_quals()`

`write_postfilt_qualplots`: writes `postfilt_qualplots` plots to disk
- merge with above into one process

#### "Infer sequence variants with DADA2" targets
- Could split fwd and rev reads into different 'streams' of processes

`tar_group_by(temp_samdf3_grouped, temp_samdf3, fcid, pcr_primers)`: groups `temp_samdf3` by `fcid` and `pcr_primers`

`tar_group_by(temp_samdf3_grouped_sample, temp_samdf3, fcid, pcr_primers, sample_id)`: groups `temp_samdf3` by `fcid`, `pcr_primers` and `sample_id`

`error_model_fwd`: runs `step_errormodel()` on forward reads per flow cell, which uses `dada2::learnErrors()` to learn base errors from data with `temp_samdf3_grouped`. Also outputs plots of the error models

`error_model_rev`: same as `error_model_fwd` but runs on reverse reads

`denoise_fwd`: run `step_dada2_single2()` on forward reads per sample, which uses `dada2::dada()` to denoise reads using error models from `error_model_fwd`

`denoise_rev`: same as `denoise_fwd` but runs on reverse reads

`priors_fwd`: pulls 'priors' from `denoise_fwd` output per sample for forward reads

`priors_rev`: same as `priors_fwd` but for reverse reads

`denoise2_fwd`: same as `denoise_fwd` but uses priors from `priors_fwd`

`denoise2_rev`: same as `denoise2_fwd` but for reverse reads

`dada`: merges overlapping denoised read pairs using `dada2::mergePairs()`; else concatenates reads pairs that don't overlap into single reads with 10 Ns in between. Outputs 'seqtab' (sequence tables) by primer sequence
- this effectively creates the ASVs: all sequences in these output files ('seqtabs') are treated as ASVs

`dada_path`: returns the file paths for `dada` output (ie. merged and concatenated reads)

#### "Filter ASVs" targets

`filtered_seqtab`: filters ASVs per flow cell/primer combo using `step_filter_asvs()`, which uses `taxreturn::get_binding_position()`, `taxreturn::subset_model()`, `dada2::removeBimeraDenovo()`, `Biostrings::DNAStringSet()`, `taxreturn::map_to_model()` and `taxreturn::codon_filter()`. Also produces some plots of the filtering output.  

`filtered_seqtab_path`: returns the file paths for `filtered_seqtab` output

`write_seqtab_summary`: writes `ASV_cleanup_summary.csv` to disk from `filtered_seqtab`

`write_seqtab_qualplots`: writes `ASV_cleanup_summary.pdf` to disk from `filtered_seqtab`

`merged_seqtab_path`: merges filtered seqtabs across all loci (ie. per primer) into a single seqtab using `dada2::mergeSequenceTables()`

#### "Assign taxonomy" targets

`tar_file(idtaxa_db_tracked)`: pulls IDTAXA database location into file 

`tar_file(ref_fasta_tracked)`: pulls reference database .fasta location info file

`tax_idtaxa`: runs `step_idtaxa()`, which uses `DECIPHER::IdTaxa()` to classify ASVs according to `idtaxa_db_tracked`

`idtaxa_path`: writes out idtaxa objects from `tax_idtaxa`

`tax_blast_path`: runs `step_blast_tophit()` , which uses `taxreturn::blast_assign_species()` to run `blastn` for each ASV and select the top hit

`joint_tax`: aggregates the output of `idtaxa_path` and `tax_blast_path` using `coalesce_tax()` and `step_join_tax_blast()` (which uses `seqateurs::na_to_unclassified()`)

`merged_tax`: merge ASVs in `joint_tax` output if they're the same, I think?

`assignment_plot`: creates an assignment plot (bar graph with colours?) using `taxreturn::blast_top_hit()`

`write_assignment_plot`: writes `assignment_plot` to file as .pdf

`tax_summary`:

#### Taxonomic analysis targets

`ps`:

`ps_summary`:

`accumulation_curve`:

`ps_filtered`:

`ps_filt_summary`:

`read_tracking`:
