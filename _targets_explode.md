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

`switching_qc`: applies `step_switching_calc()` to `temp_samdf1_grouped` to calculate index switching rates and then plot 

> `step_sample_qc()` function defined in `functions.R` but not used in pipeline. This would use FastQC on sample reads. 

> `step_multiqc` function defined in `functions.R` but not used in pipeline. This would use `ngsReports` package to write a QC report.

#### "Demultiplex and trim primers" targets

`primer_trim`: 

`primer_trim_path`:

`temp_samdf2`

#### "Filter reads" targets

`read_filter`:

`read_filter_path`: 

`temp_samdf3`:

`prefilt_qualplots`:

`write_prefilt_qualplots`:

`postfilt_qualplots`:

`write_postfilt_qualplots`:

#### "Infer sequence variants with DADA2" targets
- Could split fwd and rev reads into different 'streams' of processes

`tar_group_by(temp_samdf3_grouped, temp_samdf3, fcid, pcr_primers)`:

`tar_group_by(temp_samdf3_grouped_sample, temp_samdf3, fcid, pcr_primers, sample_id)`:

`error_model_fwd`:

`error_model_rev`:

`denoise_fwd`:

`denoise_rev`:

`priors_fwd`:

`priors_rev`:

`denoise2_fwd`:

`denoise2_rev`:

`dada`:

`dada_path`:

#### "Filter ASVs" targets

`filtered_seqtab`:

`filtered_seqtab_path`:

`write_seqtab_summary`:

`write_seqtab_qualplots`:

`merged_seqtab_path`:

#### "Assign taxonomy" targets

`tar_file(idtaxa_db_tracked)`:

`tar_file(ref_fasta_tracked)`:

`tax_idtaxa`:

`idtaxa_path`:

`tax_blast_path`:

`joint_tax`:

`merged_tax`:

`assignment_plot`:

`write_assignment_plot`:

`tax_summary`:

#### Taxonomic analysis targets

`ps`:

`ps_summary`:

`accumulation_curve`:

`ps_filtered`:

`ps_filt_summary`:

`read_tracking`:
