---
title: "Piperline style freyr submission"
output: html_document
date: "2025-01-22"
---

# Introduction

The below code is written for the Agriculture Victoria BASC computing cluster.

This workflow gives examples of how to run freyr through SLURM with a piperline style submission syntax


Install nextflow in your home directory
freyr currently requires a specific version of nextflow that is not available as a BASC module. It also uses a third-party validation plugin that doesn't work with standalone (ie. module-based) distributions of nextflow. Luckily, it is very easy to install nextflow for a specific user on BASC:

```
## this can be done in a login node

# change to your home directory
cd ~

# load Java
module load Java/17.0.6

# install nextflow in current directory
curl -s https://get.nextflow.io | bash

# make nextflow executable
chmod 777 nextflow

# make a home bin directory if you don't already have one
mkdir -p ~/bin

# move nextflow into bin to make it executable from any path
mv nextflow ~/bin

```

NOTE: This only needs to be done once before any particular user uses nextflow for the first time -- you don't need to repeat this step for subsequent runs of this pipeline, or any other nextflow pipeline.


# Clone the freyr github repository

```
# Define the directory you will be running the analysis in
working_dir=/group/pathogens/IAWS/Personal/Alexp/metabarcoding/marine_surveillance #CHANGE TO YOUR DIRECTORY

git clone https://github.com/AVR-biosecurity-bioinformatics/freyr $working_dir

```


## Updating freyr

To update to the latest version of freyr, run the below code in the terminal.

```
cd $working_dir
git pull

# if an error occurs, run:
git stash
git pull
```

# Demultiplex MiSeq run

For this workflow to run, we will need some sequencing runs to work with. If you are working with MiSeq data, it is recommended that the data is demultiplexed again using bcl2fastq, as the miseq does not put indexes in fasta headers by default which is required for the index swtiching calculation.

The below code is written for the Agriculture Victoria BASC computing cluster, and the locations will be different if you are using a different HPC cluster.

```
#load module
module load bcl2fastq2/2.20.0-foss-2018b

#raise amount of available file handles
ulimit -n 4000

###Run1

#Set up input and outputs
inputdir=/group/sequencing/210219_M03633_0489_000000000-JDYG3 #CHANGE TO YOUR SEQ RUN

fcid=$(echo $input_dir | sed 's/^.*-//')
outputdir=$working_dir/data/$fcid
samplesheet=$input_dir/SampleSheet.csv

# convert samplesheet to unix format
dos2unix $samplesheet

#Demultiplex
bcl2fastq -p 12 --runfolder-dir $input_dir \
--output-dir $outputdir \
--sample-sheet $samplesheet \
--no-lane-splitting --barcode-mismatches 1

# Copy other necessary files and move fastqs
cd $outputdir
cp -r $input_dir/InterOp $outputdir
cp $input_dir/RunInfo.xml $outputdir
cp $input_dir/[Rr]unParameters.xml $outputdir
cp $samplesheet $outputdir
mv **/*.fastq.gz $outputdir

# Append fcid to start of sample names if missing
for i in *.fastq.gz; do
  if ! [[ $i == $fcid* ]]; then
  new=$(echo ${fcid} ${i}) #append together
  new=$(echo ${new// /_}) #remove any white space
  mv -v "$i" "$new"
  fi
done

```

## Copy reference databases

Reference databases are stored in the referencedata directory on BASC

```
# Change 'folder-name' to the directory you are running the analysis in
cp /group/referencedata/mspd-db/metabarcoding/arthropod/imappests_coi_18_08_2020/* $working_dir/reference/.
```

# Submitting jobs

The parameters for the run are parsed along with the slurm script.


## Example 1: Single primer set

This is the 'classic' pipeRline style submission, where all fastq files, miseq sample sheets, and runparameters files are found automatically. 

This requires the data to be placed within the subfolders of the data directory. For example, if you are wishing to analyse two flow cells at once the data directory should look something like this:

    root/
      ├── data/
         ├── K739J/
         ├── JDYG3/
         
Paths to files can be relative to the submission location, or absolute paths     

```
# Submit slurm job - CHANGE EMAIL TO YOUR OWN EMAIL
sbatch --mail-user=your.name@email.com --account=pathogens supplementary_scripts/submit_slurm.sh \
--pcr_primers fwhF2-fwhR2n \
--for_primer_seq GGDACWGGWTGAACWGTWTAYCCHCC \
--rev_primer_seq GTRATWGCHCCDGCTARWACWGG \
--target_gene COI \
--max_primer_mismatch 0 \
--read_min_length 20 \
--read_max_length Inf \
--read_max_ee 1 \
--read_trunc_length 150 \
--read_trim_left 0 \
--read_trim_right 0 \
--asv_min_length 195 \
--asv_max_length 215 \
--high_sensitivity TRUE \
--concat_unmerged FALSE \
--genetic_code SGC4 \
--coding TRUE \
--phmm reference/folmer_fullength_model.rds \
--idtaxa_db reference/idtaxa_bftrimmed.rds \
--ref_fasta reference/insecta_hierarchial_bftrimmed.fa.gz \
--idtaxa_confidence 60 \
--run_blast TRUE \
--blast_min_identity 97 \
--blast_min_coverage 90 \
--target_kingdom Metazoa \
--target_phylum Arthropoda \
--target_class NA \
--target_order NA \
--target_family NA \
--target_genus NA \
--target_species NA \
--min_sample_reads 1000 \
--min_taxa_reads NA \
--min_taxa_ra 1e-4
```


# Example 2: Manually point to data folders

Alternatively, you can manually direct the submission script to the locations of oyur reads, sample sheets, and run parameters files. Multiple files or locations can be separated by a comma
This simple example 

if you wish to use multiple values (i.e. multiple reference databases, multiple primers per index) encapsulate them in a '' and separate them with a ;

Note, watch for the naming of the runParameters file, different versions of the miseq software output this file with a capitalized or lowercase R

```
# Submit slurm job - CHANGE EMAIL TO YOUR OWN EMAIL
sbatch --mail-user=your.name@email.com --account=pathogens supplementary_scripts/submit_slurm.sh \
--sample_sheet 'data/K739J/SampleSheet_K739J.csv,data/JDYG3/SampleSheet_JDYG3.csv'  \
--run_parameters 'data/K739J/RunParameters.xml,data/JDYG3/runParameters.xml' \
--read_dir 'data/K739J/,data/JDYG3/' \
--pcr_primers fwhF2-fwhR2n \
--for_primer_seq GGDACWGGWTGAACWGTWTAYCCHCC \
--rev_primer_seq GTRATWGCHCCDGCTARWACWGG \
--target_gene COI \
--max_primer_mismatch 0 \
--read_min_length 20 \
--read_max_length Inf \
--read_max_ee 1 \
--read_trunc_length 150 \
--read_trim_left 0 \
--read_trim_right 0 \
--asv_min_length 195 \
--asv_max_length 215 \
--high_sensitivity TRUE \
--concat_unmerged FALSE \
--genetic_code SGC4 \
--coding TRUE \
--phmm reference/folmer_fullength_model.rds \
--idtaxa_db reference/idtaxa_bftrimmed.rds \
--ref_fasta reference/insecta_hierarchial_bftrimmed.fa.gz \
--idtaxa_confidence 60 \
--run_blast TRUE \
--blast_min_identity 97 \
--blast_min_coverage 90 \
--target_kingdom Metazoa \
--target_phylum Arthropoda \
--target_class NA \
--target_order NA \
--target_family NA \
--target_genus NA \
--target_species NA \
--min_sample_reads 1000 \
--min_taxa_reads NA \
--min_taxa_ra 1e-4

```


## Example 3: Multiplexed primer sets

Some metabarcoding assays amplify a sample with multiple primer sets, either in a multiplex, or seperately and then pooled together before indexing (poolplexed). In this case each fastq file will contain multiple primers.

To handle this case, you need to provide the parameters for each primer seperated by a ; and encapsulated in a ''. This will enable the pipeline to demultiplex the multiplexed data into the seperate amplicons before further processing. Important! The parameters for each primer must be in the same order as defined in --pcr_primers.

```
# Multiple primers within each sample
sbatch --mail-user=your.name@email.com --account=pathogens supplementary_scripts/submit_slurm.sh \
--sample_sheet data/K739J/SampleSheet_K739J.csv \
--run_parameters data/K739J/RunParameters.xml \
--read_dir data/K739J/ \
--pcr_primers 'fwhF2-fwhR2nDac;EIF3LminiF4-EIF3lminiR4' \
--for_primer_seq 'GGDACWGGWTGAACWGTWTAYCCHCC;GATGCGYCGTTATGCYGATGC' \
--rev_primer_seq 'GTRATWGCHCCIGCTAADACHGG;TTRAAYACTTCYARATCRCC' \
--target_gene 'COI;EIF3L' \
--max_primer_mismatch 0 \
--read_min_length 20 \
--read_max_length Inf \
--read_max_ee 1 \
--read_trunc_length 150 \
--read_trim_left 0 \
--read_trim_right 0 \
--asv_min_length '195;207' \
--asv_max_length '215;227' \
--high_sensitivity TRUE \
--concat_unmerged FALSE \
--genetic_code 'SGC4;SGC0' \
--coding TRUE \
--phmm 'reference/Bactrocera_COI.rds;reference/Bactrocera_EIF3L.rds' \
--idtaxa_db 'reference/COI_idtaxa.rds;reference/EIF3L_idtaxa.rds' \
--ref_fasta 'reference/COI_hierarchial.fa.gz;reference/EIF3L_hierarchial.fa.gz' \
--idtaxa_confidence 60 \
--run_blast TRUE \
--blast_min_identity 97 \
--blast_min_coverage 90 \
--target_kingdom Metazoa \
--target_phylum Arthropoda \
--target_class Insecta \
--target_order Diptera \
--target_family NA \
--target_genus NA \
--target_species NA \
--min_sample_reads 1000 \
--min_taxa_reads NA \
--min_taxa_ra 1e-4

```

## Example 4: Different primer sets per sample

In other cases, each sample may only have a single primer set, but different samples within a run have different primers. For instance in this flow cell, all samples with fwh in their name use the fwhF2-fwhR2n primer set, while all samples with EIF in their name use the EIF3L primers.

In this case, you can use string matching to assign specific primers and parameters to specific samples by placing the string you wish to match within square brackets, followed by the value you wish to add to those samples. Multiple match and value combinatiosn can be seperated by a comma

For instance --pcr_primers '[COI]fwhF2-fwhR2n' will add the value fwhF2-fwhR2n to the pcr_primers parameter for only those samples with 'COI' in their sample name

NOTE: You need to ensure that all samples are covered by the strings you input, as a sample missing parameters cannot be run through freyr

```
# Match primers and parameters to samples based on string matching
sbatch --mail-user=your.name@email.com --account=pathogens supplementary_scripts/submit_slurm.sh \
--sample_sheet data/K3DVL/SampleSheet_K3DVL.csv \
--run_parameters data/K3DVL/RunParameters.xml \
--read_dir data/K3DVL/ \
--pcr_primers '[fwh]fwhF2-fwhR2nDac,[EIF]EIF3LminiF4-EIF3lminiR4' \
--for_primer_seq '[fwh]GGDACWGGWTGAACWGTWTAYCCHCC,[EIF]GATGCGYCGTTATGCYGATGC' \
--rev_primer_seq '[fwh]GTRATWGCHCCIGCTAADACHGG,[EIF]TTRAAYACTTCYARATCRCC' \
--target_gene '[fwh]COI,[EIF]EIF3L' \
--max_primer_mismatch 0 \
--read_min_length 20 \
--read_max_length Inf \
--read_max_ee 1 \
--read_trunc_length 150 \
--read_trim_left 0 \
--read_trim_right 0 \
--asv_min_length '[fwh]195,[EIF]207' \
--asv_max_length '[fwh]215,[EIF]227' \
--high_sensitivity TRUE \
--concat_unmerged FALSE \
--genetic_code '[fwh]SGC4,SGC0' \
--coding TRUE \
--phmm '[fwh]reference/Bactrocera_COI.rds,[EIF]reference/Bactrocera_EIF3L.rds' \
--idtaxa_db '[fwh]reference/COI_idtaxa.rds,[EIF]reference/EIF3L_idtaxa.rds' \
--ref_fasta '[fwh]reference/COI_hierarchial.fa.gz,[EIF]reference/EIF3L_hierarchial.fa.gz' \
--idtaxa_confidence 60 \
--run_blast TRUE \
--blast_min_identity 97 \
--blast_min_coverage 90 \
--target_kingdom Metazoa \
--target_phylum Arthropoda \
--target_class NA \
--target_order NA \
--target_family NA \
--target_genus NA \
--target_species NA \
--min_sample_reads 1000 \
--min_taxa_reads NA \
--min_taxa_ra 1e-4

```


## More complex setups

For more complex setups, you can point the submission script to a premade freyr format sample sheet and loci parameters file

**TODO**

