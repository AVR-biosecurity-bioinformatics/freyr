# Running classic pipeline with `piperline-multi` Docker container

### `test_data/single` subsampled dataset (5000 reads, 4 samples)
The following is for running the pipeline in BASC. 

    # go to dir (whatever you want)
    cd personal/piperline_tests

    # clone into new dir from classic fork 
    # (which contains multi script but is otherwise kept in sync with classic pipeline)
    git clone https://github.com/biosecurity-bioinformatics/piperline-classic.git multi_test \
    && cd multi_test

    # copy reference databases over from nextflow folder
    cp -r ../../nextflow_tests/piperline_nextflow/reference/* ./reference

    # copy minimal dataset from nextflow folder
    cp -r ../../nextflow_tests/piperline_nextflow/test_data/single/* ./data

Commands for running the pipeline through SLURM script:

    # go to directory where you cloned the repo
    cd ~/personal/piperline_tests/multi_test

    # set params
    MAIL=your.address@email.com
    THREADS=16

    sbatch \
    --mail-user=$MAIL \
    --cpus-per-task=$THREADS \
    supplementary_scripts/basc_shifter_multi.slurm \
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
    --min_taxa_ra 1e-4 \
    --threads $THREADS

1 thread = 2.613 min

16 threads = 2.593 min

There might be an issue with the `crew::crew_controller_local` code at the top of the `_targets.R` script and its interaction with the container environment.

### full `JDYG3` dataset (single flowcell)

    # go to dir (whatever you want)
    cd personal/piperline_tests

    # clone into new dir from classic fork 
    # (which contains multi script but is otherwise kept in sync with classic pipeline)
    git clone https://github.com/biosecurity-bioinformatics/piperline-classic.git multi_test_full \
    && cd multi_test_full

    # copy reference databases over from nextflow folder
    cp -r ../../nextflow_tests/piperline_nextflow/reference/* ./reference

    # copy full flowcell dataset from nextflow folder
    cp -r /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3 ./data