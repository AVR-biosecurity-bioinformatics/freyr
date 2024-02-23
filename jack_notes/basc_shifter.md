# Working notes on how to run the R pipeline in BASC using Shifter and the piperline Docker container

This will run in the `container_test` directory in `/home/js7t/personal/piperline_tests/`. It also runs interactively but will be modified to run with a slurm script.

    sinteractive

Change into the main directory you wish to make the project in

    cd /group/pathogens/IAWS/Personal/JackS/piperline_tests

    # Clone the repository from my version
    git clone https://github.com/jackscanlan/piperline.git container_test && cd container_test

Copy example demultiplexed `.fastq` and samplesheet from another metabarcoding project:

    cp -r /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3 ./data/

    # replace old samplesheet with new samplesheet (just for this example)
    rm ./data/JDYG3/SampleSheet.csv
    cp /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/SampleSheet_JDYG3.csv ./data/JDYG3/SampleSheet.csv

Pull reference database from example directoru; Zenodo download doesn't work: `/group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/reference`

    cp -r /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/reference .

Enter the container:

    module load shifter

    # pull container if you haven't already
    shifterimg pull docker:jackscanlan/piperline:0.0.2

    # enter container
    # this should be altered to the most recent version of the container
    shifter --image jackscanlan/pipeline:0.0.2


Run R script (`basc_shifter.R`) that runs the pipeline inside the container:

    Rscript ./jack_notes/basc_shifter.R
