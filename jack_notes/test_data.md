Aim: To make a test dataset that contains a minimal number of read pairs, that can be used to quickly run the pipeline for development purposes. 

Will base this on the `JDYG3` flow cell data available in the `/group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3/` directory.

> Note: Need minimum 1000 reads in a sample after filtering, so set subsampling to 5000 reads.

> Note: Also need `Undetermined` reads to correctly produce log PDFs, as well as need `RunInfo.xml` file for flowcell QC. 

Workflow:
    
    mkdir -p /group/pathogens/IAWS/Personal/JackS/piperline_tests/test_data && \
        cd /group/pathogens/IAWS/Personal/JackS/piperline_tests/test_data
    
    mkdir full_data JDYG3

    cp -r /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3/*.fastq.gz full_data

    cp full_data/JDYG3_jm00{1,2}A_* JDYG3
    cp full_data/JDYG3_Undetermined* JDYG3

    module load seqtk
    cd JDYG3
    for filename in *.fastq.gz; do
        seqtk sample $filename -s1 5000 | gzip -c > sub_${filename}
    done
    rm JDYG3_*

    for filename in *.fastq.gz; do 
        mv -- "$filename" "${filename#*_}"
    done

    cp /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/SampleSheet_JDYG3.csv .
    
    cp /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3/RunInfo.xml .

`SampleSheet_JDYG3.csv` then edited to remove non-"jm00{1,2}" samples and reuploaded as `SampleSheet_test.csv`, and converted with `dos2unix`. 


Also copy `runParameters.xml` from other directory:

    cp /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3/runParameters.xml .

