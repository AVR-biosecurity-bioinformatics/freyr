Aim: To make a test dataset that contains a minimal number of read pairs, that can be used to quickly run the pipeline for development purposes. 

Will base this on the `JDYG3` flow cell data available in the `/group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3/` directory.

Workflow:
    
    mkdir -p /group/pathogens/IAWS/Personal/JackS/piperline_tests/test_data && \
        cd /group/pathogens/IAWS/Personal/JackS/piperline_tests/test_data
    
    mkdir full_data subsamples

    cp -r /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3/*.fastq.gz full_data

    cp full_data/JDYG3_jm00{1,2}A_* subsamples

    module load seqtk
    cd subsamples
    for filename in *.fastq.gz; do
        seqtk sample $filename -s1 1000 | gzip -c > sub_${filename}
    done
    rm JDYG3_*

    for filename in *.fastq.gz; do 
        mv -- "$filename" "${filename#*_}"
    done

    cp /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/SampleSheet_JDYG3.csv .

`SampleSheet_JDYG3.csv` then edited to remove non-"jm00{1,2}" samples and reuploaded as `SampleSheet_test.csv`, and converted with `dos2unix`. 


Also copy `runParameters.xml` from other directory:

    cp /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3/runParameters.xml .

