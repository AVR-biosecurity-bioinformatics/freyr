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



Also want to select only reads from Undetermined file that match the following disallowed index combination: "ACCAATGC+CACCACTA"

    cd /home/js7t/personal/nextflow_tests/piperline_nextflow/test_data

    mkdir new_ud && cd new_ud

    # pull headers from fastq
    zcat /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3/JDYG3_Undetermined_S0_R1_001.fastq.gz | \
    awk 'sub(/^@/, "")' - > R1_headers.txt
    zcat /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3/JDYG3_Undetermined_S0_R2_001.fastq.gz | \
    awk 'sub(/^@/, "")' - > R2_headers.txt

    # EOF index combinations
    cat <<EOF > indices.txt
    GAGACGAT+GTTCTCGT
    GAGACGAT+ATGCACGA
    GAGACGAT+CGCTCTAT
    GAGACGAT+AGAGTAGC
    GACGATCT+GTTCTCGT
    GACGATCT+ATGCACGA
    GACGATCT+CGCTCTAT
    GACGATCT+AGAGTAGC
    ACGGAACA+GTTCTCGT
    ACGGAACA+ATGCACGA
    ACGGAACA+CGCTCTAT
    ACGGAACA+AGAGTAGC
    ACGTTCAG+GTTCTCGT
    ACGTTCAG+ATGCACGA
    ACGTTCAG+CGCTCTAT
    ACGTTCAG+AGAGTAGC
    EOF

    # grep headers matching index
    grep -F -f indices.txt R1_headers.txt > R1_bad.txt
    grep -F -f indices.txt R2_headers.txt > R2_bad.txt

    # use BBMap to subsample reads with disallowed headers

    module load BBMap/38.98-GCC-11.2.0

    filterbyname.sh \
    in=/group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3/JDYG3_Undetermined_S0_R1_001.fastq.gz \
    in2=/group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3/JDYG3_Undetermined_S0_R2_001.fastq.gz \
    out=./bad_R1.fastq.gz \
    out2=./bad_R2.fastq.gz \
    ow=t \
    substring=name \
    include=t \
    names=R1_bad.txt

    # concatenate "bad reads" to end of subsampled Undetermined reads
    zcat bad_R1.fastq.gz ../JDYG3/JDYG3_Undetermined_S0_R1_001.fastq.gz | gzip -c > ./JDYG3_Undetermined_S0_R1_001.fastq.gz
    zcat bad_R2.fastq.gz ../JDYG3/JDYG3_Undetermined_S0_R2_001.fastq.gz | gzip -c > ./JDYG3_Undetermined_S0_R2_001.fastq.gz

    # then moved into old directory and replaced existing files