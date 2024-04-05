Aim: To make a test dataset that contains a minimal number of read pairs, that can be used to quickly run the pipeline for development purposes. 

### new (dual locus data; two flow cells)

Based on tephritid metabarcoding data found here: `/group/pathogens/IAWS/Projects/Metabarcoding/tephritid_metabarcoding/data/`

- Two flow cells (`K77JP` and `K739J`), which are duplicate runs of the same libraries on different flow cells
- Two primer pairs/loci per sample: COI (fwhF2-fwhR2nDac) and EIF3L (EIF3LminiF4-EIF3lminiR4)

- grab four samples (ie. read pairs), for each flow cell
    - subsample each file to 5000 reads using seqtk and seed of 1
    - check indices of each sample and use them to pull all Undetermined reads that match all possible combinations of those indices (allowed and disallowed)
- copy over reference database files from `/group/pathogens/IAWS/Projects/Metabarcoding/tephritid_metabarcoding/reference` -- specifically need:
    - `/group/home/js7t/projects/Metabarcoding/tephritid_metabarcoding/diagnostic_alignments/model/Bactrocera_COI.rds`
    - `/group/home/js7t/projects/Metabarcoding/tephritid_metabarcoding/diagnostic_alignments/model/Bactrocera_EIF3L.rds`
    - `/group/home/js7t/projects/Metabarcoding/tephritid_metabarcoding/reference/COI_hierarchial.fa.gz`
    - `/group/home/js7t/projects/Metabarcoding/tephritid_metabarcoding/reference/EIF3L_hierarchial.fa.gz`
    - `/group/home/js7t/projects/Metabarcoding/tephritid_metabarcoding/reference/COI_idtaxa.rds`
    - `/group/home/js7t/projects/Metabarcoding/tephritid_metabarcoding/reference/EIF3L_idtaxa.rds`

cp /group/home/js7t/projects/Metabarcoding/tephritid_metabarcoding/diagnostic_alignments/model/Bactrocera_COI.rds /home/js7t/personal/nextflow_tests/piperline_nextflow/reference
cp /group/home/js7t/projects/Metabarcoding/tephritid_metabarcoding/diagnostic_alignments/model/Bactrocera_EIF3L.rds /home/js7t/personal/nextflow_tests/piperline_nextflow/reference
cp /group/home/js7t/projects/Metabarcoding/tephritid_metabarcoding/reference/COI_hierarchial.fa.gz /home/js7t/personal/nextflow_tests/piperline_nextflow/reference
cp /group/home/js7t/projects/Metabarcoding/tephritid_metabarcoding/reference/EIF3L_hierarchial.fa.gz /home/js7t/personal/nextflow_tests/piperline_nextflow/reference
cp /group/home/js7t/projects/Metabarcoding/tephritid_metabarcoding/reference/COI_idtaxa.rds /home/js7t/personal/nextflow_tests/piperline_nextflow/reference
cp /group/home/js7t/projects/Metabarcoding/tephritid_metabarcoding/reference/EIF3L_idtaxa.rds /home/js7t/personal/nextflow_tests/piperline_nextflow/reference



Workflow:
    
ORIG_DIR="/group/home/js7t/projects/Metabarcoding/tephritid_metabarcoding/data"
NEW_DIR="/home/js7t/personal/nextflow_tests/piperline_nextflow/test_data"

cd $NEW_DIR
mkdir -p new_dual && cd new_dual # location of new reads
mkdir -p full_data/K77JP && mkdir -p full_data/K739J # locations for complete read samples
mkdir -p K77JP K739J # ultimate locations for subsampled reads

# copy 4 reads pairs from each flow cell into full_data dir
cp -r ${ORIG_DIR}/K77JP/K77JP_{Mock1_,Mock2_,Trap1_,Trap2_}*.fastq.gz full_data/K77JP
cp -r ${ORIG_DIR}/K739J/K739J_{Mock1_,Mock2_,Trap1_,Trap2_}*.fastq.gz full_data/K739J

# copy InterOp folders, RunInfo.xml and RunParameters.xml into ultimate directories
for flowcell in K77JP K739J; do
    cp -r ${ORIG_DIR}/${flowcell}/InterOp ${NEW_DIR}/new_dual/${flowcell}
    cp -r ${ORIG_DIR}/${flowcell}/RunInfo.xml ${NEW_DIR}/new_dual/${flowcell}
    cp -r ${ORIG_DIR}/${flowcell}/RunParameters.xml ${NEW_DIR}/new_dual/${flowcell}
done

## subsample reads
module load seqtk
# subsample first flow cell
cd /home/js7t/personal/nextflow_tests/piperline_nextflow/test_data/new_dual/full_data/K77JP
for filename in *.fastq.gz; do
    seqtk sample $filename -s1 5000 | gzip -c > ../../K77JP/${filename}
done
# subsample second flow cell
cd /home/js7t/personal/nextflow_tests/piperline_nextflow/test_data/new_dual/full_data/K739J
for filename in *.fastq.gz; do
    seqtk sample $filename -s1 5000 | gzip -c > ../../K739J/${filename}
done

# also downloaded /group/home/js7t/projects/Metabarcoding/tephritid_metabarcoding/data/K77JP/SampleSheet_K77JP.csv and K739J equiv. and pruned to Mock1/2 and Trap1/2 samples MANUALLY

## sort out Undetermined reads per flowcell
# pull headers from fastq files
for flowcell in K77JP K739J; do
    zcat ${ORIG_DIR}/${flowcell}/${flowcell}_Undetermined_S0_R1_001.fastq.gz | \
    awk 'sub(/^@/, "")' - > ${NEW_DIR}/new_dual/${flowcell}/R1_headers.txt
    zcat ${ORIG_DIR}/${flowcell}/${flowcell}_Undetermined_S0_R2_001.fastq.gz | \
    awk 'sub(/^@/, "")' - > ${NEW_DIR}/new_dual/${flowcell}/R2_headers.txt
done

# EOF index combinations
cat <<EOF > ${NEW_DIR}/new_dual/indices.txt
TGGCATGT+CCTTGTAG
CTGATCGT+CCTTGTAG
ACTCGTTG+CCTTGTAG
TGCGTAGA+CCTTGTAG
TGGCATGT+GAACATCG
CTGATCGT+GAACATCG
ACTCGTTG+GAACATCG
TGCGTAGA+GAACATCG
TGGCATGT+CACCACTA
CTGATCGT+CACCACTA
ACTCGTTG+CACCACTA
TGCGTAGA+CACCACTA
TGGCATGT+TTGGTGAG
CTGATCGT+TTGGTGAG
ACTCGTTG+TTGGTGAG
TGCGTAGA+TTGGTGAG
EOF

# grep headers matching index; only take 100 reads per file
for flowcell in K77JP K739J; do
    grep -F -f ${NEW_DIR}/new_dual/indices.txt ${NEW_DIR}/new_dual/${flowcell}/R1_headers.txt | \
    head -n 100 > ${NEW_DIR}/new_dual/${flowcell}/R1_bad.txt
    grep -F -f ${NEW_DIR}/new_dual/indices.txt ${NEW_DIR}/new_dual/${flowcell}/R2_headers.txt | \
    head -n 100 > ${NEW_DIR}/new_dual/${flowcell}/R2_bad.txt
done

# use BBMap to subsample reads with disallowed headers
module load BBMap/38.98-GCC-11.2.0

for flowcell in K77JP K739J; do
    filterbyname.sh \
    in=${ORIG_DIR}/${flowcell}/${flowcell}_Undetermined_S0_R1_001.fastq.gz \
    in2=${ORIG_DIR}/${flowcell}/${flowcell}_Undetermined_S0_R2_001.fastq.gz \
    out=${NEW_DIR}/new_dual/${flowcell}/${flowcell}_Undetermined_S0_R1_001.fastq.gz \
    out2=${NEW_DIR}/new_dual/${flowcell}/${flowcell}_Undetermined_S0_R2_001.fastq.gz \
    ow=t \
    substring=name \
    include=t \
    names=${NEW_DIR}/new_dual/${flowcell}/R1_bad.txt
done

# delete header files as they're big
for flowcell in K77JP K739J; do
    rm ${NEW_DIR}/new_dual/${flowcell}/R{1,2}_*.txt
done

# make sure .csv files are okay for unix
for flowcell in K77JP K739J; do
    dos2unix ${NEW_DIR}/new_dual/${flowcell}/SampleSheet*
done

# then download and replace


Indices for both flow cell samples (same across flow cells) are:
index:
TGGCATGT
CTGATCGT
ACTCGTTG
TGCGTAGA
index2:
CCTTGTAG
GAACATCG
CACCACTA
TTGGTGAG









### old (single locus data)

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

