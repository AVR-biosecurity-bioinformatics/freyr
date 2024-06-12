#!/bin/bash
set -e
set -u
## args are the following:
# $1 = reads[0], aka. fwd read path
# $2 = reads[1], aka. rev read path
# $3 = meta.for_primer_seq, aka. fwd primer sequence
# $4 = meta.rev_primer_seq, aka. rev primer sequence
# $5 = meta.pcr_primers, aka. name of PCR primer pair
# $6 = meta.target_gene, aka. name of target gene
# $7 = meta.sample_id, aka. sample ID
# $8 = meta.fcid, aka. flowcell ID

# change "I" in primer seq to "N", if present
FWD_PRIMER=${3/I/N}
REV_PRIMER=${4/I/N}

# ## trying with BBMap (BBDuk)
# module load BBMap/38.98-GCC-11.2.0

# bbduk.sh \
# in=${1} \
# in2=${2} \
# literal=${FWD_PRIMER},${REV_PRIMER} \
# out=${7}_${6}_${5}_trim_R1.fastq.gz \
# out2=${7}_${6}_${5}_trim_R2.fastq.gz \
# outm=reject_R1.fastq.gz \
# outm2=reject_R2.fastq.gz \
# ktrim=l \
# k=10 \
# copyundefined=true \
# rcomp=t \
# tbo=f \
# minoverlap=10 \
# stats=primer_trim_stats_${7}_${6}_${5}.txt \
# lhist=primer_trim_lhist_${7}_${6}_${5}.txt

## trying with cutadapt (preferred)
# module load cutadapt/3.4-GCCcore-10.3.0

# reverse complement the primer sequences
FWD_PRIMER_RC=$(echo ${FWD_PRIMER} | \
    tr GATCRYMKSWHBVDNgatcrymkswhbvdn CTAGYRKMSWDVBHNctagyrkmswdvbhn | \
    rev)
REV_PRIMER_RC=$(echo ${REV_PRIMER} | \
    tr GATCRYMKSWHBVDNgatcrymkswhbvdn CTAGYRKMSWDVBHNctagyrkmswdvbhn | \
    rev)

# run cutadapt in "linked adapter" mode
# see here: https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-paired-end-reads
cutadapt \
    -a ^${FWD_PRIMER}...${REV_PRIMER_RC} \
    -A ^${REV_PRIMER}...${FWD_PRIMER_RC} \
    --discard-untrimmed \
    --rename="{header}" \
    --report=minimal \
    -o ${7}_${6}_${5}_trim_R1.fastq.gz \
    -p ${7}_${6}_${5}_trim_R2.fastq.gz \
    ${1} \
    ${2}

## count reads in output files
# forward reads
if [ -f ${7}_${6}_${5}_trim_R1.fastq.gz ]; then
    R1_OUT_LINES=$(zcat ${7}_${6}_${5}_trim_R1.fastq.gz | wc -l)
    R1_OUT=$(( $R1_OUT_LINES / 4 ))
else 
    R1_OUT=0
fi

# reverse reads
if [ -f ${7}_${6}_${5}_trim_R2.fastq.gz ]; then
    R2_OUT_LINES=$(zcat ${7}_${6}_${5}_trim_R2.fastq.gz | wc -l)
    R2_OUT=$(( $R2_OUT_LINES / 4 ))
else 
    R2_OUT=0
fi

# save as .csv 
echo "primer_trim,$7,$8,$5,$R1_OUT,$R2_OUT" > primer_trim_${7}_${5}_readsout.csv # columns: stage; sample_id; fcid; pcr_primers; fwd_out; fwd_out
