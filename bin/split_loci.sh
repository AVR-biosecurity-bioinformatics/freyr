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

# change "I" in primer seq to "N", if present
FWD_PRIMER=${3/I/N}
REV_PRIMER=${4/I/N}

# get lengths of primers
FWD_LEN=${#FWD_PRIMER}
REV_LEN=${#REV_PRIMER}

# set variables for BBDuk based on relative lengths of primers
if [ "$FWD_LEN" -ge "$REV_LEN" ]; then
    KMER_LEN=$FWD_LEN
    MINK_LEN=$REV_LEN
else
    KMER_LEN=$REV_LEN
    MINK_LEN=$FWD_LEN
fi

module load BBMap/38.98-GCC-11.2.0

bbduk.sh \
in=${1} \
in2=${2} \
literal=${FWD_PRIMER},${REV_PRIMER} \
out=reject1.fastq.gz \
out2=reject2.fastq.gz \
outm=${7}_${6}_${5}_R1.fastq.gz \
outm2=${7}_${6}_${5}_R2.fastq.gz \
restrictleft=${KMER_LEN} \
copyundefined=true \
k=${MINK_LEN} \
rcomp=t \
stats=split_loci_stats_${7}_${6}_${5}.txt \
lhist=split_loci_lhist_${7}_${6}_${5}.txt

## count reads in input files
# forward reads
if [[ $1 == *.gz ]]; then
    R1_IN=$(zcat $1 | wc -l)
else 
    R1_IN=$(cat $1 | wc -l)
fi
echo $(( $R1_IN / 4 )) > R1_input.txt

# reverse reads
if [[ $2 == *.gz ]]; then
    R2_IN=$(zcat $2 | wc -l)
else 
    R2_IN=$(cat $2 | wc -l)
fi
echo $(( $R2_IN / 4 )) > R2_input.txt

## count reads in output files
# forward reads
R1_OUT=$(zcat ${7}_${6}_${5}_R1.fastq.gz | wc -l)
echo $(( $R1_OUT / 4 )) > R1_output.txt

# reverse reads
R2_OUT=$(zcat ${7}_${6}_${5}_R2.fastq.gz | wc -l)
echo $(( $R2_OUT / 4 )) > R2_output.txt