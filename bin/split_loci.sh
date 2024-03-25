#!/bin/bash
set -e
set -u
# testing syntax
# echo "The name of my reads are ${1} and ${2}, while projectDir is ${3}!"

# ### use BBMap to split 

## args are the following:
# $1 = reads[0], aka. fwd read path
# $2 = reads[1], aka. rev read path
# $3 = meta.for_primer_seq, aka. fwd primer sequence
# $4 = meta.rev_primer_seq, aka. rev primer sequence
# $5 = meta.target_gene, aka. name of target gene
# $6 = meta.sample_id, aka. sample ID

# change "I" in primer seq to "N"
FWD_PRIMER=${3/I/N}
REV_PRIMER=${4/I/N}

# get length of fwd primer
FWD_LEN=${#FWD_PRIMER}
# get length of rev primer
REV_LEN=${#REV_PRIMER}

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
out=no_match1.fq.gz \
out2=no_match2.fq.gz \
outm=${6}_${5}.fastq.gz \
outm2=${6}_${5}.fastq.gz \
restrictleft=${KMER_LEN} \
copyundefined=true \
k=${MINK_LEN}


