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

module load BBMap/38.98-GCC-11.2.0

filterbysequence.sh \
in= ${1} \
in2= ${2} \
literal= ${3},${4} \
out= test1_${5}.fastq.gz \
out2= test2_${5}.fastq.gz \
include=true