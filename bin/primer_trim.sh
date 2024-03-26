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

module load BBMap/38.98-GCC-11.2.0

bbduk.sh \
in=${1} \
in2=${2} \
literal=${FWD_PRIMER},${REV_PRIMER} \
out=${7}_${6}_${5}_trim_R1.fastq.gz \
out2=${7}_${6}_${5}_trim_R2.fastq.gz \
outm=reject_R1.fastq.gz \
outm2=reject_R2.fastq.gz \
ktrim=l \
k=10 \
copyundefined=true \
rcomp=t \
tbo=f \
minoverlap=10 \
stats=primer_trim_stats_${7}_${6}_${5}.txt \
lhist=primer_trim_lhist_${7}_${6}_${5}.txt


