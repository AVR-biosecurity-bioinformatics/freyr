#!/bin/bash
# AUTHOR: alexander.piper@agriculture.vic.gov.au
# This script calculates the index switching rate using correctly determined and undetermined reads
# First the undetermined reads are split into potentially switched reads (combination of applied indexes) or "other" reads
# The switch rate is then calculated as the abundance ration of switched to correctly demultiplexed reads
# NOTE: Both the demultiplexed reads and undetermined reads must be in the same folder to run this script

# Get a list of fastq files
fq=$(ls | grep "R1_001.fastq.gz" | sort | grep -v 'Undetermined' )
undetermined=$(ls | grep 'Undetermined' | grep 'R1')

# Count correctly determined reads from all fastq files
[ -e determined_counts.txt ] && rm determined_counts.txt
touch determined_counts.txt
for f in ${fq} ;do 
echo ${f}
zcat ${f} | grep '^@M' | rev | cut -d':' -f 1 | rev | sort | uniq -c | sort -nr  | sed 's/+/ /' | sed 's/^ *//g' >> determined_counts.txt
done

# Get all potential switched combinations of used indexes
index1=$(cat determined_counts.txt | cut -d' ' -f 2)
index2=$(cat determined_counts.txt | cut -d' ' -f 3)

[ -e all_combinations.txt ] && rm all_combinations.txt
touch all_combinations.txt
for i in ${index1}
do
  for j in ${index2}
  do
    if [ "$i" \< "$j" ]
    then
     echo $i $j >> all_combinations.txt
    else
     echo 'test'
    fi
  done
done

# Count number of undetermined reads
zcat ${undetermined} | grep '^@M' | rev | cut -d':' -f 1 | rev | sort | uniq -c | sed 's/+/ /' | sed 's/^ *//g'  > undetermined_counts.txt
cat undetermined_counts.txt | cut -d' ' -f 2,3 > undetermined_index.txt

# Count number of correctly demultiplexed reads
correct_counts=$(cat determined_counts.txt | cut -d' ' -f 1 | awk '{ SUM += $1} END { print SUM }')

# Count number of switched reads
comm -12 <(sort all_combinations.txt) <(sort undetermined_index.txt) > switched_indexes.txt
switched_counts=$(grep -f "switched_indexes.txt" "undetermined_counts.txt" | cut -d' ' -f 1 | awk '{ SUM += $1} END { print SUM }')

# Count number of other reads (these can be sequencing errors, PhiX and other junk)
other_counts=$(grep -v -f "switched_indexes.txt" "undetermined_counts.txt" | cut -d' ' -f 1 | awk '{ SUM += $1} END { print SUM }')

# Calculate switch rate (in percentage)
calc(){ awk "BEGIN { print "$*" }"; }
switch_rate=$(calc $switched_counts/$correct_counts)
switch_rate_perc=$(calc $switched_counts/$correct_counts*100)

# Clean up temporary files 
rm determined_counts.txt
rm all_combinations.txt
rm undetermined_counts.txt
rm undetermined_index.txt
rm switched_indexes.txt

# Print results
echo "Correctly demultiplexed reads: ${correct_counts}"
echo "Switched reads: ${switched_counts}"
echo "Other undetermined reads: ${other_counts}"
echo "Index switching rate: ${switch_rate} (${switch_rate_perc}%)"
