#!/bin/bash

## --------------------------------------------------------------------
## 01 | Quality Control - FastQC
## --------------------------------------------------------------------

#load modules
module load FastQC/0.11.8-Java-1.8

#arguments
input=$1
output=$2
output_seq=$3

#get path
output_path=$(dirname ${output})

#quality control
fastqc -t 4 -q ${input} -o ${output_path}

#remove zipped files
file=$(echo ${output} | sed 's/html/zip/')
rm ${file}

#get the number of raw sequences
raw_lines=$(gunzip -c ${input} | wc -l)
raw_seqs=$(echo ${raw_lines} / 4 | bc) #each seq has 4 lines
echo ${raw_seqs} > ${output_seq}
