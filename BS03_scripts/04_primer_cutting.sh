#!/bin/bash

## --------------------------------------------------------------------
## 04 | primer cutting
## --------------------------------------------------------------------

#load modules
module load foss/2018b #interpreters
module load cutadapt/3.4-foss-2018b-Python-3.6.6

#arguments
input=$1
primer_forward=$2
primer_reverse=$3
output=$4

#get run and direction
run=$(basename ${input} | cut -f1 -d "_" | cut -f1 -d ".")
dir=$(basename ${input} | cut -f2 -d ".")

#create log folder
log_path="log/04_primer_cutting"
mkdir -p ${log_path} 2> /dev/null

#decide if forward or reverse
#forward
#error rate of 0.2
if [[ ${dir} == "r1" ]]; then
    cutadapt -g file:${primer_forward} \
             -o ${output} \
             -e 0.1 \
             --cores 20 \
             ${input} \
             > ${log_path}/${run}.${dir}.log #log file
fi

#reverse
if [[ ${dir} == "r2" ]]; then
    cutadapt -g file:${primer_reverse} \
             -o ${output} \
             -e 0.1 \
             --cores 20 \
             ${input} \
             > ${log_path}/${run}.${dir}.log #log file
fi
