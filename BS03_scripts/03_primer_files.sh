#!/bin/bash

## --------------------------------------------------------------------
## 03 | primer files
## --------------------------------------------------------------------

#arguments
design=$1
output_forward=$2
output_reverse=$3

#get the run
run=$(basename ${output_forward} | cut -f1 -d "_")

#create files

#microbial primers
#forward
grep ${run} ${design} |\
awk '{print $6, $7}' | sort | uniq |\
sed 's/^/>/' | sed 's/ /\n/g' \
> ${output_forward}
#reverse
grep ${run} ${design} |\
awk '{print $8, $9}' | sort | uniq |\
sed 's/^/>/' | sed 's/ /\n/g' \
> ${output_reverse}

#plant primers
#forward
nchar=$(grep ${run} ${design} | awk '{print $10}' | sort | uniq | wc -c) #empty cell = 1
if [ ${nchar} -gt 1 ]; then # there is a plant primer (cell isn't empty)
  grep ${run} ${design} |\
  awk '{print $10, $11}' | sort | uniq |\
  sed 's/^/>/' | sed 's/ /\n/g' \
  >> ${output_forward}
  #reverse
  grep ${run} ${design} |\
  awk '{print $12, $13}' | sort | uniq |\
  sed 's/^/>/' | sed 's/ /\n/g' \
  >> ${output_reverse}
fi
