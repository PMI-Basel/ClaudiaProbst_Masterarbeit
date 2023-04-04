#!/bin/bash

## --------------------------------------------------------------------
## 05 | barcode files
## --------------------------------------------------------------------

#arguments
design=$1
output=$2

#get taxa & run
taxa=$(basename ${output} | cut -f1 -d "_")
run=$(basename ${output} | cut -f2 -d "_")

 #create files
  grep ${taxa} ${design} |\
  grep ${run} |\
  awk '{print $4, $5}' | sort | uniq |\
  sed 's/^/>/' | sed 's/ /\n/g' \
  > ${output}

  #if taxa-run combination isn't exisiting, an empty file is created
