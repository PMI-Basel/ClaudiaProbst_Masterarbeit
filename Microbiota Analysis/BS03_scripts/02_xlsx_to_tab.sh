#!/bin/bash

## --------------------------------------------------------------------
## 02 | xlsx to tab - xlsx2csv
## --------------------------------------------------------------------

#load modules
module load foss/2018b #interpreters
module load xlsx2csv/0.7.4-foss-2018b-Python-3.6.6

#arguments
input=$1
output=$2

#convert the design file from xlsx to tab
xlsx2csv ${input} ${output} -d 'tab'
