#!/bin/bash

## --------------------------------------------------------------------
## 06.1 | demultiplexing
## --------------------------------------------------------------------

#load modules
module load foss/2018b #interpreters
module load cutadapt/3.4-foss-2018b-Python-3.6.6
module load pigz/2.4-GCCcore-7.3.0

#arguments
input_forward=$1
input_reverse=$2
barcodes=$3
output=$4

#get taxa & run
taxa=$(basename ${barcodes} | cut -f1 -d "_")
run=$(basename ${barcodes} | cut -f2 -d "_")

#create output folder
mkdir -p ${output} 2> /dev/null

#if barcode file is not empty, demultiplex the given barcodes
if [[ `wc -l ${barcodes} | awk '{print $1}'` -gt 0 ]]; then

  #decompress runs
  gunzip -c ${input_forward} > temp_${taxa}_${run}_r1.fastq
  gunzip -c ${input_reverse} > temp_${taxa}_${run}_r2.fastq

  #read barcode-name and barcode-sequence
  while read name;
   do read seq;
   #create output name (remove '>')
   outname=$(echo ${name} | tr -d '>')
   ###forward reads
   #grep reads with the barcode in the header, remove group-separator & save in a new file
   grep -A3 -e ^@\.\*${seq} temp_${taxa}_${run}_r1.fastq | grep -v -- "^--$" > ${output}/${run}_${outname}_F.fastq
   ###reverse reads
   grep -A3 -e ^@\.\*${seq} temp_${taxa}_${run}_r2.fastq | grep -v -- "^--$" > ${output}/${run}_${outname}_R.fastq
  done < ${barcodes}

  #rm temp files
  rm temp_${taxa}_${run}_r1.fastq temp_${taxa}_${run}_r2.fastq

  #create log file containing number of seqs per files
  mkdir -p log/06_demultiplexed 2> /dev/null

  basename ${output}/* | tr -d ".fastq.gz" > names_${taxa}_${run}.txt #filenames
  wc -l ${output}/* | sed '$d' | awk '{print $1}' > lines_${taxa}_${run}.txt #lines of file (remove last line; contains total)
  while read line; do echo $line "/4" | bc; done < lines_${taxa}_${run}.txt > seqs_${taxa}_${run}.txt #seqs of files
  paste names_${taxa}_${run}.txt seqs_${taxa}_${run}.txt > temp_${taxa}_${run}.txt #merge
  echo -e "file sequences" | cat - temp_${taxa}_${run}.txt > log/06_demultiplexed/${taxa}_${run}.log #add col headers

  #rm temp files
  rm names_${taxa}_${run}.txt lines_${taxa}_${run}.txt seqs_${taxa}_${run}.txt temp_${taxa}_${run}.txt

  #gzip files
  pigz -p 40 ${output}/*.fastq

fi
