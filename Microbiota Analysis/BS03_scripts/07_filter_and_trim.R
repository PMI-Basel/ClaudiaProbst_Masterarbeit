#!/usr/bin/env Rscript

# get command line arguments as an array
args <- commandArgs(trailingOnly = TRUE)
#saveRDS(args,"args.RDS")

#create snakemake list from passed arguments

#args files
#1    input
#2-4  output
#5-6  log
#7-12 parameters

snakemake <- list(input=args[1], output=args[2:4], log=args[5:6], parameters=args[7:15])

#redirect std out and std err to log file
log_out <- file(snakemake$log[[1]], open="wt")
log_err <- file(snakemake$log[[2]], open="wt")
sink(log_out, type="output")
sink(log_err, type="message")

#libaries
library(dada2)

#catch attributes
input_folder <- snakemake$input[[1]]
output_file <- snakemake$output[[1]]

#get taxa and run
taxa <- sapply(strsplit(basename(output_file), "_", fixed=F), `[`, 2)
run <- sapply(strsplit(basename(output_file), "_", fixed=F), `[`, 3)
run <- gsub(".RDS", "", run)

#####################
# setwd("/Users/j.waechli/Desktop/01_dadapipe")
# taxa <- "bacteria"
# run <- "BS02"
# input_folder <- "06_demultiplexed/bacteria/BS02"
# output_file <- "07_filter_and_trim/RDS/out_bacteria_BS02.RDS"
#####################

#get all taxa-run combinations
DESIGN <- read.delim("01_design/design.tab", sep = "\t")
taxa_runs <- unique(paste(DESIGN$taxa, DESIGN$run_nbr, sep="_"))

#check if taxa-run combination exist, if not skip everything
taxa_run <- paste(taxa, run, sep="_")
if (taxa_run %in% taxa_runs) {

  # ------------------------------------------------------------------------
  # get the samples
  # ------------------------------------------------------------------------

  #Forward and reverse fastq filenames have format: run_FLDnumber_F.fastq.gz and run_FLDnumber_R.fastq.gz
  path_in <- input_folder
  fnFs <- sort(list.files(path_in, pattern="_F.fastq", full.names = TRUE))
  fnRs <- sort(list.files(path_in, pattern="_R.fastq", full.names = TRUE))

  #Extract sample names
  sample.names <- sapply(strsplit(basename(fnFs), "_F.fastq", fixed=F), `[`, 1)

  #Place filtered files in filtered/ subdirectory

  path_out <- file.path(dirname(dirname(output_file)), taxa, run)
  filtFs <- file.path(path_out, paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path_out, paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names

  # ------------------------------------------------------------------------
  # Quality filtering and trimming
  # ------------------------------------------------------------------------

  print("filter and trim started")

  #filter (do each filter separated to see which filter removes how many seqs)
  if (taxa == "bacteria"){

    #create temp files
    for(i in 1:2){
      assign(paste0("temp",i, "Fs"), file.path(path_out, paste0(sample.names, "_F_temp", i, ".fastq.gz")))
      assign(paste0("temp",i, "Rs"), file.path(path_out, paste0(sample.names, "_R_temp", i, ".fastq.gz")))
    }

    #filtering parameters
    truncLen_bac=as.numeric(snakemake$parameters[1:2])
    truncQ_bac=as.numeric(snakemake$parameters[3])
    maxEE_bac=as.numeric(snakemake$parameters[4:5])

    #filtering (default: rm.phix=T)
    out_trunc <- filterAndTrim(fnFs, temp1Fs, fnRs, temp1Rs, truncLen=truncLen_bac, truncQ=truncQ_bac, compress=TRUE, multithread=TRUE, matchIDs = T)
    out_maxN <- filterAndTrim(temp1Fs, temp2Fs, temp1Rs, temp2Rs, maxN=0, compress=TRUE, multithread=TRUE, matchIDs = T)
    out_maxEE <- filterAndTrim(temp2Fs, filtFs, temp2Rs, filtRs, maxEE=maxEE_bac, compress=TRUE, multithread=TRUE, matchIDs = T)

    #remove temp files
    unlink(c(temp1Fs, temp1Rs, temp2Fs, temp2Rs))

    #combine output files
    out <- data.frame(reads.in=out_trunc[,1], trunc=out_trunc[,2], maxN=out_maxN[,2], maxEE=out_maxEE[,2])
    rm(out_trunc, out_maxN, out_maxEE)

  }

  if (taxa == "fungi"){

    #create temp files
    for(i in 1:2){
      assign(paste0("temp",i, "Fs"), file.path(path_out, paste0(sample.names, "_F_temp", i, ".fastq.gz")))
      assign(paste0("temp",i, "Rs"), file.path(path_out, paste0(sample.names, "_R_temp", i, ".fastq.gz")))
    }

    #filtering (default: rm.phix=T)
    #Because ITS region vary in length, we don't use truncLEN, but set a minLen to get rid of spurious very low-length sequences.

    #filtering parameters
    minLen_fun=as.numeric(snakemake$parameters[6])
    truncQ_fun=as.numeric(snakemake$parameters[7])
    maxEE_fun=as.numeric(snakemake$parameters[8:9])

    out_trunc <- filterAndTrim(fnFs, temp1Fs, fnRs, temp1Rs, minLen=minLen_fun, truncQ=truncQ_fun, compress=TRUE, multithread=TRUE, matchIDs = T)
    out_maxN <- filterAndTrim(temp1Fs, temp2Fs, temp1Rs, temp2Rs, maxN=0, compress=TRUE, multithread=TRUE, matchIDs = T)
    out_maxEE <- filterAndTrim(temp2Fs, filtFs, temp2Rs, filtRs, maxEE=maxEE_fun, compress=TRUE, multithread=TRUE, matchIDs = T)

    #remove temp files
    unlink(c(temp1Fs, temp1Rs, temp2Fs, temp2Rs))

    #combine output files
    out <- data.frame(reads.in=out_trunc[,1], trunc=out_trunc[,2], maxN=out_maxN[,2], maxEE=out_maxEE[,2])
    rm(out_trunc, out_maxN, out_maxEE)

  }

  #Avoid to raise an error message if no reads from a sample pass the filter
  missing_fltFs <- c()
  missing_fltRs <- c()
  for (file in filtFs) {if(!file.exists(file)) {missing_fltFs <- c(missing_fltFs, file)}} #collect names of the missing files
  for (file in filtRs) {if(!file.exists(file)) {missing_fltRs <- c(missing_fltRs, file)}}
  filtFs <- setdiff(filtFs, missing_fltFs)
  filtRs <- setdiff(filtRs, missing_fltRs)

  # ------------------------------------------------------------------------
  # output RDS files
  # ------------------------------------------------------------------------

  path_RDS <- file.path(dirname(dirname(output_file)), "RDS") #get only the first part of the path
  name_RDS <- paste0(taxa, "_", run, ".RDS")
  dir.create(path_RDS)
  saveRDS(out, file.path(path_RDS, paste0("out_", name_RDS)))
  saveRDS(filtFs, file.path(path_RDS, paste0("filtFs_", name_RDS)))
  saveRDS(filtRs, file.path(path_RDS, paste0("filtRs_", name_RDS)))

  #end script message
  print("filter and trim completed")



  # ------------------------------------------------------------------------
  # output quality profiles
  # ------------------------------------------------------------------------

  print("quality profile plotting started")

  #path_pdf <- dirname(snakemake$output[[length(snakemake$output)]]) #get only the first part of the path
  path_pdf <- "13_output/quality_check" #path
  dir.create(path_pdf, recursive = T)

  #quality profiles unfiltered and untrimmed
  pdf(paste0(path_pdf, "/", taxa, "_", run, "_reads_quality_unfilt_untrim.pdf"))
  for (i in 1:length(sample.names)) {
    try(figure <- plotQualityProfile(c(fnFs[i],fnRs[i])))
    try(print(figure))
    try(rm(figure))
  }
  dev.off()

  #quality profiles filtered and trimmed
  pdf(paste0(path_pdf, "/", taxa, "_", run, "_reads_quality_filt_trim.pdf"))
  for (i in 1:length(sample.names)) {
    try(figure <- plotQualityProfile(c(filtFs[i],filtRs[i])))
    try(print(figure))
    try(rm(figure))
  }
  dev.off()

  #end script message
  print("quality profile plotting completed")


  #set sink back to std
  sink()


} else{
  #create empty output files needed for snakemake
  empty <- vector()

  path_RDS <- file.path(dirname(dirname(output_file)), "RDS") #get only the first part of the path
  name_RDS <- paste0(taxa, "_", run, ".RDS")
  dir.create(path_RDS)
  saveRDS(empty, file.path(path_RDS, paste0("out_", name_RDS)))
  saveRDS(empty, file.path(path_RDS, paste0("filtFs_", name_RDS)))
  saveRDS(empty, file.path(path_RDS, paste0("filtRs_", name_RDS)))
}
