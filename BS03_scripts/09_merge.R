#!/usr/bin/env Rscript

# get command line arguments as an array
args <- commandArgs(trailingOnly = TRUE)

#create snakemake list from passed arguments

#args files
#1-n  input
#n+1-n+5 output
#n+6-n+7  log

n <- sum(grepl("08_error_learning", args, fixed=T)) #nbr of input files
snakemake <- list(input=args[1:n], output=args[(n+1):(n+5)], log=args[(n+6):(n+7)])

#redirect std out and std err to log file
log_out <- file(snakemake$log[[1]], open="wt")
log_err <- file(snakemake$log[[2]], open="wt")
sink(log_out, type="output")
sink(log_err, type="message")

#libaries
library(dada2)

###############################
#setwd("/Users/j.waechli/Desktop/01_dadapipe")
#files <- list.files(dirname("08_error_learning/RDS/fungi/dadaFs_BE12.RDS"), full.names = T)
###############################

#get the path from the first input file until the taxa and import all RDS
files <- list.files(dirname(snakemake$input[[1]]), full.names = T)
for (f in files) {
  temp <- readRDS(f)
  name <- gsub(".RDS", "", basename(f))
  if(length(temp)==0){files <- files[!(grepl(f,files,fixed=T))]} #remove empty files from file-vector
  if(length(temp)>0){assign(name, temp)} #only load file if file is not empty
}

#get taxa and runs
runs <- sapply(strsplit(basename(files), "_", fixed=F), `[`, 2)
runs <- gsub(".RDS", "", runs)
runs <- unique(runs)

print("merging started")

# ------------------------------------------------------------------------
# combine all bacterial or fungal samples from all runs
# ------------------------------------------------------------------------

derepFs <- vector()
for(i in runs){derepFs <- c(derepFs, get(paste0("derepFs_", i)))}
derepRs <- vector()
for(i in runs){derepRs <- c(derepRs, get(paste0("derepRs_", i)))}

dadaFs <- vector()
for(i in runs){dadaFs <- c(dadaFs, get(paste0("dadaFs_", i)))}
dadaRs <- vector()
for(i in runs){dadaRs <- c(dadaRs, get(paste0("dadaRs_", i)))}

# ------------------------------------------------------------------------
# Merge the forward and reverse reads
# ------------------------------------------------------------------------

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

###############
print("flag3")
###############

# ------------------------------------------------------------------------
# output RDS files
# ------------------------------------------------------------------------

path_RDS <- dirname(snakemake$output[[1]])
dir.create(path_RDS, recursive = T)
saveRDS(derepFs, file.path(path_RDS, "derepFs.RDS"))
saveRDS(derepRs, file.path(path_RDS, "derepRs.RDS"))
saveRDS(dadaFs, file.path(path_RDS, "dadaFs.RDS"))
saveRDS(dadaRs, file.path(path_RDS, "dadaRs.RDS"))
saveRDS(mergers, file.path(path_RDS, "mergers.RDS"))

#end script message
print("merging completed")

#set sink back to std
sink()
