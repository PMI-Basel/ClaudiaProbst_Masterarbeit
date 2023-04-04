#!/usr/bin/env Rscript

# get command line arguments as an array
args <- commandArgs(trailingOnly = TRUE)

#create snakemake list from passed arguments

#args files
#1-n     input
#n+1-n+2 output
#n+3-n+4 log

n1 <- sum(grepl("07_filter_and_trim", args, fixed=T))
n2 <- sum(grepl("08_error_learning", args, fixed=T))
n3 <- sum(grepl("09_merged", args, fixed=T))
n4 <- sum(grepl("10_table_and_taxonomy", args, fixed=T))
n <- sum(n1, n2, n3, n4)

snakemake <- list(input=args[1:n], output=args[(n+1):(n+2)], log=args[(n+3):(n+4)])


#redirect std out and std err to log file
log_out <- file(snakemake$log[[1]], open="wt")
log_err <- file(snakemake$log[[2]], open="wt")
sink(log_out, type="output")
sink(log_err, type="message")

#libaries
library(gtools)
library(dada2)
library(xlsx)

print("sequence tracking samples started")

#get all possible input files as vector
files <- unique(list.files(dirname(unlist(snakemake$input)), full.names = T))

#filter by files need for track table
files <- files[grepl("out_", files, fixed=T) | grepl("dadaFs.RDS", files, fixed=T) | grepl("dadaRs.RDS", files, fixed=T) |
               grepl("mergers.RDS", files, fixed=T) | grepl("seqtab_nochim.RDS", files, fixed=T)]


#get taxa and runs
taxa <- sapply(strsplit(basename(snakemake$input[[1]]), "_", fixed=F), `[`, 2)

#filter files by taxa
files_taxa <- files[grepl(taxa, files, fixed=T)]


#import files
for (f in files_taxa) {
  temp <- readRDS(f)
  name <- gsub(".RDS", "", basename(f))
  assign(name, temp)
  rm(name, temp, f)
}

#get runs
runs <- unique(sapply(strsplit(basename(files_taxa[grepl("out_", files_taxa, fixed=T)]), "_", fixed=F), `[`, 3))
runs <- gsub(".RDS", "", runs)

#combine out files by runs
out <- get(paste("out", taxa, runs[1], sep="_"))
if(length(runs)>1){for(i in runs[2:length(runs)]){out <- rbind(out, get(paste("out", taxa, i, sep="_")))}}


##########################
# setwd("/Users/j.waechli/Desktop/01_dadapipe")
# out <- readRDS("07_filter_and_trim/RDS/out_bacteria_BS02.RDS")
# dadaFs <- readRDS("08_error_learning/RDS/bacteria/dadaFs_BS02.RDS")
# dadaRs <- readRDS("08_error_learning/RDS/bacteria/dadaRs_BS02.RDS")
# mergers <- readRDS("09_merged/RDS/bacteria/mergers.RDS")
# seqtab_nochim <- readRDS("10_table_and_taxonomy/RDS/bacteria/seqtab_nochim.RDS")
#path_RDS <- "12_sequence_tracking/RDS"
#path_xlsx <- "13_output/sequence_tracking"
#taxa <- "bacteria"
##########################



# ------------------------------------------------------------------------
# track table
# ------------------------------------------------------------------------

#create track table
getN <- function(x) sum(getUniques(x))

rownames(out) <- sapply(strsplit(rownames(out),"_F.fastq"), `[`, 1)
n_dadaFs <- sapply(dadaFs, getN); #names(n_dadaFs) <- sapply(strsplit(names(n_dadaFs),"_F_filt.fastq.gz"), `[`, 1)
n_dadaRs <- sapply(dadaRs, getN); #names(n_dadaRs) <- sapply(strsplit(names(n_dadaRs),"_R_filt.fastq.gz"), `[`, 1)
n_mergers <- sapply(mergers, getN); names(n_mergers) <- sapply(strsplit(names(n_mergers),"_F_filt.fastq.gz"), `[`, 1)
n_seqtab_nochim <- rowSums(seqtab_nochim); names(n_seqtab_nochim) <- sapply(strsplit(names(n_seqtab_nochim),"_F_filt.fastq.gz"), `[`, 1)

track <- t(gtools::smartbind(t(out), t(n_dadaFs), t(n_dadaRs), t(n_mergers), t(n_seqtab_nochim),fill=0))
colnames(track) <- c("input", "filter_length", "filter_maxN", "filter_maxEE", "denoisedF", "denoisedR", "merged", "nonchim")


# ------------------------------------------------------------------------
# output RDS files
# ------------------------------------------------------------------------

#first path for RDS
path_RDS <- dirname(snakemake$output[[1]])
dir.create(path_RDS, recursive = T)
saveRDS(track, file.path(path_RDS, paste0(taxa, "_sample_track.RDS")))

#second path for xlsx
path_xlsx <- dirname(snakemake$output[[length(snakemake$output)]])
dir.create(path_xlsx, recursive = T)
write.xlsx(track, file.path(path_xlsx, paste0(taxa, "_sample_track.xlsx")))

#end script message
print("sequence tracking samples completed")

#set sink back to std
sink()
