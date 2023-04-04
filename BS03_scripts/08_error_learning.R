#!/usr/bin/env Rscript

# get command line arguments as an array
args <- commandArgs(trailingOnly = TRUE)


#create snakemake list from passed arguments

#args files
#1-2  input
#3-8  output
#9-10  log

snakemake <- list(input=args[1:2], output=args[3:8], log=args[9:10])


#redirect std out and std err to log file
log_out <- file(snakemake$log[[1]], open="wt")
log_err <- file(snakemake$log[[2]], open="wt")
sink(log_out, type="output")
sink(log_err, type="message")

#libaries
library(dada2)
library(grid)
library(gridExtra)

#import RDS
filtFs <- readRDS(snakemake$input[[1]])
filtRs <- readRDS(snakemake$input[[2]])

#get taxa and run
taxa <- sapply(strsplit(basename(snakemake$input[[1]]), "_", fixed=F), `[`, 2)
run <- sapply(strsplit(basename(snakemake$input[[1]]), "_", fixed=F), `[`, 3)
run <- gsub(".RDS", "", run)

#get all taxa-run combinations
DESIGN <- read.delim("01_design/design.tab", sep = "\t")
taxa_runs <- unique(paste(DESIGN$taxa, DESIGN$run_nbr, sep="_"))

#check if taxa-run combination exist, if not skip everything
taxa_run <- paste(taxa, run, sep="_")
if (taxa_run %in% taxa_runs) {
  

  # ------------------------------------------------------------------------
  # Learn the error rate
  # ------------------------------------------------------------------------
  
  print("error learning started")
  
  #error rate
  errF <- learnErrors(filtFs, multithread=TRUE)
  errR <- learnErrors(filtRs, multithread=TRUE)
  
  # ------------------------------------------------------------------------
  # Dereplicate the files
  # ------------------------------------------------------------------------
  
  #dereplicate
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  
  # ------------------------------------------------------------------------
  # Sample inference algorithm
  # ------------------------------------------------------------------------
  
  # inference algortihm
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
  dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
  
  # ------------------------------------------------------------------------
  # simplify names
  # ------------------------------------------------------------------------
  
  names(dadaFs) <- sapply(strsplit(names(dadaFs),"_F_filt.fastq.gz"), `[`, 1)
  names(dadaRs) <- sapply(strsplit(names(dadaRs),"_R_filt.fastq.gz"), `[`, 1)
  names(derepFs) <- sapply(strsplit(names(derepFs),"_F_filt.fastq.gz"), `[`, 1)
  names(derepRs) <- sapply(strsplit(names(derepRs),"_R_filt.fastq.gz"), `[`, 1)
  
  # ------------------------------------------------------------------------
  # output RDS files
  # ------------------------------------------------------------------------
  
  path_RDS <- file.path(dirname(dirname(snakemake$output[[1]])), taxa)
  name_RDS <- paste0(run, ".RDS")
  dir.create(path_RDS, recursive = T)
  saveRDS(errF, file.path(path_RDS, paste0("errF_", name_RDS)))
  saveRDS(errR, file.path(path_RDS, paste0("errR_", name_RDS)))
  saveRDS(derepFs, file.path(path_RDS, paste0("derepFs_", name_RDS)))
  saveRDS(derepRs, file.path(path_RDS, paste0("derepRs_", name_RDS)))
  saveRDS(dadaFs, file.path(path_RDS, paste0("dadaFs_", name_RDS)))
  saveRDS(dadaRs, file.path(path_RDS, paste0("dadaRs_", name_RDS)))
  
  # ------------------------------------------------------------------------
  # output error rate learing
  # ------------------------------------------------------------------------
  
  #path_pdf <- dirname(snakemake$output[[length(snakemake$output)]]) #get only the first part of the path
  path_pdf <- "13_output/quality_check" #path
  
  pdf(paste0(path_pdf, "/", taxa, "_", run, "_error_rate.pdf"))
    plot(plotErrors(errF, nominalQ=TRUE))
    grid.text(paste(run, "forward"),hjust=-1.2, vjust = -27.5, rot = 90)
    plot(plotErrors(errR, nominalQ=TRUE))
    grid.text(paste(run, "reverse"),hjust=-1.2, vjust = -27.5, rot = 90)
  dev.off()
  
  
  #end script message
  print("error learning completed")
  
  #set sink back to std
  sink()
  
} else{
  #create empty output files needed for snakemake
  empty <- vector()
  
  path_RDS <- file.path(dirname(dirname(snakemake$output[[1]])), taxa)
  name_RDS <- paste0(run, ".RDS")
  dir.create(path_RDS, recursive = T)
  saveRDS(empty, file.path(path_RDS, paste0("errF_", name_RDS)))
  saveRDS(empty, file.path(path_RDS, paste0("errR_", name_RDS)))
  saveRDS(empty, file.path(path_RDS, paste0("derepFs_", name_RDS)))
  saveRDS(empty, file.path(path_RDS, paste0("derepRs_", name_RDS)))
  saveRDS(empty, file.path(path_RDS, paste0("dadaFs_", name_RDS)))
  saveRDS(empty, file.path(path_RDS, paste0("dadaRs_", name_RDS)))
}

