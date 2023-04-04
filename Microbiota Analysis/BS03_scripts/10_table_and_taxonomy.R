#!/usr/bin/env Rscript

# get command line arguments as an array
args <- commandArgs(trailingOnly = TRUE)

#create snakemake list from passed arguments

#args files
#1    input
#2-6  output
#7-8  log
#9-15 parameters

snakemake <- list(input=args[1], output=args[2:6], log=args[7:8], parameters=args[9:15])

#redirect std out and std err to log file
log_out <- file(snakemake$log[[1]], open="wt")
log_err <- file(snakemake$log[[2]], open="wt")
sink(log_out, type="output")
sink(log_err, type="message")

#libaries
library(dada2)
library(vegan)
library(Biostrings)
library(DECIPHER)

#get taxa
taxa <- basename(dirname(snakemake$input[[1]]))

#import RDS
mergers <- readRDS(snakemake$input[[1]])

print("table and taxonomy started")


########################################
# setwd("/Users/j.waechli/Desktop/01_dadapipe")
# taxa <- "bacteria"
# mergers <- readRDS("09_merged/RDS/bacteria/mergers.RDS")
# classifier <- "both"
# taxa_database_RDP <- "db/silva_nr_v132_train_set.fa.gz"
# taxa_database_IDTAXA <- load("db/SILVA_SSU_r138_2019.rdata")
# path_RDS <- "10_table_and_taxonomy/RDS/bacteria"
# 
# setwd("/Users/j.waechli/Desktop/01_dadapipe")
# taxa <- "fungi"
# mergers <- readRDS("09_merged/RDS/fungi/mergers.RDS")
# classifier <- "both"
# taxa_database_RDP <- "db/sh_general_release_dynamic_10.05.2021.fasta"
# taxa_database_IDTAXA <- load("db/UNITE_v2020_February2020.rdata")
# path_RDS <- "10_table_and_taxonomy/RDS/fungi"
# path_pdf <- "13_output/fungi/quality_check"
# minBoot=50

########################################

# ------------------------------------------------------------------------
# ASV table
# ------------------------------------------------------------------------

#create table and remove bimera
seqtab <- makeSequenceTable(mergers)
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# ------------------------------------------------------------------------
# Assign taxanomy
# ------------------------------------------------------------------------
#choose classifier
classifier <- snakemake$parameters[1]

#get the database for bacteria or fungi
if(taxa == "bacteria" & classifier %in% c("RDP", "both")) {taxa_database_RDP <-  snakemake$parameters[4]}
if(taxa == "bacteria" & classifier %in% c("IDTAXA", "both")) {load(snakemake$parameters[5])}
if(taxa == "fungi" & classifier %in% c("RDP", "both")) {taxa_database_RDP <-  snakemake$parameters[6]}
if(taxa == "fungi" & classifier %in% c("IDTAXA", "both")) {load(snakemake$parameters[7])}

#RDP classification
if(classifier %in% c("RDP", "both")){
  
  minBoot <- as.numeric(snakemake$parameters[2])
  TAXA_RDP <- assignTaxonomy(seqtab_nochim, taxa_database_RDP, multithread=TRUE, minBoot = minBoot)
  TAXA_RDP <- as.data.frame(TAXA_RDP)
  
  #cleaner taxa for fungi
  if(taxa == "fungi"){for (i in 1:ncol(TAXA_RDP)) {TAXA_RDP[,i] <- gsub(".__", "", TAXA_RDP[,i])}}
  
  if(classifier=="RDP"){TAXA <- TAXA_RDP} #rename if only RDP
}

#IDTAXA classification
if(classifier %in% c("IDTAXA", "both")){
  #assign
  dna <- DNAStringSet(getSequences(seqtab_nochim)) # Create a DNAStringSet from the ASVs
  minBoot <- as.numeric(snakemake$parameters[3])
  ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE, threshold = minBoot) # use all processors
  if(taxa == "bacteria"){ranks <- c("domain", "phylum", "class", "order", "family", "genus")} # ranks of interest
  if(taxa == "fungi"){ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")} # ranks of interest
  
  # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
  TAXA_IDTAXA <- t(sapply(ids, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
  }))
  
  TAXA_IDTAXA <- as.data.frame(TAXA_IDTAXA)
  rownames(TAXA_IDTAXA) <- getSequences(seqtab_nochim)
  colnames(TAXA_IDTAXA) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  
  #overwrite unidentified taxa and taxa whre only species is assigned with NA
  unidentified <- apply(TAXA_IDTAXA, 2, function(x) grepl("unidentified",x))
  TAXA_IDTAXA[unidentified] <- NA
  only_sp <- is.na(TAXA_IDTAXA$Genus) & !(is.na(TAXA_IDTAXA$Species))
  try(TAXA_IDTAXA$Species[only_sp] <- NA, silent = T)
  
  if(classifier=="IDTAXA"){TAXA <- TAXA_IDTAXA}   #rename if only IDTAXA
}

#combine classifications
if(classifier=="both"){

  #empty df
  TAXA <- TAXA_RDP
  TAXA[,] <- NA
  TAXA <- as.data.frame(TAXA)
  TAXA$classifier <- NA
  
  #choose classifier with more assigned TAXA-levels
  for(i in 1:nrow(TAXA)){
    n_RDP <- sum(!(is.na(TAXA_RDP[i,])))
    n_IDTAXA <- sum(!(is.na(TAXA_IDTAXA[i,])))
    if(n_RDP >= n_IDTAXA){tax <- TAXA_RDP[i,]; classifier <- "RDP"} #same or more assigned TAXA-levels in RDP
    else{tax <- TAXA_IDTAXA[i,]; classifier <- "IDTAXA"} #more assigned TAXA-levels in IDTAXA
    TAXA[i,1:length(tax)] <- tax
    TAXA$classifier[i] <- classifier

  }
}

# ------------------------------------------------------------------------
# Header names
# ------------------------------------------------------------------------

#change headers from sequence to ASV numbers
seq <- colnames(seqtab_nochim)
ASV <- c()
for (i in 1:length(seq)) {ASV[i] <- paste("ASV",i , sep="") }
ASV_seq <- data.frame(ASV=ASV, seq=seq)
ASV_seq[,1] <- paste(">", ASV_seq[,1], sep="")

#replace sequence headers by ASV numbers in seqtab_nochim
colnames(seqtab_nochim) <- ASV

#replace the sequence headers by the correct ASV number in TAXA
for (i in 1:nrow(TAXA)){
  for (j in 1:nrow(ASV_seq)){
    if(rownames(TAXA)[i] == ASV_seq[j,2]){
      rownames(TAXA)[i] <- ASV_seq[j,1]
    }
  }
}
rownames(TAXA) <- gsub(">", "", rownames(TAXA))

# ------------------------------------------------------------------------
# output RDS files
# ------------------------------------------------------------------------

path_RDS <- dirname(snakemake$output[[1]])
dir.create(path_RDS, recursive = T)
saveRDS(seqtab, file.path(path_RDS, "seqtab.RDS"))
saveRDS(seqtab_nochim, file.path(path_RDS, "seqtab_nochim.RDS"))
saveRDS(TAXA, file.path(path_RDS, "taxa.RDS"))
saveRDS(ASV_seq, file.path(path_RDS, "sequences.RDS"))

# ------------------------------------------------------------------------
# output rarefaction plot
# ------------------------------------------------------------------------

path_pdf <- dirname(snakemake$output[[length(snakemake$output)]]) #get only the first part of the path
dir.create(path_pdf, recursive = T)
runs <- unique(sapply(strsplit(rownames(seqtab_nochim), "_", fixed=F), `[`, 1))

pdf(paste0(path_pdf, "/", taxa, "_rarefaction.pdf"))
for (run in runs) {
  DAT <- seqtab_nochim[grep(run, rownames(seqtab_nochim)),]
  vegan::rarecurve(DAT, step = 20, label=F, main = paste("Rarefaction", run), ylab="ASVs", xlab="Sequencing depth")
}
dev.off()

# path_pdf <- dirname(snakemake$output[[length(snakemake$output)]]) #get only the first part of the path
# dir.create(path_pdf, recursive = T)
# 
# pdf(paste0(path_pdf, "/", taxa, "_rarefaction.pdf"))
#   vegan::rarecurve(seqtab_nochim, step = 20, label=F, main = "Rarefaction Plot", ylab="ASVs", xlab="Sequencing depth")
# dev.off()

#end script message
print("table and taxonomy completed")

#set sink back to std
sink()
