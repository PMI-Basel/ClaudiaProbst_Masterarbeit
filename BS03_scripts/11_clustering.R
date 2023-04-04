#!/usr/bin/env Rscript

# get command line arguments as an array
args <- commandArgs(trailingOnly = TRUE)

#create snakemake list from passed arguments

#args files
#1-3  input
#4-8  output
#9-10 log
#11   parameters

snakemake <- list(input=args[1:3], output=args[4:8], log=args[9:10], parameters=args[11])


#redirect std out and std err to log file
log_out <- file(snakemake$log[[1]], open="wt")
log_err <- file(snakemake$log[[2]], open="wt")
sink(log_out, type="output")
sink(log_err, type="message")

#libaries
library(Biostrings)
library(DECIPHER)
library(dplyr)
library(seqRFLP)

#get taxa
taxa <- basename(dirname(snakemake$input[[1]]))

#import RDS
seqtab_nochim <- readRDS(snakemake$input[[1]])
TAXA <- readRDS(snakemake$input[[2]])
ASV_seq <- readRDS(snakemake$input[[3]])

print("clustering started")

###############
# setwd("/Users/j.waechli/Desktop/Jan_Waelchli/01_Projects/01_Support/Charlotte/BS01/Version_6/10_table_and_taxonomy/RDS/bacteria")
# seqtab_nochim <- readRDS("seqtab_nochim.RDS")
# TAXA <- readRDS("taxa.RDS")
# ASV_seq <- readRDS("sequences.RDS")
# taxa <- "bacteria"
# nproc <- 4
# percent <- 97
###############


# ------------------------------------------------------------------------
# Clustering
# ------------------------------------------------------------------------

nproc <- as.numeric(snakemake$parameters[1])

#change ASV names to corresponding sequences
seqtab_nochim_ASV <- seqtab_nochim
colnames(seqtab_nochim_ASV) <- ASV_seq$seq

#align sequences and compute distance matrix
ASV_seqs <- Biostrings::DNAStringSet(ASV_seq$seq) #convert to DNA string set
#ASV_seqs_aln <- DECIPHER::AlignSeqs(ASV_seqs, processors = nproc)
#ASV_dist <- DECIPHER::DistanceMatrix(ASV_seqs_aln, processors = nproc)

#get cutoff
percent <- as.numeric(basename(dirname(snakemake$output[[1]])))
cutoff <- (100-percent)/100

#########################
#save(ASV_seqs, nproc, cutoff, file = paste0("clustering_", percent))
#########################

#OTUs
#DECIPHER::Clusterize() was DECIPHER::IdClusters() in DECIPHER versions < 2.24.0
ASV_cluster <- DECIPHER::Clusterize(ASV_seqs, processors = nproc, cutoff = cutoff) #cutoff = 0.03 -> 97% OTU

#get the data
ASV_cluster$ASV <- gsub(">","",ASV_seq$ASV) #add ASV name
ASV_cluster$ASV_seq <- ASV_seq$seq #add sequences

# Number of sequences per ASV
ASV_cluster$ASV_abu <- NA
for (seq in ASV_cluster$ASV_seq) {
  abu <- sum(seqtab_nochim_ASV[,colnames(seqtab_nochim_ASV) == seq])
  ASV_cluster$ASV_abu[ASV_cluster$ASV_seq == seq] <- abu
}

#most abundant sequences for each OTU
ASV_cluster$OTU_seq <- NA
ASV_cluster$OTU_abu <- NA

for (OTU in ASV_cluster$cluster) {
  ASV_cluster_subset <- ASV_cluster[ASV_cluster$cluster == OTU,] #all entries for current OTU
  n_seq <- sum(ASV_cluster_subset$ASV_abu)
  most_abu_seq <- ASV_cluster_subset$ASV_seq[ASV_cluster_subset$ASV_abu == max(ASV_cluster_subset$ASV_abu)]# OTU(s) with most seqs
  most_abu_seq <- most_abu_seq[1] #take the first entry in case multiple OTUs have max number of seqs
  ASV_cluster$OTU_seq[ASV_cluster$cluster == OTU] <- most_abu_seq #add seq to table
  ASV_cluster$OTU_abu[ASV_cluster$cluster == OTU] <- n_seq #add num to table
}

#sort by OTU abundance
ASV_cluster_sort <- ASV_cluster[order(ASV_cluster$OTU_abu, decreasing = T),] #sort by OTU name

#rename OTU by abundance
clusters_sort <- unique(ASV_cluster_sort$cluster)
ASV_cluster$OTU <- NA
for(i in seq(length(clusters_sort))){
  ASV_cluster$OTU[ASV_cluster$cluster == clusters_sort[i]] <- paste0("OTU", i) #rename
}

### FINAL OTU FILES

#OTU ASV_seq
OTU_ASV_seq_all <- ASV_cluster[order(ASV_cluster$OTU_abu, decreasing = T),]
OTU_ASV_seq_all <- data.frame(OTU = OTU_ASV_seq_all$OTU, seq = OTU_ASV_seq_all$OTU_seq)
OTU_ASV_seq <- data.frame(OTU = unique(OTU_ASV_seq_all$OTU), seq = NA)
for (OTU in OTU_ASV_seq$OTU) {
  seq <- OTU_ASV_seq_all$seq[OTU_ASV_seq_all$OTU == OTU][1] #get the sequence
  OTU_ASV_seq$seq[OTU_ASV_seq$OTU == OTU] <- seq #add to seqtab_nochima frame
}


#OTU seqtab_nochim
OTU_seqtab_nochim <- seqtab_nochim_ASV %>% t %>% rowsum(ASV_cluster$cluster) %>% t #cluster
for(i in seq(length(clusters_sort))){
  colnames(OTU_seqtab_nochim)[colnames(OTU_seqtab_nochim) == clusters_sort[i]] <- paste0("OTU", i) #rename
}
OTU_seqtab_nochim <- OTU_seqtab_nochim[,OTU_ASV_seq$OTU] #sort


#OTU TAXA
OTU_TAXA <- data.frame(matrix(nrow = nrow(OTU_ASV_seq), ncol = ncol(TAXA)))
rownames(OTU_TAXA) <- OTU_ASV_seq$OTU
colnames(OTU_TAXA) <- colnames(TAXA)

for(OTU in OTU_ASV_seq$OTU){
  seq <- OTU_ASV_seq$seq[OTU_ASV_seq$OTU == OTU] #OTU seq
  ASV <- ASV_cluster$ASV[ASV_cluster$ASV_seq == seq] #ASV belong to OTU seq
  
  TAXA_ASV <- TAXA[rownames(TAXA) == ASV[1], ] #TAXA belong to the first ASV
  OTU_TAXA[rownames(OTU_TAXA) == OTU, ] <- TAXA_ASV #transfer to OTU TAXA
}

#OTU_ASV
OTU_ASV <- data.frame(OTU = ASV_cluster$OTU, ASV = ASV_cluster$ASV, 
                      OTU_seq = ASV_cluster$OTU_seq, ASV_seq = ASV_cluster$ASV_seq,
                      OTU_abu = ASV_cluster$OTU_abu, ASV_abu = ASV_cluster$ASV_abu)

OTU_ASV <- OTU_ASV[order(ASV_cluster$OTU_abu, decreasing = T),] #sort

# ------------------------------------------------------------------------
# output RDS files
# ------------------------------------------------------------------------

path_RDS <- dirname(snakemake$output[[1]])
if(!(dir.exists(path_RDS))){dir.create(path_RDS, recursive = T)} #create folder if it does not exist

saveRDS(OTU_seqtab_nochim, file.path(path_RDS, paste0("OTU_seqtab_nochim_", percent, ".RDS")))
saveRDS(OTU_ASV_seq, file.path(path_RDS, paste0("OTU_ASV_seq_", percent, ".RDS")))
saveRDS(OTU_TAXA, file.path(path_RDS, paste0("OTU_taxa_", percent, ".RDS")))
saveRDS(OTU_ASV, file.path(path_RDS, paste0("OTU_ASV_", percent, ".RDS")))

# ------------------------------------------------------------------------
# output DAT, TAXA, SEQ
# ------------------------------------------------------------------------

path_pdf <- dirname(snakemake$output[[length(snakemake$output)]]) #get path
if(!(dir.exists(path_pdf))){dir.create(path_pdf, recursive = T)} #create folder if it does not exist

write.table(OTU_seqtab_nochim, paste0(path_pdf, "/", taxa, "_COUNT", percent, ".tab"), sep="\t")
write.table(OTU_TAXA, paste0(path_pdf, "/", taxa, "_TAXA", percent, ".tab"), sep="\t")
#write.table(OTU_ASV_seq, paste0(path_pdf, "/", taxa, "_SEQ", percent, ".tab"), sep="\t")
invisible(seqRFLP::dataframe2fas(OTU_ASV_seq, paste0(path_pdf, "/", taxa, "_SEQ", percent, ".fasta")))
if (percent < 100) {write.table(OTU_ASV, paste0(path_pdf, "/", taxa, "_ABU", percent, ".tab"), sep="\t", row.names = F)}

#end script message
print("clustering completed")

#set sink back to std
sink()
