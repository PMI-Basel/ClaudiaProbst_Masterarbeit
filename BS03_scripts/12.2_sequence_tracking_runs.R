#!/usr/bin/env Rscript

# get command line arguments as an array
args <- commandArgs(trailingOnly = TRUE)

#create snakemake list from passed arguments

#args files
#1-n     input
#n+1-n+3 output
#n+4-n+5 log

n1 <- sum(grepl("02_run_qc", args, fixed=T))
n2 <- sum(grepl("12_sequence_tracking", args, fixed=T)) -3 #-3 output/log-files
n <- sum(n1, n2)

snakemake <- list(input=args[1:n], output=args[(n+1):(n+3)], log=args[(n+4):(n+5)])


#redirect std out and std err to log file
log_out <- file(snakemake$log[[1]], open="wt")
log_err <- file(snakemake$log[[2]], open="wt")
sink(log_out, type="output")
sink(log_err, type="message")

#libaries
library(tidyverse)
library(xlsx)

print("sequence tracking runs started")

#get taxa and runs
DESIGN <- read.delim("01_design/design.tab", sep = "\t")
taxa <- unique(DESIGN$taxa)
taxa <- taxa[nchar(taxa)>0]

DESIGN <- read.delim("01_design/design.tab", sep = "\t")
run <- unique(DESIGN$run_nbr)
run <- run[nchar(run)>0]

input <- list(paste0("02_run_qc/", run, ".r1_sequences.txt"), paste0("02_run_qc/", run, ".r2_sequences.txt"), 
              paste0("12_sequence_tracking/RDS/", taxa, "_sample_track.RDS"))

#all input files
files <- unlist(input)

#import files
for (f in files) {
  if (grepl("RDS", f, fixed=T)) {
    temp <- readRDS(f)
    name <- gsub(".RDS", "", basename(f))
  }
  if (grepl("txt", f, fixed=T)) {
    temp <- scan(f)
    name <- gsub("_sequences.txt", "", basename(f))
  }
  assign(name, temp)
  rm(name, temp, f)
}

#taxa and run data.frame
runs_bac <- vector()
runs_fun <- vector()
try(runs_bac <- unique(sapply(strsplit(rownames(bacteria_sample_track), "_", fixed=F), `[`, 1)), silent=T)
try(runs_fun <- unique(sapply(strsplit(rownames(fungi_sample_track), "_", fixed=F), `[`, 1)), silent=T)

runs_all <- unique(c(runs_bac, runs_fun))
runs <- data.frame(runs=runs_all, bac=F, fun=F)
try(runs$bac <- runs_bac %in% runs$runs, silent = T)
try(runs$fun <- runs_fun %in% runs$runs, silent = T)

# ------------------------------------------------------------------------
# Loop over all runs
# ------------------------------------------------------------------------

#track summary table
for (run in runs$runs) {
  
  # ------------------------------------------------------------------------
  # sequences tracking summary table
  # ------------------------------------------------------------------------
  
  #summary of read tracking
  run_seqs <- get(paste0(run, ".r1"))
  #run_dem_seqs_bac <- seqs[,5][seqs$V1 == run]
  #run_dem_seqs_fun <- seqs[,7][seqs$V1 == run]
  ifelse(runs$bac[runs$runs==run], 
         run_track_bac <- colSums(bacteria_sample_track[grepl(run,rownames(bacteria_sample_track)),]), 
         run_track_bac <- rep(NA, 8))
  ifelse(runs$fun[runs$runs==run], 
         run_track_fun <- colSums(fungi_sample_track[grepl(run,rownames(fungi_sample_track)),]), 
         run_track_fun <- rep(NA, 8))
  
  ifelse("bacteria" %in% taxa, names <- names(run_track_bac), names <- names(run_track_fun))
  
  run_track_sum <- data.frame(run=run,
                              step=c("raw", names[1:4], "denoised", names[7:8]),
                              bacteria=c(NA, run_track_bac[1:4], mean(run_track_bac[5:6]), run_track_bac[7:8]),
                              fungi=c(NA, run_track_fun[1:4], mean(run_track_fun[5:6]), run_track_fun[7:8])) #mean of denoising
  
  run_track_sum$step <- factor(run_track_sum$step, levels = run_track_sum$step)
  run_track_sum$total <- rowSums(run_track_sum[,3:4], na.rm = T)
  run_track_sum$total[1] <- run_seqs #add raw reads
  if(is.na(run_track_sum$bacteria[2])){run_track_sum$fungi[1] <- run_track_sum$total[1]} #if only fungi we do now the total of fungal raw reads
  if(is.na(run_track_sum$fungi[2])){run_track_sum$bacteria[1] <- run_track_sum$total[1]} #if only bacteria we do now the total of bacterial raw reads
  run_track_sum$percent <- paste(round(run_track_sum$total / run_seqs * 100, 1), "%")
  
  # ------------------------------------------------------------------------
  # Sequences tracking plot
  # ------------------------------------------------------------------------
  
  #long format for plot
  run_track_sum_long <- pivot_longer(run_track_sum, c(total, bacteria, fungi), names_to = "taxa", values_to = "seqs")
  run_track_sum_long$percent[run_track_sum_long$taxa != "total"] <- NA
  run_track_sum_long$taxa <- factor(run_track_sum_long$taxa, levels = c("total", "bacteria", "fungi"))
  
  #plot
  run_track_plot <- ggplot(run_track_sum_long,aes(step, seqs, group=taxa))+
    geom_point(position=position_dodge(0.2))+
    geom_line(aes(col=taxa),position=position_dodge(0.2))+
    #geom_line(aes(step, bacteria), col="red")+
    geom_text(aes(label=percent), hjust=0.2, vjust=-.5)+
    scale_y_continuous(trans="log2",limits = c(min(run_track_sum$fungi)*0.5,max(run_track_sum$total)*1.05))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90))+
    ggtitle(paste("Sequences tracking:", run))
  
  
  #store for each run
  assign(paste0(run, "_track_sum"), run_track_sum)
  assign(paste0(run, "_track_plot"), run_track_plot)
  
}


# ------------------------------------------------------------------------
# output RDS files
# ------------------------------------------------------------------------

#combine runs
#tables
for (run in runs$runs) {
  if(run == runs$runs[1]){track_sum <- get(paste0(run,"_track_sum"))}
  else{track_sum <- rbind(track_sum, get(paste0(run,"_track_sum")))}
}

saveRDS(track_sum, snakemake$output[[1]])
write.xlsx(track_sum, snakemake$output[[2]])



#plots
pdf(snakemake$output[[3]])
for (run in runs$runs) {
  plot(get(paste0(run, "_track_plot")))
}
dev.off()


#end script message
print("sequence tracking runs completed")

#set sink back to std
sink()
