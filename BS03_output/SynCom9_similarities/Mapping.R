library("Biostrings")
library("DECIPHER")
library("xlsx")
library("tidyr")
library("tibble")
library("ggplot2")
library("ggrepel")
library("scales")


## set source to file location
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

#import
ASV <- Biostrings::readDNAStringSet("bacteria_SEQ100.fasta")
ISO <- Biostrings::readDNAStringSet("16S_syn9.fasta")

### MAPPING ###

#convert DNAStringSet as data frame
ASV_df <- as.data.frame(ASV)
colnames(ASV_df) <- "seqs"

ISO_df <- as.data.frame(ISO)
colnames(ISO_df) <- "seqs"
rownames(ISO_df) <- sapply(strsplit(names(ISO),"-"), `[`, 1)

#create new data.frame to save distances
distances <- data.frame(matrix(NA, nrow = nrow(ASV_df), ncol = nrow(ISO_df)), row.names = rownames(ASV_df))
colnames(distances) <- rownames(ISO_df)

#mapping function
map <- function(ASV_seq) {
  seqs_to_map <- data.frame(seqs=c(as.character(ASV_seq),as.character(ISO_df$seqs[i]))) #create data.frame with ASV and ISO to map
  seqs_to_map <- Biostrings::DNAStringSet(seqs_to_map$seqs) #convert to DNA string set
  seqs_aln <- DECIPHER::AlignSeqs(seqs_to_map, processors = 1, verbose = F) #align, processors = NULL uses the max available number of processors
  dist_all <- DECIPHER::DistanceMatrix(seqs_aln, processors = 1, verbose = F) #similarity table
  return(dist_all[1,2])
}

#calculate distances
#this loop may take a while
for(i in 1:nrow(ISO_df)) {
  #loop over each ASV
  distances[,i] <- apply(ASV_df, 1, map) #mapping and save values in dist data.frame
  print(paste0("seq ", i, "/",  nrow(ISO_df), " mapped"))
}

#create similarity- and top-hits table
similarities <- 1 - distances #similarity
hits <- apply(similarities, 2, function(x) {c(rownames(similarities)[which(x==max(x))][1], max(x))}) #most similar ASV and the similarity
hits <- t(as.data.frame(hits)) #convert to data.frame
hits[,2] <- round(as.numeric(hits[,2]), 3) * 100 #percent identical
colnames(hits) <- c("ASV", "similarity")

similarities_long <- add_column(similarities, ASV=rownames(similarities), .before = 1)
similarities_long <- pivot_longer(similarities_long, -"ASV", names_to = "ISO", values_to = "values")


#plot
plot1 <- ggplot(similarities_long[similarities_long$values>=0.95,], aes(x=ISO,y=values,color=values))+
            geom_point(size=3, alpha=0.5)+
            geom_label_repel(aes(label = ASV), na.rm = TRUE , show.legend = F)+
            scale_color_gradient(low = "blue", high = "red")+
            scale_y_continuous(labels = percent)+
            theme_bw()+
            ggtitle("Similarities")

plot1
ggsave("Similarities.png", height=12, width=11) 

# ggplot(similarities_long[], aes(x=ASV,y=values,color=values))+
#   geom_point(size=3, alpha=0.5)+
#   geom_label_repel(aes(label = ISO), na.rm = TRUE , show.legend = F)+
#   scale_color_gradient(low = "blue", high = "red")+
#   scale_y_continuous(labels = percent)+
#   theme_bw()+
#   ggtitle("Similarities")


#write.files
dir.create("output")
write.xlsx(similarities, "output/similarities.xlsx")
png("output/Top-Similarities.png",8000, 6000, res=600)
  plot1
dev.off()
