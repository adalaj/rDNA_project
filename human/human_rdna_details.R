##Aim: Find the total nucleotide distribution percent in KY962518_2018 NCBI entry that is for human rDNA locus. 
## this calculation as been previously done in my earliar code. 
##however, after alot of back and forth we redefined the promoter and IGS boundary. 
## we simply took 2202 from IGS end side, because these were matching with UCSC hg38 nucleotides and added them to front of 5'ETS. We call them as promoter. 
## so manually i added 3500 bp form IGS to upstream of 5'ETS 

#in future if you want to avoid manual adding you can use stringr which i did for mouse rDNA refer code mouse_rdna_details.R





library(Biostrings)
library(seqinr)
library(stringr)
library(data.table)
library(tidyverse)


##I wanted to add new promoter and redefined IGS to rDNA
hg38_rDNA<- read.fasta(file = "KY962518_added_3500nt_IGS_upstream_nontemplate.fasta", forceDNAtolower = FALSE, as.string = TRUE) #belongs to package seqinr

getwd#forceDNAtolower false is used to keep dna seq in uppercase, as. string to entire seq as one string
hg38_rDNA_seq<- hg38_rDNA[[1]]
nchar(hg38_rDNA_seq)
#48338


#Now i want to separate the rdna1 into subsections described by the authors
#1..3657- 5’ external transcribed spacer
#3658..5526 - 18S ribosomal RNA
#5527..6596 - internal transcribed spacer 1
#6597..6753 - 5.8S ribosomal RNA
#6754..7920 - internal transcribed spacer 2
#7921..12971- 28S ribosomal RNA
#12972..13332 - external transcribed spacer
#I assume 13333 ..44838 - intergenic spacer; IGS


promoter<- str_sub(hg38_rDNA_seq, start= 1299, end = 3500) # my defined promoter
ets5<- str_sub(hg38_rDNA_seq, start = 3501, end = 7157)
s18<- str_sub(hg38_rDNA_seq, start= 7158, end = 9026)
its1<- str_sub(hg38_rDNA_seq, start= 9027, end = 10096)
s5.8 <- str_sub(hg38_rDNA_seq, start= 10097, end = 10253)
its2<- str_sub(hg38_rDNA_seq, start= 10254, end = 11420)
s28<- str_sub(hg38_rDNA_seq, start= 11421, end = 16471)
ets3<- str_sub(hg38_rDNA_seq, start= 16472, end = 16832)
igs<- str_sub(hg38_rDNA_seq, start= 16833, end = 46137) # my defined IGS 
entire_rdna<- str_sub(hg38_rDNA_seq, start= 1299, end=46136)
coding_rdna<- str_sub(hg38_rDNA_seq, start= 3501, end=16832)
entire_rdna_no_igs<- str_sub(hg38_rDNA_seq, start= 1299, end=16832)


Name<- c("Promoter_KY962518", "5'ETS_KY962518", "18S_KY962518", "ITS1_KY962518", "5.8S_KY962518", "ITS2_KY962518", 
         "28S_KY962518", "3'ETS_KY962518", "IGS_KY962518", "entire_rdna_KY962518", "coding_rdna", "no_igs_KY962518")


Sequences <- c(promoter,ets5, s18, its1, s5.8, its2, s28, ets3,igs, entire_rdna,coding_rdna, entire_rdna_no_igs)

Details <- c("1299_3500","3501_7157", "7158_9026", "9027_10096", "10097_10253",
             "10254_11420", "11421_16471", "16472_16832","16833_48338", "1299_46136", "3501_16832", "1299_16832")



rdna_hg38_dataset<- data.frame(Name, Sequences, Details)

attempt1<- data.frame((matrix(nrow = 0, ncol=5)))
for (j in 1: nrow(rdna_hg38_dataset)){
  a<- nchar (rdna_hg38_dataset[j,2])
  b<- str_count(rdna_hg38_dataset[j,2], "A")
  c<- str_count(rdna_hg38_dataset[j,2], "T")
  d<- str_count(rdna_hg38_dataset[j,2], "G")
  e<- str_count(rdna_hg38_dataset[j,2], "C")
  f <- c(a, b, c, d, e)
  attempt1[nrow(attempt1)+1,] <- f
}

colnames(attempt1) <- c("Total_nucleotides", "A", "T", "G", "C")
rdna_hg38_dataset_v2<- cbind(rdna_hg38_dataset, attempt1)


rdna_hg38_dataset_details_v2 <- rdna_hg38_dataset_v2 %>% 
  mutate(across(all_of(c(5,6,7,8)), function(x) x/Total_nucleotides*100, .names = "{col}%")) %>% 
  mutate(across(all_of(c(9,10,11,12)), function (x) round(x, 2))) %>% 
  select(1,2,3,4,5,9,6,10,7,11,8,12)


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/output")

fwrite(rdna_hg38_dataset_details_v2, 
       file = "rdna_hg38_chr21_2018_dataset_details_v2.csv")

##filter rows that contain a certain string

rdna_hg38_dataset_sequences<- rdna_hg38_dataset_details_v2 %>% select(1,2)


##create FASTA format
for (j in 1: nrow(rdna_hg38_dataset_sequences)[1:9]){
  
  write(paste(">",rdna_hg38_dataset_sequences[j,1], sep=''),                                            
        file = "rDNA_KY962518_2018_sequence_fasta_format.txt",
        append = TRUE)
  
  write(rdna_hg38_dataset_sequences[j,2],                                            
        file = "rDNA_KY962518_2018_sequence_fasta_format.txt",
        append = TRUE)
  
  write("\n",                                            
        file = "rDNA_KY962518_2018_sequence_fasta_format.txt",
        append = TRUE)
  
}



#refer human_rdna_gc_skew.R code to see calculation of GC percent, GC skew and AT percent

rdna_human<- fread("rdna_hg38_chr21_2018_dataset_details_v3.csv", sep = ",", header = TRUE)

#here i will be adding GC and AT percent 

rdna_2018<-rdna_human[c(1:9), c("Name", "A%", "G%", "C%", "T%", "GC_perc", "AT_perc")]
rdna_2018_no_igs<- rdna_human[c(1:8), c("Name", "A%", "G%", "C%", "T%", "GC_perc", "AT_perc")]

datasets<- list("rdna_2018"= rdna_2018,
                "rdna_2018_no_igs" = rdna_2018_no_igs)


for (nm in names(datasets)) {
  
  i <- datasets[[nm]]   # get the dataset by name
  
  nucleotide_new <- separate(i, Name, "Name", sep = "_")
  
  nucleotide_reshape <- pivot_longer(nucleotide_new, -Name, 
                                     names_to = "Nucleotide", 
                                     values_to = "Percent")
  
  # save with the dataset name in the filename
  fwrite(nucleotide_reshape,
         file = paste0("nucleotide_", nm, "_graph_input2.csv"))
  
  nucleotide_reshape$Name <- factor(
    nucleotide_reshape$Name,
    levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
               "ITS2","28S", "3'ETS", "IGS")
  )
  
  nucleotide_reshape$Nucleotide <- factor(
    nucleotide_reshape$Nucleotide,
    levels = c("A%", "T%", "G%", "C%", "AT_perc", "GC_perc")  # nucleotides first
  )
  
  nucleotide_labels <- c("A%", "T%", "G%", "C%", "AT %", "GC %")
  
  p <- ggplot(nucleotide_reshape, aes(x = Name, y = Percent, fill = Nucleotide)) + 
    geom_bar(stat = 'identity', color = "black") +
    theme(axis.text.x = element_text(angle = 45, size = 30, hjust = 1)) +
    labs(title = "Nucleotide distribution percent in human rDNA",
         subtitle = "KY962518",
         x = "rDNA region",
         y = "Nucleotide distribution percent") +
    geom_text(aes(label = Percent), position = position_stack(vjust = 0.5), size = 12) +
    scale_fill_manual(
      values = c(
        "A%" = "#FF9999",      # soft red
        "T%" = "#E7C318",      # orange
        "G%" = "#99CCFF",      # blue
        "C%" = "#66CC66",      # green
        "AT_perc" = "#FFCC99", # gold for AT
        "GC_perc" = "#b194c2"  # purple for GC
      ),
      breaks = c("A%", "T%", "G%", "C%", "AT_perc", "GC_perc"),
      labels = c("A", "T", "G", "C", "AT%", "GC%")
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 30),
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 30),
          panel.background = element_blank(),   # removes grey background
          plot.background  = element_blank())
  
  ggsave(paste0(nm, "_entire_sequences_nucleotide_distribution.tiff"),
         plot = p, width = 15, height = 10, dpi = 150)

  atgc <- nucleotide_reshape |>
    dplyr::filter(Nucleotide %in% c("AT_perc", "GC_perc")) |>
    dplyr::mutate(Nucleotide = factor(Nucleotide, levels = c("AT_perc", "GC_perc")))
  
  
  
 p_atgc<- ggplot(atgc, aes(x = Name, y = Percent, fill = Nucleotide)) + 
    geom_bar(stat = 'identity', color = "black") +
    theme(axis.text.x = element_text(angle = 45, size = 30, hjust = 1)) +
    labs(title = "AT vs GC percent in human rDNA",
         subtitle = "KY962518",
         x = "rDNA region",
         y = "Percent") +
   scale_fill_manual(
     values = c(
       "AT_perc" = "#FFCC99", # gold for AT
       "GC_perc" = "#b194c2"  # purple for GC
     ),
     breaks = c("AT_perc", "GC_perc"),
     labels = c("AT%", "GC%")
   ) +
    geom_text(aes(label = Percent), position = position_stack(vjust = 0.5), size = 12) + #this will increase the font of the text inside the bar
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 30),
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 30),
          panel.background = element_blank(),   # removes grey background
          plot.background  = element_blank())
  
  ggsave(paste0(nm, "_AT_vs_GC_sequences_nucleotide_distribution.tiff"),
         plot = p_atgc, width = 15, height = 10, dpi = 150)
  
  
  
  atgc_only <- nucleotide_reshape |>
    dplyr::filter(Nucleotide %in% c("A%", "T%", "G%", "C%")) |>
    dplyr::mutate(Nucleotide = factor(Nucleotide, levels = c("A%", "T%", "G%", "C%")))
  
  
  p_atgc_only<- ggplot(atgc_only, aes(x = Name, y = Percent, fill = Nucleotide)) + 
    geom_bar(stat = 'identity', color = "black") +
    theme(axis.text.x = element_text(angle = 45, size = 30, hjust = 1)) +
    labs(title = "Nucleotide percent in human rDNA",
         subtitle = "KY962518",
         x = "rDNA region",
         y = "Percent") +
    scale_fill_manual(
      values = c(
        "A%" = "#FF9999",      # soft red
        "T%" = "#E7C318",      # orange
        "G%" = "#99CCFF",      # blue
        "C%" = "#66CC66"    # green

      ),
      breaks = c("A%", "T%", "G%", "C%"),
      labels = c("A%", "T%", "G%", "C%")
    ) +
    geom_text(aes(label = Percent), position = position_stack(vjust = 0.5), size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 30),
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 30),
          panel.background = element_blank(),   # removes grey background
          plot.background  = element_blank())
  
  ggsave(paste0(nm, "_ATGC_only_sequences_nucleotide_distribution.tiff"),
         plot = p_atgc_only, width = 15, height = 10, dpi = 150)
  }



