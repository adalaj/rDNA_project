# ------------------------------------------------------------------------------
# This code is part of paper: In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus.
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
#   This script extracts individual regions from the human rDNA reference
#   sequence (KY962518, with added upstream 3.5 kb IGS extension) and computes
#   nucleotide composition (A, T, G, C counts and percentages) for each region.
#   The script also generates publication-quality bar plots comparing AT vs GC
#   content and individual base distribution across annotated rDNA segments.
#
# Major Steps:
#   1. Load extended human rDNA sequence (KY962518_added_3500nt_IGS_upstream).
#   2. Define functional subregions (promoter, ETS, ITS, 18S, 5.8S, 28S, IGS).
#   3. Compute total nucleotide counts (A, T, G, C) and percentages.
#   4. Save composition tables as CSV and FASTA files for downstream analyses.
#   5. Generate bar plots:
#        - AT vs GC content across rDNA subregions (Fig. 3)
#        - Individual nucleotide distribution (A, T, G, C; Fig. 4)
#        - Entire rDNA composition summary (Fig. 2)
#
# Input:
#   - KY962518_added_3500nt_IGS_upstream_nontemplate.fasta
#       Extended human rDNA reference sequence including upstream IGS.
#
# Output:
#   - rdna_hg38_chr21_2018_dataset_details_v2.csv
#       Nucleotide counts and percentages for all defined regions.
#   - rDNA_KY962518_2018_sequence_fasta_format.txt
#       FASTA-formatted sequences for promoter, ETS, ITS, rRNAs, and IGS.
#   - rdna_2018_no_igs_AT_vs_GC_sequences_nucleotide_distribution.tiff
#       Bar plot of AT vs GC content.
#   - rdna_2018_no_igs_ATGC_only_sequences_nucleotide_distribution.tiff
#       Bar plot of individual nucleotide percentages.
#   - entire_rdna_sequences_nucleotide_distribution.tiff
#       Summary of overall rDNA nucleotide composition.
#
# Notes:
#   - Sequence boundaries follow KY962518 annotations with added upstream IGS.
#   - Base composition is calculated using `stringr::str_count`.
#   - Percentages are normalized to total nucleotides per region.
#   - Visualization performed using ggplot2 (part of tidyverse).
#   - Output figures are in TIFF format for publication use.
#
#------------------------------------------------------------------------------




#load libraries
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



rdna_human<- fread ("rdna_hg38_chr21_2018_dataset_details_v2.csv", sep = ",", header = TRUE)
attempt1<- data.frame((matrix(nrow = 0, ncol=6))) #6 coloumns because Gc skew has 6 columns

#calculate GC skew
source("./gc_skew_function.R")

for ( i in 1:nrow(rdna_human)){
  gc_skew_data<- gc_skew(rdna_human$Sequences[i])
  attempt1<- rbind(attempt1, gc_skew_data)
}

selected_gc_skew_data <- attempt1 %>% select(GC_skew_value) 
rdna_human<- cbind(rdna_human, selected_gc_skew_data)


#calculated GC skew for template strand as well
attempt2<- data.frame((matrix(nrow = 0, ncol=6)))
for ( i in 1:nrow(rdna_human)){
  rev_seq<- as.character(reverseComplement(DNAString(rdna_human$Sequences[i])))
  temp_gc_skew<- gc_skew(rev_seq)
  attempt2<- rbind(attempt2, temp_gc_skew)
  
}

template_GC_skew<- attempt2 %>% select(window_seq, G_count, C_count, GC_skew_value)
colnames(template_GC_skew)<- c("template_seq", "temp_G_count", "temp_C_count", "temp_GC_skew_value")
rdna_human<- cbind(rdna_human, template_GC_skew)


# wanted to add cumulative sum of total nucleotides for future 
rdna_human<- rdna_human %>% mutate(norm_GC_skew = GC_skew_value/sum(GC_skew_value))
rdna_human<- rdna_human %>% mutate(temp_norm_GC_skew = temp_GC_skew_value/sum(GC_skew_value))

rdna_human$x_axis <- cumsum(rdna_human$Total_nucleotides)

#calculated GC and AT perc
rdna_human<- rdna_human %>% mutate(GC_perc = round(((G+C)/(Total_nucleotides))*100,2))
rdna_human<- rdna_human %>% mutate(AT_perc = round(((A+T)/(Total_nucleotides))*100,2))


#also wanted to calculate CpG or CG dinucleotide

rdna_human <- rdna_human %>%
  mutate(
    CpG_count = str_count(Sequences, "CG"),
    CpG_perc = (CpG_count/Total_nucleotides)*100,
    CpG_OE = (CpG_count * Total_nucleotides) / (G * C))



fwrite(rdna_human, "rdna_hg38_chr21_2018_dataset_details_v3.csv")


rdna_human<- fread("rdna_hg38_chr21_2018_dataset_details_v3.csv", sep = ",", header = TRUE)

#here i will be adding GC and AT percent 

rdna_2018<-rdna_human[c(1:9), c("Name", "A%", "G%", "C%", "T%", "GC_perc", "AT_perc")]
rdna_2018_no_igs<- rdna_human[c(1:8), c("Name", "A%", "G%", "C%", "T%", "GC_perc", "AT_perc")]

datasets<- list(#"rdna_2018"= rdna_2018,
                "rdna_2018_no_igs" = rdna_2018_no_igs)


for (nm in names(datasets)) {
  
  i <- datasets[[nm]]   # get the dataset by name 
  #nucleotide_new <- separate(datasets[["rdna_2018_no_igs"]], Name, "Name", sep = "_")
  
  nucleotide_new <- separate(i, Name, "Name", sep = "_")
  
  nucleotide_reshape <- pivot_longer(nucleotide_new, -Name, 
                                     names_to = "Nucleotide", 
                                     values_to = "Percent")
  
  # save with the dataset name in the filename
  #fwrite(nucleotide_reshape,
         #file = paste0("nucleotide_", nm, "_graph_input2.csv"))
  
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
  
 
  
  atgc <- nucleotide_reshape |>
    dplyr::filter(Nucleotide %in% c("AT_perc", "GC_perc")) |>
    dplyr::mutate(Nucleotide = factor(Nucleotide, levels = c("AT_perc", "GC_perc")))
  
  
 
 p_atgc<- ggplot(atgc, aes(x = Name, y = Percent, fill = Nucleotide)) + 
    geom_bar(stat = 'identity', color = "black") +
    labs(#title = "AT vs GC percent comparison",
         x = "Human rDNA region",
         y = "AT and GC Content (%)", 
         fill = NULL) +
   scale_fill_manual(
     values = c(
       "AT_perc" = "#FFCC99", # gold for AT
       "GC_perc" = "#b194c2"  # purple for GC
     ),
     breaks = c("AT_perc", "GC_perc"),
     labels = c("AT%", "GC%")
   ) +
    geom_text(aes(label = Percent), position = position_stack(vjust = 0.5), size = 18)+ #fontface ="bold") + #this will increase the font of the text inside the bar
    theme(plot.title = element_text(hjust = 0.5), #face = "bold", size = 20),
          #plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 80),#face = "bold"),
          axis.line = element_line(color = "black", linewidth = 4),
          panel.grid = element_blank(),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, margin = margin (r=20)),
          axis.ticks = element_line(color = "black", linewidth = 4),
          axis.ticks.length = unit(30, "pt"),
          axis.text.x = element_text(angle = 45, hjust = 1.0, color = "black"),
          axis.text.y  = element_text(color = "black"),
          panel.background = element_blank(),   # removes grey background
          plot.background  = element_blank(), 
          legend.position = "top") 
 
  ggsave("rdna_2018_no_igs_AT_vs_GC_sequences_nucleotide_distribution.png",
         plot = p_atgc, width = 20, height = 15.5, dpi = 600)
  
  
  
  atgc_only <- nucleotide_reshape |>
    dplyr::filter(Nucleotide %in% c("A%", "T%", "G%", "C%")) |>
    dplyr::mutate(Nucleotide = factor(Nucleotide, levels = c("A%", "T%", "G%", "C%")))
  
  

  p_atgc_only<- ggplot(atgc_only, aes(x = Name, y = Percent, fill = Nucleotide)) + 
    geom_bar(stat = 'identity', color = "black") +
    #theme(axis.text.x = element_text(angle = 45, size = 30, hjust = 1)) +
    labs(#title = "Nucleotide distribution in human rDNA",
         x = "Human rDNA region",
         y = "Nucleotide Percent (%)",
         fill = NULL) +
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
    geom_text(aes(label = Percent), position = position_stack(vjust = 0.5), size = 18)+# fontface = "bold") +
    theme(plot.title = element_text(hjust = 0.5), #face = "bold", size=20),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 80), #face = "bold"),
          axis.line = element_line(color = "black", linewidth = 4),
          panel.grid = element_blank(),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, margin = margin(r=20)),
          axis.ticks = element_line(color = "black", linewidth = 4),
          axis.ticks.length = unit(30, "pt"),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y  = element_text(color = "black"),
          panel.background = element_blank(),   # removes grey background
          plot.background  = element_blank(),
          legend.position = "top") 
  
  ggsave( "rdna_2018_no_igs_ATGC_only_sequences_nucleotide_distribution.png",
         plot = p_atgc_only, width = 20, height = 16, dpi = 600)
  
  
  
  
  # dataset 1: individual bases
  df1 <- data.frame(
    Name = "Bases",
    Nucleotide = c("A%", "T%", "G%", "C%"),
    Percent = c(12.01, 15.55, 36.88, 35.56)
  )
  
  df1$Nucleotide <- factor(
    df1$Nucleotide,
    levels = c("A%", "T%", "G%", "C%")  # nucleotides first
  )
  
  
  # dataset 2: combined groups
  df2 <- data.frame(
    Name = "Groups",
    Nucleotide = c("AT%", "GC%"),
    Percent = c(27.56, 72.44)
  )
  
  df2$Nucleotide <- factor(
    df2$Nucleotide,
    levels = c("AT%", "GC%")  # nucleotides first
  )
  
  
  # combine both
  
  entire_only <- rbind(df1, df2)
  
 
  
  
  # plot
  entire_rdna<- ggplot(entire_only, aes(x = Name, y = Percent, fill = Nucleotide)) +
    geom_bar(stat = "identity", color = "black") +
    labs(
      #title = "Nucleotide Distribution",
      x = "Human rDNA",
      y = "Nucleotide Percent (%)",
      fill = NULL
    ) +
    scale_fill_manual(
      values = c(
        "A%" = "#FF9999",   # soft red
        "T%" = "#E7C318",   # orange
        "G%" = "#99CCFF",   # blue
        "C%" = "#66CC66",   # green
        "AT%" = "#FFCC99",   # purple 
        "GC%" = "#b194c2"    # dark red
      ),
      breaks = c("A%", "T%", "G%", "C%", "AT%", "GC%"),
      labels = c("A%", "T%", "G%", "C%", "AT%", "GC%")
      #guide = guide_legend(reverse = TRUE)
    ) +
    geom_text(aes(label = Percent), position = position_stack(vjust = 0.5), size = 18)+ #, fontface = "bold") +
    theme(
      plot.title    = element_text(hjust = 0.5), #face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0.5),
      text          = element_text(size = 80), #face = "bold"),
      axis.line     = element_line(color = "black", linewidth = 4),
      panel.grid    = element_blank(),
      axis.title.y  = element_text(angle = 90, vjust = 0.5, hjust = 0.5, margin = margin(r = 20)),
      axis.ticks    = element_line(color = "black", linewidth = 4),
      axis.ticks.length = unit(30, "pt"),
      axis.text.x   = element_text(angle = 45, hjust = 1, color = "black"), # keep straight under bars
      axis.text.y   = element_text(color = "black"),
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.position  = "right"
    )
  
  ggsave( "entire_rdna_sequences_nucleotide_distribution.png",
          plot = entire_rdna, width = 12, height = 16, dpi = 600)
  
  
  }


