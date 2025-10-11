# ------------------------------------------------------------------------------
# This code is part of paper: In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus.
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
#   This script computes and visualizes GC content (%) and GC skew [(G−C)/(G+C)]
#   across the human rDNA transcription unit (promoter to 3′ETS) using custom
#   functions (`gc_content()` and `gc_skew()`). Sliding-window analysis highlights
#   regional variation in base composition and strand bias along the rDNA locus.
#
# Major Steps:
#   1. Source the GC content and GC skew functions.
#   2. Load required libraries (stringr, tidyverse, data.table, Biostrings).
#   3. Read the extended human rDNA FASTA sequence (KY962518 + upstream IGS).
#   4. Extract the 5′ETS–3′ETS region (positions 1299–16832).
#   5. Compute GC content and GC skew using a user-defined sliding window (e.g., 100 bp).
#   6. Export sliding-window data tables as CSV files.
#   7. Generate publication-ready plots:
#        - Fig. 1E: GC content (%) profile
#        - Fig. 1F: GC skew profile
#
# Input:
#   - KY962518_added_3500nt_IGS_upstream_nontemplate.fasta
#       Extended human rDNA reference sequence (non-template strand).
#
# Output:
#   - KY962518_5ETS_TO_3ETS_gc_content_sliding_<window>bp.csv
#       GC content values for each window across the rDNA.
#   - KY962518_5ETS_TO_3ETS_gc_skew_sliding_<window>bp.csv
#       GC skew values for each window across the rDNA.
#   - KY962518_5ETS_TO_3ETS_gc_content_sliding_<window>bp.png  (Fig. 1E)
#   - KY962518_5ETS_TO_3ETS_gc_skew_sliding_<window>bp.png     (Fig. 1F)
#
# Notes:
#   - Sliding windows are overlapping; step size = 1 bp.
#   - Reference positions correspond to KY962518 numbering (start = 1299).
#   - Dashed horizontal lines indicate genome-average GC (~40.9%) or neutral skew (0).
#   - Figures are exported at 600 dpi for publication quality.
#
# Dependencies:
#   gc_content_function.R
#   gc_skew_function.R
#
# ------------------------------------------------------------------------------



# I made a loop that will calculate GC skew and GC content in 40, 70, 100 overlapping and non-overlapping window size for rdna locus (as a whole not compartmentalized) 
#starting from 5'ETS and ending at 3'ETS).

#load GC skew and GC content function
#calculate GC skew
source("./gc_skew_function.R")

#calculate GC content
source("./GC_content_function.R")


#load libraries
library(stringr) #needed for GC skew function
library(tidyverse)
library(data.table)
library(Biostrings)






# i did GC content of entire rdna 

gc_content(entire_rdna) #entire rdna promoter to igs (1299 to 46136)
#G_count C_count gc_content_value gc_content_perc
#1   11682   14354        0.5806682        58.06682

gc_content(entire_rdna_no_igs) #entire rdna promoter to igs (1299 to 16832)
#G_count C_count gc_content_value gc_content_perc
#1    5535    5446         0.706901         70.6901




rdna_human_seq<- readDNAStringSet("KY962518_added_3500nt_IGS_upstream_nontemplate.fasta")
seq<- as.character(rdna_human_seq[[1]]) 
seq2<- str_sub(seq, 1299, 16832) #just decided to have 5'ETS to 3'ETS
nchar(seq2) 
#[1] 15534

#bin_size=c(30, 70, 100)
bin_size = c(100)
for ( i in bin_size){
  
  gc_skew_data<- gc_skew(seq2, window_size = i)
  sliding_skew<- gc_skew_data$sliding_window_results
  fwrite(sliding_skew, paste0("KY962518_5ETS_TO_3ETS_gc_skew_sliding_data_", i,"bp.csv"))
  
  
  #if you want to run GC skew separately: 
  #for ( i in bin_size){
  #filename1<- paste0("KY962518_5ETS_TO_3ETS_gc_skew_sliding_data_", i,"bp.csv")
  #sliding_skew<- fread(filename1, sep = ",", header = TRUE)
  
  #filename2<- paste0("KY962518_5ETS_TO_3ETS_gc_skew_non_sliding_data_", i,"bp.csv")
  #non_sliding_skew<- fread(filename2, sep = ",", header = TRUE)
  
  gc_content_data<- gc_content(seq2, window_size = i)
  sliding_content<- gc_content_data$sliding_window_results
  fwrite(sliding_content, paste0("KY962518_5ETS_TO_3ETS_gc_content_sliding_data_", i,"bp.csv"))
  
  
  #Fig1E
  sliding_content_graph<- ggplot(sliding_content, aes(x = start, y = gc_content_perc)) +
    geom_line(color = "#E21515", size = 1.5) +
    #geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth =2.0) +
    geom_hline(yintercept = 40.89, linetype = "dashed", color="darkgreen", linewidth =3.0)+
    labs(#title = "GC Content across human rDNA",
         #subtitle= paste0("Sliding window of 100 nt"),
         x = "Position (bp)",
         y = "GC Content (%)") +
    scale_y_continuous(limits= c(0,100), breaks=seq(0,100, by=25))+
    
    scale_x_continuous(breaks = c(0, 2202, 5858, 7726, 8795, 8951, 15167, 15527), 
                       labels =c("Pro","5'ETS", "18S","ITS1", "5.8S", "ITS2", "28S", "3'ETS"))+
    
    theme(plot.title = element_text(hjust = 0.5), #, face = "bold"), #at this font size its giving illusion that its bold when its not.
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 80), #face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 6),  # <-- the rectangle
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, margin = margin(r=20)), # Center Y-axis title
          axis.ticks.y = element_line(color = "black", linewidth = 4),
          axis.ticks.length = unit(30, "pt"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20, color = "black"),
          axis.text.y = element_text(color = "black"))
  
  ggsave(paste0("KY962518_5ETS_TO_3ETS_gc_content_sliding_", i, "bp.png"), 
         plot =  sliding_content_graph, width=18, height = 10, dpi = 600)
  
  
  #Fig1F
  sliding_skew_graph<- ggplot(sliding_skew, aes(x = start, y = GC_skew_value)) +
    geom_line(color = "#E21515", size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth =3.0) +
    labs(#title = "GC Skew across human rDNA",
         #subtitle= paste0("Sliding window of 100 nt"),
         x = "Position (bp)", #ggplot bydefault need x and y axis label
        y = "GC Skew") +
    scale_x_continuous(breaks = c(0, 2202, 5858, 7726, 8795, 8951, 15167, 15527), 
                       labels =c("Pro","5'ETS", "18S","ITS1", "5.8S", "ITS2", "28S", "3'ETS"))+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 80), #face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 6),  # <-- the rectangle
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.7, margin = margin(r=20)), # Center Y-axis title
          axis.ticks.y = element_line(color = "black", linewidth = 4),
          axis.ticks.length = unit(30, "pt"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20, color = "black"),
          axis.text.y = element_text(color = "black"))
  
  ggsave(paste0("KY962518_5ETS_TO_3ETS_gc_skew_sliding_", i, "bp.png"), 
         plot = sliding_skew_graph, width=18, height = 10, dpi = 600)
  

 
  
} 




 





#######All are extra stuff that is no longer needed
{
  #calculate CpG content
  #source("./CpG_dinucleotides_function.R")
  rdna_human<- fread("rdna_hg38_chr21_2018_dataset_details_v3.csv", sep = ",", header = TRUE)
  
  full_length_gc_skew<- ggplot(rdna_human, aes(x = x_axis, y = norm_GC_skew)) +
    geom_line(color = "#E21515") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "GC Skew Across rDNA Sequence KY962518",
         x = "rDNA",
         y = "Normalized GC Skew") +
    scale_x_continuous(breaks = c(2202, 5859, 7728,8798, 8955, 10122, 15173, 15534, 47040 ),
                       labels = c("promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS", "IGS"))+
    
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20, color = "black"),
          axis.text.y = element_text(color = "black"))
  
  ggsave("KY962518_full_length_normalised_gc_skew.tiff", 
         plot = full_length_gc_skew, width=18, height = 10, dpi = 150)
  
  
  full_length_gc_skew_excld_igs<- ggplot(rdna_human[1:8,], aes(x = x_axis, y = norm_GC_skew)) +
    geom_line(color = "#E21515") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "GC Skew Across rDNA Sequence KY962518 (Exculding IGS)",
         x = "rDNA",
         y = "Normalised GC Skew") +
    scale_x_continuous(breaks = c(2202, 5859, 7728,8798, 8955, 10122, 15173, 15534 ),
                       labels = c("promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS"))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),#if you want to add a rectangle box or you can use theme_minimal()
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20, color = "black"),
          axis.text.y = element_text(color = "black"))
  
  ggsave("KY962518_full_length_excld_IGS_normalised_gc_skew.tiff", 
         plot = full_length_gc_skew_excld_igs, width=18, height = 10, dpi = 150)
  
  #plotting GC content
  full_length_gc_content_excld_igs<- ggplot(rdna_human[1:8,], aes(x = x_axis, y = GC_perc)) +
    geom_line(color = "#E21515")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 40.89, linetype = "dashed", color="darkgreen")+
    labs(title = "GC Content Across rDNA Sequence KY962518 (Exculding IGS)",
         x = "rDNA",
         y = "GC Content") +
    scale_x_continuous(breaks = c(2202, 5859, 7728,8798, 8955, 10122, 15173, 15534 ),
                       labels = c("promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS"))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),#if you want to add a rectangle box or you can use theme_minimal()
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20, color = "black"),
          axis.text.y = element_text(color = "black"))
  
  ggsave("KY962518_full_length_excld_IGS_GC_content.tiff", 
         plot = full_length_gc_content_excld_igs, width=18, height = 10, dpi = 150)
  
  
  
  #plotting CpG dinucleotide content
  full_length_CpG_excld_igs<- ggplot(rdna_human[1:8,], aes(x = x_axis, y = CpG_OE)) +
    geom_line(color = "#E21515")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0.6, linetype = "dashed", color="darkorange")+
    labs(title = "CpG dinucleotides O/E ratio Across rDNA Sequence KY962518 (Exculding IGS)",
         x = "rDNA",
         y = "CpG dinucleotides O/E ratio") +
    scale_x_continuous(breaks = c(2202, 5859, 7728,8798, 8955, 10122, 15173, 15534 ),
                       labels = c("promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS"))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),#if you want to add a rectangle box or you can use theme_minimal()
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20, color = "black"),
          axis.text.y = element_text(color = "black"))
  
  ggsave("KY962518_full_length_excld_IGS_CpG_dinucleotides_OE_ratio.tiff", 
         plot = full_length_CpG_excld_igs, width=18, height = 10, dpi = 150)
  
  
  full_length_CpG_perc_excld_igs<- ggplot(rdna_human[1:8,], aes(x = x_axis, y = CpG_perc)) +
    geom_line(color = "#E21515")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "CpG dinucleotides rDNA Sequence KY962518 (Exculding IGS)",
         x = "rDNA",
         y = "CpG dinucleotides percent") +
    scale_x_continuous(breaks = c(2202, 5859, 7728,8798, 8955, 10122, 15173, 15534 ),
                       labels = c("promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS"))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),#if you want to add a rectangle box or you can use theme_minimal()
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20, color = "black"),
          axis.text.y = element_text(color = "black"))
  
  ggsave("KY962518_full_length_excld_IGS_CpG_dinucleotides_perc.tiff", 
         plot = full_length_CpG_perc_excld_igs, width=18, height = 10, dpi = 150)
  
  
# there is no need to make template ones because it will be mirror image of non -template

# Above result showed that there is a GC skew rDNA region. the maximum is seen in promoter region, followed by 28s region than 5'ETS. 
# lowest is seen in 18S.

 #CpG dinucleotides 
# as per Gardiner-Garden and Frommer: paper title : CpG Islands in vertebrate genomes (1987) teh CpG has to have a length of 200, GC >50%, and CpG O/E >0.6
# but in 2002 a new paper published takai and Jones (Comprehensive analysis of CpG islands in human chromosomes 21 and 22)
# they said that earlier formula was made when sequencing was not performed, contain repetitive and alu elements. 
# This demanded new formula. as per which length of 500bp, GC% >55%, and CpG O/E >0.65 is needed

bin_size = c(200, 500)
for ( i in bin_size){

  CpG_dinucleotides_data<- CpG_content(seq2, window_size = i)
  sliding_CpG<- CpG_dinucleotides_data$sliding_window_results
  fwrite(sliding_CpG, paste0("KY962518_5ETS_TO_3ETS_CpG_dinucleotides_sliding_data_", i,"bp.csv"))
  
  
  sliding_CpG_graph<- ggplot(sliding_CpG, aes(x = start, y = CpG_OE)) +
    geom_line(color = "#E21515") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0.65, linetype = "dashed", color="orange")+
    labs(title = "CpG dinucleotides O/E ratio across rDNA Sequence KY962518",
         subtitle= paste0(i, " bp sliding window size"),
         x = "Position (bp)",
         y = "CpG dinucleotides O/E ratio") +
    scale_x_continuous(breaks = c(0, 2202, 5859, 7728,8798, 8955, 10122, 15173, 15534), 
                       labels =c("Promoter","5'ETS", "18S","ITS1", "5.8S", "ITS2", "28S", "3'ETS", "IGS"))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20, color = "black"),
          axis.text.y = element_text(color = "black"))
  
  ggsave(paste0("KY962518_5ETS_TO_3ETS_sliding_CpG_graph_", i, "bp.tiff"), 
         plot =  sliding_CpG_graph, width=18, height = 10, dpi = 150)
  
  
  
}  

# if you want to run non-sliding skew
non_sliding_skew<- gc_skew_data$fixed_window_results
fwrite(non_sliding_skew, paste0("KY962518_5ETS_TO_3ETS_gc_skew_non_sliding_data_", i,"bp.csv"))

  non_sliding_graph<- ggplot(non_sliding_skew, aes(x = start, y = GC_skew_value)) +
    geom_line(color = "#E21515") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "GC Skew Across rDNA Sequence KY962518",
         subtitle= paste0(i, " bp non sliding window size"),
         x = "Position (bp)",
        y = "GC Skew (G−C / G+C)") +
    scale_x_continuous(breaks = c(0, 3657,5526, 6596, 6753, 7920, 12971, 13333), 
                       labels =c("5'ETS", "18S","ITS1", "5.8S", "ITS2", "28S", "3'ETS", "IGS"))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20, color = "black"),
          axis.text.y = element_text(color = "black"))
  
  ggsave(paste0("KY962518_5ETS_TO_3ETS_gc_skew_non_sliding_", i, "bp.tiff"), 
         plot = non_sliding_graph, width=18, height = 10, dpi = 150)
  
  non_sliding<- paste0("KY962518_5ETS_TO_3ETS_gc_content_non_sliding_data_", i,"bp.csv")
  non_sliding_content<- fread(non_sliding, sep = ",", header = TRUE)
  
  
  #if you want to re-run only GC content thing
  
  #for ( i in bin_size){
    
    
  
  #if you want to run non sliding GC content
  non_sliding_content<- gc_content_data$fixed_window_results
  fwrite(non_sliding_content, paste0("KY962518_5ETS_TO_3ETS_gc_content_non_sliding_data_", i,"bp.csv"))
  
  non_sliding_content_graph<- ggplot(non_sliding_content, aes(x = start, y = gc_content_perc)) +
    geom_line(color = "#E21515") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 40.89, linetype = "dashed", color="darkgreen")+
    labs(title = "GC Content Across rDNA Sequence KY962518",
         subtitle= paste0(i, " bp non sliding window size"),
         x = "Position (bp)",
        y = "GC Fraction (G+C / A+T+G+C)") +
    scale_y_continuous(limits= c(0,100), breaks=seq(0,100, by=20))+
    scale_x_continuous(breaks = c(0, 3657,5526, 6596, 6753, 7920, 12971, 13333), 
                       labels =c("5'ETS", "18S","ITS1", "5.8S", "ITS2", "28S", "3'ETS", "IGS"))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20, color = "black"),
          axis.text.y = element_text(color = "black"))
  
  ggsave(paste0("KY962518_5ETS_TO_3ETS_gc_content_non_sliding_", i, "bp.tiff"), 
         plot = non_sliding_content_graph, width=18, height = 10, dpi = 150)
  
}


  
  



