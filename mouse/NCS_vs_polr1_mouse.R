# ------------------------------------------------------------------------------
# This code is part of paper: In silico Mapping of Non-Canonical DNA Structures Across the mouse Ribosomal DNA Locus.
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
#   This script quantifies and visualizes the relationship between transcriptional
#   regulators (POLR1A) and non-canonical DNA structures (R-loops, G4s, i-motifs)
#   across the modified mouse rDNA locus (*chrR*, 1–44,838 bp).
#   It generates binned correlation maps, zoomed promoter-to-3′ETS profiles,
#   and pairwise correlation plots for Figure 5.
#
# Major Steps:
#   1. Define 100-bin segmentation of *chrR* region (rDNA locus from hg38-rDNA_v1.0.bed).
#   2. Integrate POLR1A ChIP-seq counts with RLFS, G4FS, and iMFS annotations.
#   3. Normalize all datasets (0–1 scale) for cross-signal comparison.
#   4. Assess directional agreement (increase/decrease per bin) between POLR1A and each structure.
#   5. Compute Pearson and Spearman correlations and plot:
#        - POLR1A vs RLFS, G4FS, iMFS line plots (normalized trends)
#        - Corresponding scatterplots with correlation coefficients.
#   6. Generate genome schematic for rDNA context (full and zoomed black-and-white maps).
#
# Input:
#   - hg38-rDNA_v1.0.bed (from Vikram et al., 2020)
#   - ChIP_entire_rDNA_100bins_summary_raw_counts.tab (POLR1A ChIP-seq)
#   - RLFS_mouse_modified_[template/nontemplate].bed
#   - G4FS_mouse_modified_[template/nontemplate].bed
#   - imotif_mouse_modified_[template/nontemplate].bed
#
# Output:
#   - *_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv (graph-ready merged data)
#   - *_RLFS_with_POLR1A.png / *_G4FS_with_POLR1A.png / *_iMFS_with_POLR1A.png (Fig 6 and Supplementary table 6)
#   - *_pearson.png (scatterplots with r and P values) (Fig 6)
#   - rdna_black_and_white.png, rdna_zoom_black_and_white.png (rDNA region schematics)
#
# Notes:
#   - RLFS = R-loop forming sequence; G4FS = predicted G4-forming sequence; iMFS = i-motif forming sequence.
#   - Normalized data (0–1) used to compare relative enrichment profiles.
#   - Directional inverse trends (POLR1A↑ vs RLFS↓, etc.) are quantified as percentage mismatches per bin.
#   - Figures correspond to correlation analyses in Figure 5 of the manuscript.
#
# ------------------------------------------------------------------------------


#step1: identity chrR region of interest that where i will count the read count of UBF and polr1

library(data.table)
library(tidyverse)
library(ggtext)
library(karyoploteR)
library(Biostrings)

#set working directory

modified_rdna<- fread("mm39-rDNA_v1.0.bed", sep = "\t", header = FALSE)
#mm39-rDNA_v1.0.bed file contains the bed format of modified mouse rdna used in 
#Vikram paper title: Construction and validation of customized genomes for mouse and mouse ribosomal DNA mapping    deposited in github

filt_modified_rdna<- modified_rdna %>% filter(V1 == "chrR")
filt_modified_rdna<- filt_modified_rdna %>% select(1:6) #bed format
colnames(filt_modified_rdna)<- c("chr", "start", "end", "name", "score", "strand")
filt_modified_rdna<- filt_modified_rdna %>% mutate(length = (end-start))

#here it starts from 1 which part of IGS and end at 44838 which is also part of IGS
# plan is to make 100 bins in this range 1:44838 map RLFS, imotif, G4s and then a zoom image from promoter to 3'ETS+ 2000 


#step2 look whether 100 is enough
bin_size <- c(100) #70,30, 50, 75, 100, 200, 500, 1000
for (i in bin_size){
  bin_size_new<- i +1
  bin_break <- seq(1, 45306, length.out= bin_size_new)
  bin_width<- diff(bin_break)#[1]
  bin_data<- data.frame(
    bin_lower_value= round(head(bin_break,-1)), #rounded because it was decimel
    bin_upper_value= round(tail(bin_break,-1)),
    bin_midpoints = round((head(bin_break,-1) + tail(bin_break,-1))/2)
  )
}


bin_data$chr<- "chrR"

entire_rdna_100bins<- bin_data %>% select(4,1,2)
fwrite(entire_rdna_100bins, "chrR_mouse_entire_rdna_100bins.bed", sep = "\t")
#manually removed header 

#also saving bin_data because it has binpoints that will be needed for plotting
#fwrite(bin_data, "entire_rdna_100bins.csv")


#Last login: Sun Apr  6 23:50:07 on ttys000
#(base) jyotiadala@Mac ~ % cd Downloads 
#(base) jyotiadala@Mac Downloads % conda activate python3.11
#(python3.11) jyotiadala@Mac Downloads % python g4_canonical_finder_3.11python.py Mouse_BK000964.3_Modified_Snapgene.fasta > Output_G4FS_Mouse_BK000964.3_Modified_Snapgene.txt  

#(python3.11) jyotiadala@Mac Downloads % conda activate python2.7
#(python2.7) jyotiadala@Mac Downloads % python QmRLFS-finder.py -bed -i Mouse_BK000964.3_Modified_Snapgene.fasta -o Mouse_BK000964.3_Modified_Snapgene_qmRLFS      
#QmRLFS-finder.py (version v1.5)
#run on Mon Apr 07 2025 00:00:37 
#command line: python QmRLFS-finder.py -bed -i Mouse_BK000964.3_Modified_Snapgene.fasta -o Mouse_BK000964.3_Modified_Snapgene_qmRLFS
#Time used: 0.35 mins

RLFS<- fread("Mouse_BK000964.3_Modified_Snapgene_qmRLFS.out.bed", sep = "\t", header = FALSE) #68
G4FS<- fread("Output_G4FS_Mouse_BK000964.3_Modified_Snapgene.txt", sep = "\t", header = FALSE) #87



# i have to add strand colur 
RLFS <- RLFS %>% mutate(V9 = ifelse(V6 == "+", "65,105,225", "205,50,120"))
RLFS$V1<- "chrR"
fwrite(RLFS, "RLFS_Mouse_modified_BK000964_Modified.bed", sep = "\t")
#removed header manually


G4FS<- G4FS %>% mutate(V4= 0) # chnage sequece information to zero.
G4FS<- G4FS %>% select(V1, V2, V3, V5, V4, V6) # change the order as per bed file
G4FS<- G4FS %>% mutate(V7 = V2)
G4FS<- G4FS %>% mutate(V8 =V3)
G4FS <- G4FS %>% mutate(V9 = ifelse(V6 == "+", "65,105,225", "205,50,120"))
G4FS<- G4FS %>% mutate(V10 = 2)
G4FS<- G4FS %>% mutate(V11 = V2)
G4FS<-G4FS %>% mutate(V12 = V3)
G4FS$V1<- "chrR"
fwrite(G4FS, "G4FS_Mouse_modified_BK000964_Modified.bed", sep = "\t")
#removed header manually


#prepare fasta file for imotif detection, because the snapgene fasta is in lowercase and iM finder cannot predict
mouse_rDNA<- readDNAStringSet(file = "Mouse_BK000964.3_Modified_Snapgene.fasta") #belongs to package BIOSTRINGS
#this modified fasta is from paralkar paper
mouse_rDNA_seq<- mouse_rDNA[[1]]
nchar(mouse_rDNA_seq)
# 45306


template_mouse_rdna<- reverseComplement(mouse_rDNA_seq)
nchar(template_mouse_rdna)
#45306



# Make uppercase and keep as DNAString
nontemplate_upper <- DNAString(toupper(mouse_rDNA_seq))
template_upper <- DNAString(toupper(template_mouse_rdna))


nontemplate_upper_set<- DNAStringSet(nontemplate_upper)
template_upper_set<- DNAStringSet(template_upper)

# Name the sequence (important for FASTA headers)
names(nontemplate_upper_set) <- "Nontemplate_BK000964"
names(template_upper_set) <- "Template_BK000964"


# Write to FASTA
writeXStringSet(nontemplate_upper_set, filepath = "Mouse_BK000964.3_Modified_nontemplate.fasta")
writeXStringSet(template_upper_set, filepath = "Mouse_BK000964.3_Modified_template.fasta")


#for i motif 
# I added mouse_BK000964.1_Modified_nontemplate.fasta and mouse_BK000964.1_Modified_template.fasta to predict>upload DNA FILE (FASTA FORMAT) (https://im-seeker.org/predict)>end to end with default model>predict
# output saved as output_imotif_mouse_default_end_to_end_prediction_nontemplate_modified.csv (csv is the default saving option)
# and output_imotif_mouse_default_end_to_end_prediction_template_modified.csv

#i am deleting both teh fasta (template and nontemplate) files because its teh same information in uppercase


imotif_nontemplate<- fread("output_imotif_mouse_default_end_to_end_prediction_nontemplate_modified.csv", header = TRUE, sep = ",") #28
imotif_template<- fread("output_imotif_mouse_default_end_to_end_prediction_template_modified.csv", header = TRUE, sep = ",") #23
imotif_nontemplate$strand<- "+"
imotif_template$strand<- "-"

imotif_master<- rbind(imotif_nontemplate, imotif_template)
imotif_master$chr<- as.character(imotif_master$chr)
imotif_master$chr<- "chrR"


imotif_bed<- imotif_master %>% select(chr, beg, end, strand) 
imotif_bed$name<- "imotif"
imotif_bed$score<- 0
imotif_bed<- imotif_bed %>% select(chr, beg, end, name, score, strand) 

fwrite(imotif_bed, "iMFS_Mouse_BK000964_Modified.bed", sep = "\t")


###########################################################
#prepare template RLFS, G4s, imotif  vs POL I for making graphs

#count count of RLFS, iMFS, and G4FS in the rdna 100 bins
# Simple unnamed list
list_of_names <- list("iMFS", "G4FS", "RLFS")

# Define binning parameters
bin_size <- 100
bin_size_new <- bin_size + 1
bin_break <- seq(1, 45306, length.out = bin_size_new)
bin_width <- diff(bin_break)#[1]




# Prepare empty lists to store counts
all_counts_nontemplate <- list()
all_counts_template <- list()

# Loop through each structure
for (label in list_of_names) {
  
  # === Load and process non-template ===
  ncs_file_name <- paste0(label, "_Mouse_BK000964_Modified.bed")
  ncs_file <- fread(ncs_file_name, sep = "\t", header = FALSE)
  ncs_file$V1<- "chrR"
  nontemplate<- ncs_file %>% filter(ncs_file$V6=="+")
  nontemplate$actual_start <- nontemplate$V2
  nontemplate$actual_end <- nontemplate$V3
  nt_hist <- hist(nontemplate$actual_start, breaks = bin_break, plot = FALSE)
  all_counts_nontemplate[[label]] <- nt_hist$counts
  
  
  template<- ncs_file %>% filter(ncs_file$V6=="-")
  template$actual_start <- template$V3
  template$actual_end <- template$V2
  t_hist <- hist(template$actual_start, breaks = bin_break, plot = FALSE)
  all_counts_template[[label]] <- t_hist$counts
  
}

# Shared bin info
bin_info <- data.frame(
  bin_lower_value = round(head(bin_break, -1)),
  bin_upper_value = round(tail(bin_break, -1)),
  bin_midpoints   = round((head(bin_break, -1) + tail(bin_break, -1)) / 2)
)

# Construct final data.frames
combined_count_nontemplate <- cbind(bin_info, data.frame(
  G4FS_counts = all_counts_nontemplate[["G4FS"]],
  RLFS_counts  = all_counts_nontemplate[["RLFS"]],
  iMFS_counts = all_counts_nontemplate[["iMFS"]]
))

combined_count_template <- cbind(bin_info, data.frame(
  G4FS_counts = all_counts_template[["G4FS"]],
  RLFS_counts  = all_counts_template[["RLFS"]],
  iMFS_counts = all_counts_template[["iMFS"]]
))


##refer Deeptools_NCS_vs_mouse_rdna_pol_I_chip_seq.txt file to see commands i did for average bin summary using deeptools multibigwigsummary
chip<- fread("ChIP_mouse_rDNA_100bins_summary_raw_counts.tab", sep = "\t", header = TRUE)

# 1. Rename chip columns
colnames(chip) <- c("chr", "bin_start", "bin_end", "POLR1A")

# 2. Merge chip with structural counts (column-wise, assuming same number of rows)
chip_with_template <- cbind(chip, combined_count_template[, c("bin_midpoints","G4FS_counts", "RLFS_counts", "iMFS_counts")])
chip_with_nontemplate <- cbind(chip, combined_count_nontemplate[, c("bin_midpoints","G4FS_counts", "RLFS_counts", "iMFS_counts")])



data_of_interest<- list(
  Nontemplate = chip_with_nontemplate,
  Template = chip_with_template
)

for (i in names(data_of_interest)){
  
  #need to scale both the axis as between zero to 1 for all UBF, POLR1A, RLFS, G4FS, iMFS
  data_of_interest[[i]]<- data_of_interest[[i]] %>% mutate(norm_POLR1A= (POLR1A- min(POLR1A))/(max(POLR1A)-min(POLR1A)),
                                                           norm_G4FS_counts= (G4FS_counts- min(G4FS_counts))/(max(G4FS_counts)-min(G4FS_counts)),
                                                           norm_RLFS_counts= (RLFS_counts- min(RLFS_counts))/(max(RLFS_counts)-min(RLFS_counts)),
                                                           norm_iMFS_counts= (iMFS_counts- min(iMFS_counts))/(max(iMFS_counts)-min(iMFS_counts)))
  
  
  #When POLR1A increases from one bin to the next, does RLFS tend to decrease, and vice versa?
  #what percent of times, we see this scenario
  
  #idea is to the difference of current value subtracted by previous value. 
  
  #you are not measuring whether the value itself is high or low, but rather whether the signal changed upward or downward compared to the previous bin.
  
  #Then sign(diff) tells you:
  #+1 → increase compared to previous bin
  #-1 → decrease compared to previous bin
  # 0 → no change
  
  #example: Case 1: POLR1A value = 0.6, previous = 0.8
  #diff = -0.2, sign = -1 → interpreted as a decrease, even though the absolute signal is still positive.
  
  
  #This measures whether the signal increased or decreased compared to the previous bin.
  
  data_of_interest[[i]] <- data_of_interest[[i]] %>%
    mutate(
      #diff_norm_UBF = norm_UBF - lag(norm_UBF, default = 0),
      diff_norm_POLR1A = norm_POLR1A - lag(norm_POLR1A),
      diff_norm_G4FS_counts = norm_G4FS_counts - lag(norm_G4FS_counts),
      diff_norm_RLFS_counts = norm_RLFS_counts - lag(norm_RLFS_counts),
      diff_norm_iMFS_counts = norm_iMFS_counts - lag(norm_iMFS_counts)
    )
  
  
  #If one goes up and the other goes down (opposite signs), you mark it as 1 (a “mismatch” in direction)
  data_of_interest[[i]] <- data_of_interest[[i]] %>%
    mutate(POLR1_RLFS_sign_match = ifelse(sign(diff_norm_POLR1A) != sign(diff_norm_RLFS_counts), 1, 0))
  
  data_of_interest[[i]] <- data_of_interest[[i]] %>%
    mutate(POLR1_G4FS_sign_match = ifelse(sign(diff_norm_POLR1A) != sign(diff_norm_G4FS_counts), 1, 0))
  
  data_of_interest[[i]] <- data_of_interest[[i]] %>%
    mutate(POLR1_iMFS_sign_match = ifelse(sign(diff_norm_POLR1A) != sign(diff_norm_iMFS_counts), 1, 0))
  
  
  fwrite(data_of_interest[[i]], paste0(i, "_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv"))
  
}




###################################################
#to make graphs

#####read again

zoom_start <- 6000 #round off of promoter
zoom_end <- 23000 #round off of end of 3'ETS

#Nontemplate<- fread("Nontemplate_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv", sep = ",", header = TRUE)
Template<- fread("Template_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv", sep = ",", header = TRUE)

#Entire_zoom <- Entire %>% filter(bin_midpoints >= zoom_start & bin_midpoints <= zoom_end)
#Nontemplate_zoom <- Nontemplate %>% filter(bin_midpoints >= zoom_start & bin_midpoints <= zoom_end)
Template_zoom <- Template %>% filter(bin_midpoints >= zoom_start & bin_midpoints <= zoom_end)


region_marks <- c(1, 6307, 9307, 13314,      15184, 16184, 16341, 17429, 22159,22709)
region_labels <- c("","Pro", "5′ET", "18S", "ITS1", "5.8S", "ITS2", "28S", "3′ET","")




data_of_interest<- list(
  #Nontemplate = Nontemplate,
  #Nontemplate_zoom = Nontemplate_zoom,
  #Template = Template,
  Template_zoom = Template_zoom
)


for (i in names(data_of_interest)){
  cols <- c("RLFS" = "#aa2a85", "POLR1A" = "#D2691E")
  
  # legend-only data (NA coordinates => nothing is drawn on the panel)
  legend_df <- data.frame(
    bin_midpoints = NA_real_,
    y = NA_real_,
    grp = factor(c("POLR1A","RLFS"), levels = c("POLR1A","RLFS"))
  )
  
  
  ##plot the RLFS with POLR1A
  plot_rlfs_polr1a <-
    ggplot() +
    #geom_line(data = data, aes(x = bin_midpoints, y = UBF, color = "UBF"), size = 2.0) +
    geom_line(data = data_of_interest[[i]], aes(x = bin_midpoints, y = norm_POLR1A, color = "POLR1A"), size = 3.0,show.legend = FALSE) +
    
    geom_line(data = data_of_interest[[i]],
              aes(x = bin_midpoints, y = norm_RLFS_counts, color = "RLFS"), size = 3.0, show.legend = FALSE) +
    #geom_point(data = chip_with_nontemplate,
    #aes(x = bin_midpoints, y = norm_RLFS_counts, shape = "RLFS", color = "RLFS"), size = 3) +
    
    # legend-only filled squares
    geom_point(data = legend_df,
               aes(x = bin_midpoints, y = y, color = grp),
               shape = 15, size = 10, inherit.aes = FALSE, na.rm = TRUE) +
    
    scale_color_manual(name = NULL, 
                       values = cols)+
    #labels = c(
    #"RLFS"   = "<span style='color:#aa2a85;'>RLFS</span>",
    #"POLR1A" = "<span style='color:#D2691E;'>POLR1A</span>" #"UBF" = "#f609b4" #CD32A0, #9e277b
    #),
    
    guides(color = guide_legend(override.aes = list(shape = 15, size = 40, linetype = NA))) +
    
    #scale_shape_manual(name = "Non-canonical structure", values = c("RLFS" = 16)) +
    scale_x_continuous(breaks = region_marks, labels = region_labels) +
    scale_y_continuous(
      name = "POLR1A ChIP-seq signal",
      limits = c(0,1),
      breaks = seq(0, 1, 0.2),
      sec.axis = sec_axis(~., breaks = seq(0, 1, 0.2),
                          name = "RLFS counts")
    ) +
    labs(
      #title = paste0(i," RLFS and POLR1A"),
      x = "Mouse rDNA region"
    ) +
    theme(axis.ticks = element_line(color = "black", linewidth = 4),
          axis.ticks.length = unit(50, "pt"),
          #plot.title = element_text(hjust = 0.5, face = "bold", size=20),
          text = element_text(size = 150),
          panel.border = element_rect(color = "black",fill=NA, linewidth = 6),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),#if you want to add a rectangle box or you can use theme_minimal()
          axis.text.x = element_text(angle = 90, hjust = 1, size=30, color = "black", margin = margin(b=20)),
          axis.title.y.left  = element_text(color = "#D2691E", margin = margin (r=20)), #size = 100),  # left axis orange
          axis.title.y.right = element_text(color = "#aa2a85", margin = margin (l=20)), #size = 100),  # right axis purple
          axis.text.y.left   = element_text(color = "#D2691E"), #size = 150),
          axis.text.y.right  = element_text(color = "#aa2a85"), #size = 150),
          legend.position = c(0.05,0.98),
          legend.justification = c("left", "top"),
          legend.key = element_blank(),
          #legend.title = element_text(size=100), 
          legend.text = ggtext::element_markdown(size=100))
  #legend.key.size = unit(3, "cm"))
  
  
  
  
  ggsave(paste0(i, "_RLFS_with_POLR1A.png"), plot_rlfs_polr1a, width = 40, height = 31, dpi = 300)
  
  
  # correlation test
  cor_test <- cor.test(data_of_interest[[i]]$norm_POLR1A,
                       data_of_interest[[i]]$norm_RLFS_counts,
                       method = "pearson")
  
  # round values
  r_val <- round(cor_test$estimate, 1)
  p_val <- signif(cor_test$p.value, 3)
  
  
  # plot
  pearson_rlfs_polr1a<- ggplot(data_of_interest[[i]],
                               aes(x = norm_POLR1A, y = norm_RLFS_counts)) +
    geom_point(alpha = 0.6, size = 8, color = "black") +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    annotate("text", 
             x = 0.1, y = 0.9, hjust = 0,
             label = paste0("Pearson r = ", r_val," ", "P = ", p_val),
             size = 40) +   # text size scales differently than theme text
    labs(#title = paste0(i, "Correlation between POLR1A and RLFS"),
      x = "POLR1A ChIP-seq signal",
      y = "RLFS counts") +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    scale_y_continuous(limits = c(-0.25,1), breaks = seq(0,1,0.2)) +
    theme(#plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      text = element_text(size = 150, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 4),
      axis.ticks.length = unit(50, "pt"),
      axis.title.x = element_text(color = "#D2691E", margin = margin (t=20)), #size = 50),
      axis.title.y = element_text(color = "#aa2a85", hjust= 0.7, margin = margin (r=20)), #size = 50),
      axis.text.x = element_text(color = "#D2691E", hjust = 0.5), #size = 50),
      axis.text.y = element_text(color = "#aa2a85", hjust = 0.5), #size = 50),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 6),
      axis.line = element_blank())
  
  
  
  
  ggsave(paste0(i, "_RLFS_with_POLR1A_pearson.png"), pearson_rlfs_polr1a, width = 30, height = 20, dpi = 300)
  
  
  
  ############ plot G4FS
  
  cols <- c("G4FS" = "#228B22", "POLR1A" = "#D2691E")
  
  # legend-only data (NA coordinates => nothing is drawn on the panel)
  legend_df <- data.frame(
    bin_midpoints = NA_real_,
    y = NA_real_,
    grp = factor(c("POLR1A","G4FS"), levels = c("POLR1A","G4FS"))
  )
  
  plot_g4s_polr1a <-ggplot() +
    #geom_line(data = data, aes(x = bin_midpoints, y = UBF, color = "UBF"), size = 1.2) +
    geom_line(data = data_of_interest[[i]], aes(x = bin_midpoints, y = norm_POLR1A, color = "POLR1A"), size = 3.0,show.legend = FALSE) +
    
    geom_line(data = data_of_interest[[i]],
              aes(x = bin_midpoints, y = norm_G4FS_counts, color = "G4FS"), size = 3.0, show.legend = FALSE) +
    #geom_point(data = chip_with_nontemplate,
    #aes(x = bin_midpoints, y = norm_G4FS_counts, shape = "G4FS", color = "G4FS"), size = 3) +
    
    # legend-only filled squares
    geom_point(data = legend_df,
               aes(x = bin_midpoints, y = y, color = grp),
               shape = 15, size = 10, inherit.aes = FALSE, na.rm = TRUE) +
    
    scale_color_manual(name = NULL, 
                       values = cols,
                       labels = c("G4FS" = "G4FS", "POLR1A" = "POLR1A"))+
    #labels = c(
    #"G4FS"   = "<span style='color:#228B22;'>G4FS</span>",
    #"POLR1A" = "<span style='color:#D2691E;'>POLR1A</span>" #"UBF" = "#f609b4" #CD32A0, #9e277b
    #),
    
    guides(color = guide_legend(override.aes = list(shape = 15, size = 40, linetype = NA))) +
    
    #scale_shape_manual(name = "Non-canonical structure", values = c("G4FS" = 16)) +
    scale_x_continuous(breaks = region_marks, labels = region_labels) +
    scale_y_continuous(
      name = "POLR1A ChIP-seq signal",
      limits = c(0,1),
      breaks = seq(0, 1, 0.2),
      sec.axis = sec_axis(~., breaks = seq(0, 1, 0.2),
                          name = "G4FS counts")
    ) +
    labs(
      #title = paste0(i," G4FS and POLR1A"),
      x = "Mouse rDNA region"
    ) +
    theme(axis.ticks = element_line(color = "black", linewidth = 4),
          axis.ticks.length = unit(50, "pt"),
          #plot.title = element_text(hjust = 0.5, face = "bold", size=20),
          text = element_text(size = 150),
          panel.border = element_rect(color = "black",fill=NA, linewidth = 6),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),#if you want to add a rectangle box or you can use theme_minimal()
          axis.text.x = element_text(angle = 90, hjust = 1, size=30, color = "black", margin = margin(b=20)),
          axis.title.y.left  = element_text(color = "#D2691E", margin = margin (r=20)), #size = 50),  # left axis orange
          axis.title.y.right = element_text(color = "#228B22", margin = margin (l=20)), #size = 50),  # right axis purple
          axis.text.y.left   = element_text(color = "#D2691E"),# size = 50),
          axis.text.y.right  = element_text(color = "#228B22"), #size = 50),
          legend.position = c(0.05,0.98),
          legend.justification = c("left", "top"),
          legend.key = element_blank(),
          #legend.title = element_text(size=80), 
          legend.text = ggtext::element_markdown(size=100))
  
  
  
  
  ggsave(paste0(i, "_G4FS_with_POLR1A.png"), plot_g4s_polr1a, width = 40, height = 31, dpi = 300)
  
  
  
  
  
  
  ###################### correlation test
  cor_test <- cor.test(data_of_interest[[i]]$norm_POLR1A,
                       data_of_interest[[i]]$norm_G4FS_counts,
                       method = "pearson")
  
  # round values
  r_val <- round(cor_test$estimate, 1)
  p_val <- signif(cor_test$p.value, 3)
  
  
  # plot
  pearson_g4s_polr1a<- ggplot(data_of_interest[[i]],
                              aes(x = norm_POLR1A, y = norm_G4FS_counts)) +
    geom_point(alpha = 0.6, size = 8, color = "black") +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    annotate("text", 
             x = 0.1, y = 0.9, hjust = 0,
             label = paste0("Pearson r = ", r_val," ", "P = ", p_val),
             size = 40) +   # text size scales differently than theme text
    labs(#title = paste0(i, "Correlation between POLR1A and G4FS"),
      x = "POLR1A ChIP-seq signal",
      y = "G4FS counts") +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    scale_y_continuous(limits = c(-0.25,1), breaks = seq(0,1,0.2)) +
    theme(#plot.title = element_text(hjust = 0.5, face = "bold"),
      text = element_text(size = 150, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 4),
      axis.ticks.length = unit(50, "pt"),
      axis.title.x = element_text(color = "#D2691E", margin = margin (t=20)), #size = 50),
      axis.title.y = element_text(color = "#228B22", hjust= 0.7, margin = margin (r=20)), #size = 50),
      axis.text.x = element_text(color = "#D2691E", hjust = 0.5), #size = 50),
      axis.text.y = element_text(color = "#228B22", hjust = 0.5), #size = 50),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 6),
      axis.line = element_blank())
  
  
  
  ggsave(paste0(i, "_G4FS_with_POLR1A_pearson.png"), pearson_g4s_polr1a, width = 30, height = 20, dpi = 300)
  
  
  ###############plot imotif
  cols <- c("iMFS" = "#32A0CD", "POLR1A" = "#D2691E")
  
  # legend-only data (NA coordinates => nothing is drawn on the panel)
  legend_df <- data.frame(
    bin_midpoints = NA_real_,
    y = NA_real_,
    grp = factor(c("POLR1A","iMFS"), levels = c("POLR1A","iMFS"))
  )
  
  plot_imfs_polr1a <-
    ggplot() +
    #geom_line(data = data, aes(x = bin_midpoints, y = UBF, color = "UBF"), size = 2.0) +
    geom_line(data = data_of_interest[[i]], aes(x = bin_midpoints, y = norm_POLR1A, color = "POLR1A"), size = 3.0, show.legend = FALSE) +
    
    geom_line(data = data_of_interest[[i]], aes(x = bin_midpoints, y = norm_iMFS_counts, color = "iMFS"), size = 3.0,show.legend = FALSE) +
    #geom_point(data = chip_with_nontemplate,
    #aes(x = bin_midpoints, y = norm_iMFS_counts, shape = "iMFS", color = "iMFS"), size = 3) +
    
    # legend-only filled squares
    geom_point(data = legend_df,
               aes(x = bin_midpoints, y = y, color = grp),
               shape = 15, size = 10, inherit.aes = FALSE, na.rm = TRUE) +
    
    scale_color_manual(name = NULL,
                       values = cols)+
    #values = c("iMFS" = "#32A0CD","POLR1A" = "#D2691E"),
    #labels = c(
    #"iMFS"   = "<span style='color:#32A0CD;'>iMFS</span>",
    #"POLR1A" = "<span style='color:#D2691E;'>POLR1A</span>"
    #)) +#"UBF" = "#f609b4"
    guides(color = guide_legend(override.aes = list(shape = 15, size = 40, linetype = NA))) +
    
    #scale_shape_manual(name = "Non-canonical structure", values = c("iMFS" = 16)) +
    scale_x_continuous(breaks = region_marks, labels = region_labels) +
    scale_y_continuous(
      name = "POLR1A ChIP-seq signal",
      limits = c(0,1),
      breaks = seq(0, 1, 0.2),
      sec.axis = sec_axis(~., breaks = seq(0, 1, 0.2),
                          name = "iMFS counts")
    ) +
    labs(
      #title = paste0(i," iMFS and POLR1A"),
      x = "Mouse rDNA region"
    ) +
    theme(axis.ticks = element_line(color = "black", linewidth = 4),
          axis.ticks.length = unit(50, "pt"),
          #plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
          text = element_text(size = 150),
          panel.border = element_rect(color = "black",fill=NA, linewidth = 6),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),#if you want to add a rectangle box or you can use theme_minimal()
          axis.text.x = element_text(angle = 90, hjust = 1, size=30, color = "black", margin = margin(b=20)),
          axis.title.y.left  = element_text(color = "#D2691E", margin = margin (r=20)), #size = 50),  # left axis orange
          axis.title.y.right = element_text(color = "#32A0CD", margin = margin (l=20)), #size = 50),  # right axis purple
          axis.text.y.left   = element_text(color = "#D2691E"), #size = 50),
          axis.text.y.right  = element_text(color = "#32A0CD"), #size = 50),
          legend.position = c(0.05,0.98),
          legend.justification = c("left", "top"),
          legend.key = element_blank(),
          #legend.title = element_text(size=100), 
          legend.text = ggtext::element_markdown(size=100))
  
  
  
  ggsave(paste0(i, "_iMFS_with_POLR1A.png"), plot_imfs_polr1a, width = 40, height = 31, dpi = 300)
  
  
  
  
  
  
  ###################### correlation test
  cor_test <- cor.test(data_of_interest[[i]]$norm_POLR1A,
                       data_of_interest[[i]]$norm_iMFS_counts,
                       method = "pearson")
  
  # round values
  r_val <- round(cor_test$estimate, 1)
  p_val <- signif(cor_test$p.value, 3)
  
  
  # plot
  pearson_imfs_polr1a<- ggplot(data_of_interest[[i]],
                               aes(x = norm_POLR1A, y = norm_iMFS_counts)) +
    geom_point(alpha = 0.6, size = 3, color = "black") +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    annotate("text", 
             x = 0.1, y = 0.9, hjust = 0,
             label = paste0("Pearson r = ", r_val," ", "P = ", p_val),
             size = 40) +   # text size scales differently than theme text
    labs(#title = paste0(i, "Correlation between POLR1A and iMFS"),
      x = "POLR1A ChIP-seq signal",
      y = "iMFS counts") +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    scale_y_continuous(limits = c(-0.25,1), breaks = seq(0,1,0.2)) +
    theme(#plot.title = element_text(hjust = 0.5, face = "bold"),
      text = element_text(size = 150, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 4),
      axis.ticks.length = unit(50, "pt"),
      axis.title.x = element_text(color = "#D2691E", margin = margin (t=20)), #size = 50),
      axis.title.y = element_text(color = "#32A0CD", hjust= 0.7, margin = margin (r=20)), #size = 50),
      axis.text.x = element_text(color = "#D2691E", hjust = 0.5), #size = 50),
      axis.text.y = element_text(color = "#32A0CD", hjust = 0.5), #size = 50),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 6),
      axis.line = element_blank()
    )
  
  
  
  ggsave(paste0(i, "_iMFS_with_POLR1A_pearson.png"), pearson_imfs_polr1a, width = 30, height = 20, dpi = 300)
  
}

data_of_interest<- list(
  #Entire = Entire,
  Nontemplate = Nontemplate,
  Template = Template
)

for (i in names(data_of_interest)){
  print(i)
  print("RLFS")
  #data_of_interest[[i]][data_of_interest[[i]] == ""] <- NA
  proportion<- (sum(data_of_interest[[i]]$POLR1_RLFS_sign_match == 1, na.rm = TRUE)/99)*100 #99 because first entry of this column should be NA
  print(paste0(proportion, "% of the time, the directions are opposite."))
  
  
  print("G4FS")
  proportion2<- (sum(data_of_interest[[i]]$POLR1_G4FS_sign_match == 1, na.rm = TRUE)/99)*100
  print(paste0(proportion2, "% of the time, the directions are opposite."))
  
  
  
  print("iMFS")
  proportion3<- (sum(data_of_interest[[i]]$POLR1_iMFS_sign_match == 1, na.rm = TRUE)/99)*100
  print(paste0(proportion3, "% of the time, the directions are opposite."))
  
}


#[1] "Nontemplate"
#[1] "RLFS"
#[1] "87.8787878787879% of the time, the directions are opposite."
#[1] "G4FS"
#[1] "94.949494949495% of the time, the directions are opposite."
#[1] "iMFS"
#[1] "87.8787878787879% of the time, the directions are opposite."

#[1] "Template"
#[1] "RLFS"
#[1] "91.9191919191919% of the time, the directions are opposite."
#[1] "G4FS"
#[1] "91.9191919191919% of the time, the directions are opposite."
#[1] "iMFS"
#[1] "88.8888888888889% of the time, the directions are opposite."



data_of_interest2<- list(
  Nontemplate_zoom = Nontemplate_zoom,
  Template_zoom = Template_zoom
)

for (i in names(data_of_interest2)){
  data_of_interest2[[i]]$POLR1_RLFS_sign_match#[1]<- NA
  data_of_interest2[[i]]$POLR1_G4FS_sign_match#[1]<- NA
  data_of_interest2[[i]]$POLR1_iMFS_sign_match#[1]<- NA
  print(i)
  print("RLFS")
  proportion<- (sum(data_of_interest2[[i]]$POLR1_RLFS_sign_match == 1, na.rm = TRUE)/34)*100
  print(paste0(proportion, "% of the time, the directions are opposite."))
  
  
  print("G4FS")
  proportion2<- (sum(data_of_interest2[[i]]$POLR1_G4FS_sign_match == 1, na.rm = TRUE)/34)*100
  print(paste0(proportion2, "% of the time, the directions are opposite."))
  
  
  
  print("iMFS")
  proportion3<- (sum(data_of_interest2[[i]]$POLR1_iMFS_sign_match == 1, na.rm = TRUE)/34)*100
  print(paste0(proportion3, "% of the time, the directions are opposite."))
  
}



#[1] "Nontemplate_zoom"
#[1] "RLFS"
#[1] "73.5294117647059% of the time, the directions are opposite."
#[1] "G4FS"
#[1] "97.0588235294118% of the time, the directions are opposite."
#[1] "iMFS"
#[1] "79.4117647058823% of the time, the directions are opposite."

#[1] "Template_zoom"
#[1] "RLFS"
#[1] "88.2352941176471% of the time, the directions are opposite."
#[1] "G4FS"
#[1] "91.1764705882353% of the time, the directions are opposite."
#[1] "iMFS"
#[1] "100% of the time, the directions are opposite."



png("mouse_rdna_zoom_black_and_white.png", width = 15, height= 10, units= "in", res = 600)
##plotting begins
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=6307, end=22709))

kp <- plotKaryotype(genome=custom_genome, plot.type = 2)

kpRect(kp, chr = 'rDNA_locus', x0 = 6307, x1 =9339 , y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks 3000 bp of  promoter

kpRect(kp, chr = 'rDNA_locus', x0 = 9307, x1 = 13313 , y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks 5'ETS 

kpRect(kp, chr = 'rDNA_locus', x0 = 13314, x1 = 15183, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 18S

kpRect(kp, chr = 'rDNA_locus', x0 = 15184, x1 = 16183, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks ITS1

kpRect(kp, chr = 'rDNA_locus', x0 = 16184, x1 = 16340, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 5.8S

kpRect(kp, chr = 'rDNA_locus', x0 = 16341, x1 = 17428, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks ITS2

kpRect(kp, chr = 'rDNA_locus', x0 = 17429, x1 = 22158, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 28S

kpRect(kp, chr = 'rDNA_locus', x0 = 22159, x1 = 22709, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks 3'ETS

dev.off()















