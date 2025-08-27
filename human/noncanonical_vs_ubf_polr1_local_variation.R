#Idea is to see trend of UBF and polr1 in presence of R-loop, G4s, imotif

#step1: identity chrR region of interest that where i will count the read count of UBF and polr1

library(data.table)
library(tidyverse)
install.packages("ggtext")
library(ggtext)
library(karyoploteR)

#set working directory
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/IGV/human/IGV_files/bedfiles")

modified_rdna<- fread("hg38-rDNA_v1.0.bed", sep = "\t", header = FALSE)
#hg38-rDNA_v1.0.bed file contains the bed format of modified rdna used in 
#Vikram paper title: Construction and validation of customized genomes for human and mouse ribosomal DNA mapping    deposited in github

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
  bin_break <- seq(1, 44838, length.out= bin_size_new)
  bin_width<- diff(bin_break)[1]
  bin_data<- data.frame(
    bin_lower_value= round(head(bin_break,-1)), #rounded because it was decimel
    bin_upper_value= round(tail(bin_break,-1)),
    bin_midpoints = round((head(bin_break,-1) + tail(bin_break,-1))/2)
)
}


bin_data$chr<- "chrR"

entire_rdna_100bins<- bin_data %>% select(4,1,2)
fwrite(entire_rdna_100bins, "chrR_entire_rdna_100bins.bed", sep = "\t")

#also saving bin_data because it has binpoints that will be needed for plotting
fwrite(bin_data, "entire_rdna_100bins.csv")



#go to ssh wisconsin 
#get the outraw files in this region
#out raw files are now in bedfiles folder that has rawount in all these regions which i will be using for mapping all non-canonical structures.


#here it the outfile

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/IGV/human/IGV_files/bedfiles")

chip<- fread("ChIP_entire_rDNA_100bins_summary_raw_counts.tab", sep = "\t", header = TRUE)



###########################################################
#prepare entire, nontemplate, template RLFS, G4s, imotif  for making graphs
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/IGV/human/IGV_files/Bedgraphs/input_bedfiles")


# Simple unnamed list
list_of_names <- list("imotif", "RLFS", "pG4CS")

# Define binning parameters
bin_size <- 100
bin_size_new <- bin_size + 1
bin_break <- seq(1, 44838, length.out = bin_size_new)
bin_width <- diff(bin_break)[1] #448.37

# Prepare empty lists to store counts
all_counts_nontemplate <- list()
all_counts_template <- list()
all_counts_entire <- list()

# Loop through each structure
for (label in list_of_names) {
  
  # === Load and process non-template ===
  nt_file <- paste0(label, "_human_modified_nontemplate.bed")
  nontemplate <- fread(nt_file, sep = "\t", header = FALSE)
  nontemplate$actual_start <- nontemplate$V2
  nontemplate$actual_end <- nontemplate$V3
  nt_hist <- hist(nontemplate$actual_start, breaks = bin_break, plot = FALSE)
  all_counts_nontemplate[[label]] <- nt_hist$counts
  
  # === Load and process template ===
  t_file <- paste0(label, "_human_modified_template.bed")
  template <- fread(t_file, sep = "\t", header = FALSE)
  template$actual_start <- template$V3
  template$actual_end <- template$V2
  t_hist <- hist(template$actual_start, breaks = bin_break, plot = FALSE)
  all_counts_template[[label]] <- t_hist$counts
  
  # === Combine into "entire" ===
  entire <- rbind(nontemplate, template)
  colnames(entire) <- c("chr", "start", "end", "name", "score", "strand", "actual_start", "actual_end")
  assign(paste0(label, "_entire"), entire)
  
  entire_hist <- hist(entire$actual_start, breaks = bin_break, plot = FALSE)
  all_counts_entire[[label]] <- entire_hist$counts
}

# Shared bin info
bin_info <- data.frame(
  bin_lower_value = round(head(bin_break, -1)),
  bin_upper_value = round(tail(bin_break, -1)),
  bin_midpoints   = round((head(bin_break, -1) + tail(bin_break, -1)) / 2)
)

# Construct final data.frames
combined_count_nontemplate <- cbind(bin_info, data.frame(
  pG4CS_counts = all_counts_nontemplate[["pG4CS"]],
  RLFS_counts  = all_counts_nontemplate[["RLFS"]],
  imotif_counts = all_counts_nontemplate[["imotif"]]
))

combined_count_template <- cbind(bin_info, data.frame(
  pG4CS_counts = all_counts_template[["pG4CS"]],
  RLFS_counts  = all_counts_template[["RLFS"]],
  imotif_counts = all_counts_template[["imotif"]]
))

combined_count_entire <- cbind(bin_info, data.frame(
  pG4CS_counts = all_counts_entire[["pG4CS"]],
  RLFS_counts  = all_counts_entire[["RLFS"]],
  imotif_counts = all_counts_entire[["imotif"]]
))



# 1. Rename chip columns
colnames(chip) <- c("chr", "bin_start", "bin_end", "UBF", "POLR1A")

# 2. Merge chip with structural counts (column-wise, assuming same number of rows)
chip_with_entire <- cbind(chip, combined_count_entire[, c("bin_midpoints","pG4CS_counts", "RLFS_counts", "imotif_counts")])
chip_with_template <- cbind(chip, combined_count_template[, c("bin_midpoints","pG4CS_counts", "RLFS_counts", "imotif_counts")])
chip_with_nontemplate <- cbind(chip, combined_count_nontemplate[, c("bin_midpoints","pG4CS_counts", "RLFS_counts", "imotif_counts")])


####Save these files as graph input

fwrite(chip_with_entire, "Entire_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv")
fwrite(chip_with_nontemplate, "Nontemplate_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv")
fwrite(chip_with_template, "Template_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv")


###################################################


Entire<- fread("Entire_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv", sep = ",", header = TRUE)
Nontemplate<- fread("Nontemplate_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv", sep = ",", header = TRUE)
Template<- fread("Template_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv", sep = ",", header = TRUE)

data_of_interest<- list(
  Entire = Entire,
  Nontemplate = Nontemplate,
  Template = Template
)

for (i in names(data_of_interest)){
colnames(data_of_interest[[i]])[9]<- "iMFS_counts" #previously it was imotif_counts

#need to scale both the axis as between zero to 1 for all UBF, POLR1A, RLFS, pG4CS, iMFS
data_of_interest[[i]]<- data_of_interest[[i]] %>% mutate(norm_UBF= (UBF- min(UBF))/(max(UBF)-min(UBF)),
                                                         norm_POLR1A= (POLR1A- min(POLR1A))/(max(POLR1A)-min(POLR1A)),
                                                         norm_pG4CS_counts= (pG4CS_counts- min(pG4CS_counts))/(max(pG4CS_counts)-min(pG4CS_counts)),
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
    diff_norm_UBF = norm_UBF - lag(norm_UBF, default = 0),
    diff_norm_POLR1A = norm_POLR1A - lag(norm_POLR1A),
    diff_norm_pG4CS_counts = norm_pG4CS_counts - lag(norm_pG4CS_counts),
    diff_norm_RLFS_counts = norm_RLFS_counts - lag(norm_RLFS_counts),
    diff_norm_iMFS_counts = norm_iMFS_counts - lag(norm_iMFS_counts)
  )


#If one goes up and the other goes down (opposite signs), you mark it as 1 (a “mismatch” in direction)
data_of_interest[[i]] <- data_of_interest[[i]] %>%
  mutate(POLR1_RLFS_sign_match = ifelse(sign(diff_norm_POLR1A) != sign(diff_norm_RLFS_counts), 1, 0))

data_of_interest[[i]] <- data_of_interest[[i]] %>%
  mutate(POLR1_pG4CS_sign_match = ifelse(sign(diff_norm_POLR1A) != sign(diff_norm_pG4CS_counts), 1, 0))

data_of_interest[[i]] <- data_of_interest[[i]] %>%
  mutate(POLR1_iMFS_sign_match = ifelse(sign(diff_norm_POLR1A) != sign(diff_norm_iMFS_counts), 1, 0))


fwrite(data_of_interest[[i]], paste0(i, "_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv"))

}


#####read again
zoom_start <- 7000
zoom_end <- 23000

Entire<- fread("Entire_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv", sep = ",", header = TRUE)
Nontemplate<- fread("Nontemplate_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv", sep = ",", header = TRUE)
Template<- fread("Template_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv", sep = ",", header = TRUE)
Entire_zoom <- Entire %>% filter(bin_midpoints >= zoom_start & bin_midpoints <= zoom_end)
Nontemplate_zoom <- Nontemplate %>% filter(bin_midpoints >= zoom_start & bin_midpoints <= zoom_end)
Template_zoom <- Template %>% filter(bin_midpoints >= zoom_start & bin_midpoints <= zoom_end)



region_marks <- c(1, 7137, 9339, 12996, 14865, 15935, 16092, 17259, 22310,22671, 44838)
region_labels <- c("","Pro", "5′ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3′ETS","", "")


data_of_interest<- list(
  Entire = Entire,
  Entire_zoom = Entire_zoom,
  Nontemplate = Nontemplate,
  Nontemplate_zoom = Nontemplate_zoom,
  Template = Template,
  Template_zoom = Template_zoom
)


for (i in names(data_of_interest)){
  ##plot the RLFS with POLR1A
  plot_rlfs_polr1a <-
    ggplot() +
    #geom_line(data = data, aes(x = bin_midpoints, y = UBF, color = "UBF"), size = 1.2) +
    geom_line(data = data_of_interest[[i]], aes(x = bin_midpoints, y = norm_POLR1A, color = "POLR1A"), size = 1.2) +
    
    geom_line(data = data_of_interest[[i]],
              aes(x = bin_midpoints, y = norm_RLFS_counts, color = "RLFS"), size = 1.2) +
    #geom_point(data = chip_with_nontemplate,
    #aes(x = bin_midpoints, y = norm_RLFS_counts, shape = "RLFS", color = "RLFS"), size = 3) +
    
    scale_color_manual(name = NULL, 
                       values = c("RLFS" = "#aa2a85","POLR1A" = "#E68532"),
                       labels = c(
                         "RLFS"   = "<span style='color:#aa2a85;'>RLFS</span>",
                         "POLR1A" = "<span style='color:#E68532;'>POLR1A</span>"
                       )) +#"UBF" = "#f609b4" #CD32A0, #9e277b
    #scale_shape_manual(name = "Non-canonical structure", values = c("RLFS" = 16)) +
    scale_x_continuous(breaks = region_marks, labels = region_labels) +
    scale_y_continuous(
      name = "Normalized POLR1A ChIP-seq signal(0-1)",
      limits = c(0,1),
      breaks = seq(0, 1, 0.2),
      sec.axis = sec_axis(~., breaks = seq(0, 1, 0.2),
                          name = "Normalized RLFS counts(0-1)")
    ) +
    labs(
      title = paste0(i," RLFS and POLR1A"),
      x = "Human rDNA region (100 bins)"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 40),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),#if you want to add a rectangle box or you can use theme_minimal()
          axis.text.x = element_text(angle = 90, hjust = 1, size=20, color = "black"),
          axis.title.y.left  = element_text(color = "#E68532"),  # left axis orange
          axis.title.y.right = element_text(color = "#aa2a85"),  # right axis purple
          axis.text.y.left   = element_text(color = "#E68532"),
          axis.text.y.right  = element_text(color = "#aa2a85"),
          legend.position = "top", 
          legend.title = element_text(size=40), 
          legend.text = ggtext::element_markdown(size=40))
  
  ggsave(paste0(i, "_RLFS_with_POLR1A.png"), plot_rlfs_polr1a, width = 20, height = 14, dpi = 300)
  
  
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
    geom_point(alpha = 0.6, size = 3, color = "black") +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    annotate("text", 
             x = 0.1, y = 0.9, hjust = 0,
             label = paste0("Pearson r = ", r_val," ", 
                            "\n P = ", p_val),
             size = 12) +   # text size scales differently than theme text
    labs(title = paste0(i, "Correlation between POLR1A and RLFS"),
         x = "Normalized POLR1A ChIP-seq signal(0-1)",
         y = "Normalized RLFS counts(0-1)") +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    scale_y_continuous(limits = c(-0.25,1), breaks = seq(0,1,0.2)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 40, color = "black"),
          axis.title.x = element_text(color = "#E68532"),
          axis.title.y = element_text(color = "#aa2a85"),
          axis.text.x = element_text(color = "#E68532"),
          axis.text.y = element_text(color = "#aa2a85"),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          legend.position = "top", 
          legend.title = element_text(size = 40), 
          legend.text = element_text(size = 40))
  
  ggsave(paste0(i, "_RLFS_with_POLR1A_pearson.png"), pearson_rlfs_polr1a, width = 15, height = 14, dpi = 300)

  
  
  ############ plot pG4CS
  plot_g4s_polr1a <-
    ggplot() +
    #geom_line(data = data, aes(x = bin_midpoints, y = UBF, color = "UBF"), size = 1.2) +
    geom_line(data = data_of_interest[[i]], aes(x = bin_midpoints, y = norm_POLR1A, color = "POLR1A"), size = 1.2) +
    
    geom_line(data = data_of_interest[[i]],
              aes(x = bin_midpoints, y = norm_pG4CS_counts, color = "pG4CS"), size = 1.2) +
    #geom_point(data = chip_with_nontemplate,
    #aes(x = bin_midpoints, y = norm_pG4CS_counts, shape = "pG4CS", color = "pG4CS"), size = 3) +
    
    scale_color_manual(name = NULL, 
                       values = c("pG4CS" = "#228B22","POLR1A" = "#E68532"),
                       labels = c(
                         "pG4CS"   = "<span style='color:#228B22;'>pG4CS</span>",
                         "POLR1A" = "<span style='color:#E68532;'>POLR1A</span>"
                       )) +#"UBF" = "#f609b4"
    #scale_shape_manual(name = "Non-canonical structure", values = c("pG4CS" = 16)) +
    scale_x_continuous(breaks = region_marks, labels = region_labels) +
    scale_y_continuous(
      name = "POLR1A Signal Intensity",
      limits = c(0,1),
      breaks = seq(0, 1, 0.2),
      sec.axis = sec_axis(~., breaks = seq(0, 1, 0.2),
                          name = "pG4CS structure Count")
    ) +
    labs(
      title = paste0(i," pG4CS and POLR1A"),
      x = "Human rDNA region (100 bins)"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 40),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),#if you want to add a rectangle box or you can use theme_minimal()
          axis.text.x = element_text(angle = 90, hjust = 1, size=20, color = "black"),
          axis.title.y.left  = element_text(color = "#E68532"),  # left axis orange
          axis.title.y.right = element_text(color = "#228B22"),  # right axis purple
          axis.text.y.left   = element_text(color = "#E68532"),
          axis.text.y.right  = element_text(color = "#228B22"),
          legend.position = "top", 
          legend.title = element_text(size=40), 
          legend.text = ggtext::element_markdown(size=40))


  ggsave(paste0(i, "_pG4CS_with_POLR1A.png"), plot_g4s_polr1a, width = 20, height = 14, dpi = 300)
  
    
  
 
  
  
  ###################### correlation test
  cor_test <- cor.test(data_of_interest[[i]]$norm_POLR1A,
                       data_of_interest[[i]]$norm_pG4CS_counts,
                       method = "pearson")
  
  # round values
  r_val <- round(cor_test$estimate, 1)
  p_val <- signif(cor_test$p.value, 3)


  # plot
  pearson_g4s_polr1a<- ggplot(data_of_interest[[i]],
         aes(x = norm_POLR1A, y = norm_pG4CS_counts)) +
    geom_point(alpha = 0.6, size = 3, color = "black") +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    annotate("text", 
             x = 0.1, y = 0.9, hjust = 0,
             label = paste0("Pearson r = ", r_val," ", 
                            "\n P = ", p_val),
             size = 12) +   # text size scales differently than theme text
    labs(title = paste0(i, "Correlation between POLR1A and pG4CS"),
         x = "Normalized POLR1A ChIP-seq signal(0-1)",
         y = "Normalized pG4CS counts(0-1)") +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    scale_y_continuous(limits = c(-0.25,1), breaks = seq(0,1,0.2)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 40, color = "black"),
          axis.title.x = element_text(color = "#E68532"),
          axis.title.y = element_text(color = "#228B22"),
          axis.text.x = element_text(color = "#E68532"),
          axis.text.y = element_text(color = "#228B22"),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          legend.position = "top", 
          legend.title = element_text(size = 40), 
          legend.text = element_text(size = 40))
  
  ggsave(paste0(i, "_pG4CS_with_POLR1A_pearson.png"), pearson_g4s_polr1a, width = 15, height = 14, dpi = 300)
  
  
  ###############plot imotif
  plot_imfs_polr1a <-
    ggplot() +
    #geom_line(data = data, aes(x = bin_midpoints, y = UBF, color = "UBF"), size = 1.2) +
    geom_line(data = data_of_interest[[i]], aes(x = bin_midpoints, y = norm_POLR1A, color = "POLR1A"), size = 1.2) +
    
    geom_line(data = data_of_interest[[i]],
              aes(x = bin_midpoints, y = norm_iMFS_counts, color = "iMFS"), size = 1.2) +
    #geom_point(data = chip_with_nontemplate,
    #aes(x = bin_midpoints, y = norm_iMFS_counts, shape = "iMFS", color = "iMFS"), size = 3) +
    
    scale_color_manual(name = NULL, 
                       values = c("iMFS" = "#32A0CD","POLR1A" = "#E68532"),
                       labels = c(
                         "iMFS"   = "<span style='color:#32A0CD;'>iMFS</span>",
                         "POLR1A" = "<span style='color:#E68532;'>POLR1A</span>"
                       )) +#"UBF" = "#f609b4"
    #scale_shape_manual(name = "Non-canonical structure", values = c("iMFS" = 16)) +
    scale_x_continuous(breaks = region_marks, labels = region_labels) +
    scale_y_continuous(
      name = "POLR1A Signal Intensity",
      limits = c(0,1),
      breaks = seq(0, 1, 0.2),
      sec.axis = sec_axis(~., breaks = seq(0, 1, 0.2),
                          name = "iMFS structure Count")
    ) +
    labs(
      title = paste0(i," iMFS and POLR1A"),
      x = "Human rDNA region (100 bins)"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 40),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),#if you want to add a rectangle box or you can use theme_minimal()
          axis.text.x = element_text(angle = 90, hjust = 1, size=20, color = "black"),
          axis.title.y.left  = element_text(color = "#E68532"),  # left axis orange
          axis.title.y.right = element_text(color = "#32A0CD"),  # right axis purple
          axis.text.y.left   = element_text(color = "#E68532"),
          axis.text.y.right  = element_text(color = "#32A0CD"),
          legend.position = "top", 
          legend.title = element_text(size=40), 
          legend.text = ggtext::element_markdown(size=40))
  
  
  ggsave(paste0(i, "_iMFS_with_POLR1A.png"), plot_imfs_polr1a, width = 20, height = 14, dpi = 300)
  
  
  
  
  
  
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
             label = paste0("Pearson r = ", r_val," ", 
                            "\n P = ", p_val),
             size = 12) +   # text size scales differently than theme text
    labs(title = paste0(i, "Correlation between POLR1A and iMFS"),
         x = "Normalized POLR1A ChIP-seq signal(0-1)",
         y = "Normalized iMFS counts(0-1)") +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    scale_y_continuous(limits = c(-0.25,1), breaks = seq(0,1,0.2)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 40, color = "black"),
          axis.title.x = element_text(color = "#E68532"),
          axis.title.y = element_text(color = "#32A0CD"),
          axis.text.x = element_text(color = "#E68532"),
          axis.text.y = element_text(color = "#32A0CD"),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          legend.position = "top", 
          legend.title = element_text(size = 40), 
          legend.text = element_text(size = 40))
  
  ggsave(paste0(i, "_iMFS_with_POLR1A_pearson.png"), pearson_imfs_polr1a, width = 15, height = 14, dpi = 300)
  
}

data_of_interest1<- list(
  Entire = Entire,
  Nontemplate = Nontemplate,
  Template = Template
)

for (i in names(data_of_interest1)){
  print(i)
  print("RLFS")
  proportion<- (sum(data_of_interest[[i]]$POLR1_RLFS_sign_match == 1, na.rm = TRUE)/99)*100 #99 because first entry of this column should be NA
  print(paste0(proportion, "% of the time, the directions are opposite."))
  #88.89%  (polr1a vs rloop)
  #~88.9% of the time, the directions are opposite.
  #This is a strong indicator of an inverse relationship.
  
  
  print("pG4CS")
  proportion2<- (sum(data_of_interest[[i]]$POLR1_pG4CS_sign_match == 1, na.rm = TRUE)/99)*100
  print(paste0(proportion2, "% of the time, the directions are opposite."))
  
  
  
  print("iMFS")
  proportion3<- (sum(data_of_interest[[i]]$POLR1_iMFS_sign_match == 1, na.rm = TRUE)/99)*100
  print(paste0(proportion3, "% of the time, the directions are opposite."))
  
}

#[1] "Entire"
#[1] "RLFS"
#[1] "82.8282828282828% of the time, the directions are opposite."
#[1] "pG4CS"
#[1] "83.8383838383838% of the time, the directions are opposite."
#[1] "iMFS"
#[1] "72.7272727272727% of the time, the directions are opposite."



#[1] "Nontemplate"
#[1] "RLFS"
#[1] "88.8888888888889% of the time, the directions are opposite."
#[1] "pG4CS"
#[1] "93.9393939393939% of the time, the directions are opposite."
#[1] "iMFS"
#[1] "84.8484848484848% of the time, the directions are opposite."

#[1] "Template"
#[1] "RLFS"
#[1] "82.8282828282828% of the time, the directions are opposite."
#[1] "pG4CS"
#[1] "85.8585858585859% of the time, the directions are opposite."
#[1] "iMFS"
#[1] "84.8484848484848% of the time, the directions are opposite."




data_of_interest2<- list(
  Entire_zoom = Entire_zoom,
Nontemplate_zoom = Nontemplate_zoom,
  Template_zoom = Template_zoom
)

for (i in names(data_of_interest2)){
  data_of_interest2[[i]]$POLR1_RLFS_sign_match[1]<- NA
  data_of_interest2[[i]]$POLR1_pG4CS_sign_match[1]<- NA
  data_of_interest2[[i]]$POLR1_iMFS_sign_match[1]<- NA
  print(i)
  print("RLFS")
  proportion<- (sum(data_of_interest2[[i]]$POLR1_RLFS_sign_match == 1, na.rm = TRUE)/34)*100
  print(paste0(proportion, "% of the time, the directions are opposite."))
  #88.89%  (polr1a vs rloop)
  #~88.9% of the time, the directions are opposite.
  #This is a strong indicator of an inverse relationship.
  
  
  print("pG4CS")
  proportion2<- (sum(data_of_interest2[[i]]$POLR1_pG4CS_sign_match == 1, na.rm = TRUE)/34)*100
  print(paste0(proportion2, "% of the time, the directions are opposite."))
  
  
  
  print("iMFS")
  proportion3<- (sum(data_of_interest2[[i]]$POLR1_iMFS_sign_match == 1, na.rm = TRUE)/34)*100
  print(paste0(proportion3, "% of the time, the directions are opposite."))
  
}
  


#[1] "Entire_zoom"
#[1] "RLFS"
#[1] "76.4705882352941% of the time, the directions are opposite."
#[1] "pG4CS"
#[1] "85.2941176470588% of the time, the directions are opposite."
#[1] "iMFS"
#[1] "76.4705882352941% of the time, the directions are opposite."

#[1] "Nontemplate_zoom"
#[1] "RLFS"
#[1] "76.4705882352941% of the time, the directions are opposite."
#[1] "pG4CS"
#[1] "85.2941176470588% of the time, the directions are opposite."
#[1] "iMFS"
#[1] "82.3529411764706% of the time, the directions are opposite."

#[1] "Template_zoom"
#[1] "RLFS"
#[1] "76.4705882352941% of the time, the directions are opposite."
#[1] "pG4CS"
#[1] "91.1764705882353% of the time, the directions are opposite."
#[1] "iMFS"
#[1] "88.2352941176471% of the time, the directions are opposite."



png("rdna_black_and_white.png", width = 15, height= 10, units= "in", res = 600)
##plotting begins
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=44838))
kp <- plotKaryotype(genome=custom_genome, plot.type = 2)

kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =9339 , y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #IGS plus promoter

kpRect(kp, chr = 'rDNA_locus', x0 = 9340, x1 = 12996 , y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks 5'ETS (3501+(3657-1))
#3501+(3657-1) = 7157

kpRect(kp, chr = 'rDNA_locus', x0 = 12997, x1 = 14865, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 18S
#7158+ (1869-1) = 9026

kpRect(kp, chr = 'rDNA_locus', x0 = 14866, x1 = 15935, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks ITS1S
#9027+ (1070-1) = 10096 

kpRect(kp, chr = 'rDNA_locus', x0 = 15936, x1 = 16092, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 5.8S
#10097+ (157-1) = 10253

kpRect(kp, chr = 'rDNA_locus', x0 = 16093, x1 = 17259, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks ITS2
#10254+(1167-1) = 11420

kpRect(kp, chr = 'rDNA_locus', x0 = 17260, x1 = 22310, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 28S
#11421+(5051-1) = 16471

kpRect(kp, chr = 'rDNA_locus', x0 = 22311, x1 = 22671, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks 3'ETS
#16472+(361-1) = 16832

kpRect(kp, chr = 'rDNA_locus', x0 = 22672, x1 = 44838, y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks IGS
#16833+ (31506-1)= 48338
dev.off()


png("rdna_zoom_black_and_white.png", width = 15, height= 10, units= "in", res = 600)
##plotting begins
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=7137, end=22671))

kp <- plotKaryotype(genome=custom_genome, plot.type = 2)

kpRect(kp, chr = 'rDNA_locus', x0 = 7137, x1 =9339 , y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks 2202 bp of  promoter

kpRect(kp, chr = 'rDNA_locus', x0 = 9340, x1 = 12996 , y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks 5'ETS (3501+(3657-1))
#3501+(3657-1) = 7157

kpRect(kp, chr = 'rDNA_locus', x0 = 12997, x1 = 14865, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 18S
#7158+ (1869-1) = 9026

kpRect(kp, chr = 'rDNA_locus', x0 = 14866, x1 = 15935, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks ITS1S
#9027+ (1070-1) = 10096 

kpRect(kp, chr = 'rDNA_locus', x0 = 15936, x1 = 16092, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 5.8S
#10097+ (157-1) = 10253

kpRect(kp, chr = 'rDNA_locus', x0 = 16093, x1 = 17259, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks ITS2
#10254+(1167-1) = 11420

kpRect(kp, chr = 'rDNA_locus', x0 = 17260, x1 = 22310, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 28S
#11421+(5051-1) = 16471

kpRect(kp, chr = 'rDNA_locus', x0 = 22311, x1 = 22671, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks 3'ETS
#16472+(361-1) = 16832

dev.off()



cor.test(chip_with_entire$POLR1A, chip_with_entire$RLFS_counts, method = "spearman")

#Spearman's rank correlation rho

#data:  chip_with_entire$POLR1A and chip_with_entire$RLFS_counts
#S = 72996, p-value = 1.174e-09
#alternative hypothesis: true rho is not equal to 0
#sample estimates:rho 0.5619818 


#cor.test(chip_with_entire$UBF, chip_with_entire$RLFS_counts, method = "spearman")
#Spearman's rank correlation rho

#data:  chip_with_entire$UBF and chip_with_entire$RLFS_counts
#S = 77361, p-value = 9.233e-09
#alternative hypothesis: true rho is not equal to 0
#sample estimates:rho 0.5357884 

                     
                      

######shape information



#| Shape    | Description         | Shape Code |
#  | -------- | ------------------- | ---------- |
#  | ●        | Solid circle        | `16`       |
#  | ■        | Solid square        | `15`       |
#  | ▲        | Solid triangle up   | `17`       |
#  | ◀ / ▶    | Left/right triangle | `24`/`25`  |
#  | ◆        | Solid diamond       | `18`       |
#  | ⬣        | Plus (cross)        | `3`        |
#  | ×        | Cross               | `4`        |
#  | ⬤ (open) | Hollow circle       | `1`        |
#  | ☐ (open) | Hollow square       | `0`        |
  

##wanted to plot only G4s, imotif, and RLFS and do correlation
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/IGV/human/IGV_files/bedfiles")

list_of_files<- list(
  entire= "Entire_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv",
  Nontemplate = "Nontemplate_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv",
  Template = "Template_non_canonical_with_Chip_UBF_POLR1A_graph_input.csv"
)

for (i in names(list_of_files)) {
  
  data <- fread(list_of_files[[i]], sep = ",", header = TRUE)
  
  ggplot() +
    geom_density(data = data,
                 aes(x = bin_midpoints,
                     y = ..density.. * sum(data$pG4CS_counts) * bin_width/ sum(data$pG4CS_counts > 0),
                     weight = pG4CS_counts,
                     color = "pG4CS"),
                 size = 1, bw = 2000, kernel = "gaussian") +
    
    geom_density(data = data,
                 aes(x = bin_midpoints,
                     y = ..density.. * sum(data$RLFS_counts) * bin_width / sum(data$RLFS_counts > 0),
                     weight = RLFS_counts,
                     color = "RLFS"),
                 size = 1, bw = 2000, kernel = "gaussian") +
    
    geom_density(data = data,
                 aes(x = bin_midpoints,
                     y = ..density.. * sum(data$imotif_counts) * bin_width / sum(data$imotif_counts > 0),
                     weight = imotif_counts,
                     color = "imotif"),
                 size = 1, bw = 2000, kernel = "gaussian") +
    
    scale_color_manual(name = "Non-canonical structures",
                       values = c("pG4CS" = "#E3A81C", 
                                  "RLFS" = "#1CE3A8", 
                                  "imotif" = "#A81CE3")) +
    
    scale_x_continuous(breaks = seq(0, 50000, by = 5000),
                       labels = seq(0, 50000, by = 5000)) +
    
    labs(title = "Distribution of non-canonical structures in the human rDNA",
         x = "Human rDNA region (100 bins)",
         y = "Bin Density (Smooth)") +
    
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.x = element_line(color = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
          axis.ticks.y = element_line(color = "black"),
          legend.key = element_rect(color = NA))
  
  
  #print(p)
}










