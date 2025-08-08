#Idea is to see trend of UBF and polr1 in presence of R-loop, G4s, imotif

#step1: identity chrR region of interest that where i will count the read count of UBF and polr1

library(data.table)
library(tidyverse)

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

#plotting begins

plot_rlfs_with_chip <- function(data, label_suffix) {
  ggplot() +
    geom_line(data = data, aes(x = bin_midpoints, y = UBF, color = "UBF"), size = 1.2) +
    geom_line(data = data, aes(x = bin_midpoints, y = POLR1A, color = "POLR1A"), size = 1.2) +
    
    geom_line(data = data,
              aes(x = bin_midpoints, y = RLFS_counts * 1000),
              color = "black", size = 0.8) +
    geom_point(data = data,
               aes(x = bin_midpoints, y = RLFS_counts * 1000, shape = "RLFS"),
               color = "black", size = 3) +
    
    scale_color_manual(name = "Signal", values = c("UBF" = "#f609b4", "POLR1A" = "#75E11E")) +
    scale_shape_manual(name = "Non-canonical structure", values = c("RLFS" = 16)) +
    scale_x_continuous(breaks = seq(0, 50000, by = 5000)) +
    scale_y_continuous(
      name = "Signal Intensity (UBF / POLR1A)",
      breaks = c(0, 5000, 10000, 15000, 20000),
      sec.axis = sec_axis(~ . / 1000, 
                          breaks = c(0,4,10,15,20),
                          name = "Non-canonical structure Count")
    ) +
    labs(
      title = paste("RLFS (", label_suffix, "), UBF, and POLR1A ChIP signal over human rDNA", sep = ""),
      x = "Human rDNA region with 100 bins"
    ) +
  theme_minimal() +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 30),
      axis.text.y.left = element_text(size = 30),
      axis.text.y.right = element_text(size = 30),
      axis.title.y.left = element_text(size = 30),
      axis.title.y.right = element_text(size = 30),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      text = element_text(size = 16),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.position = "bottom"
    )
  
}




# Create plots
plot_rlfs_entire <- plot_rlfs_with_chip(chip_with_entire, "entire")
plot_rlfs_template <- plot_rlfs_with_chip(chip_with_template, "template")
plot_rlfs_nontemplate <- plot_rlfs_with_chip(chip_with_nontemplate, "nontemplate")

# Save if needed
ggsave("RLFS_entire_with_UBF_POLR1A.tiff", plot_rlfs_entire, width = 18, height = 10, dpi = 150)
ggsave("RLFS_template_with_UBF_POLR1A.tiff", plot_rlfs_template, width = 18, height = 10, dpi = 150)
ggsave("RLFS_nontemplate_with_UBF_POLR1A.tiff", plot_rlfs_nontemplate, width = 18, height = 10, dpi = 150)





############ zoom in version

zoom_plot_rlfs_with_chip <- function(zoom_data, label_suffix) {
  zoom_start <- 7000
  zoom_end <- 23000
  
  chip_zoom <- zoom_data %>% filter(bin_midpoints >= zoom_start & bin_midpoints <= zoom_end)
  
  region_marks <- c(7137, 9339, 12996, 14865, 15935, 16092, 17259, 22310,22671)
  region_labels <- c("Promoter", "5′ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3′ETS","IGS")
  
  ggplot() +
    # ChIP signal
    geom_line(data = chip_zoom, aes(x = bin_midpoints, y = UBF, color = "UBF"), size = 1.2) +
    geom_line(data = chip_zoom, aes(x = bin_midpoints, y = POLR1A, color = "POLR1A"), size = 1.2) +
    geom_line(data = chip_zoom,
              aes(x = bin_midpoints, y = RLFS_counts * 1000),
              color = "black", size = 0.8) +
    geom_point(data = chip_zoom,
               aes(x = bin_midpoints, y = RLFS_counts * 1000, shape = "RLFS"),
               color = "black", size = 3) +
    
    scale_color_manual(name = "Signal", values = c("UBF" = "#f609b4", "POLR1A" = "#75E11E")) +
    scale_shape_manual(name = "Non-canonical structure", values = c("RLFS" = 16)) +
    scale_x_continuous(breaks = region_marks, labels = region_labels) +
    scale_y_continuous(
      name = "Signal Intensity (UBF / POLR1A)",
      breaks = c(0, 5000, 10000, 15000, 20000),
      sec.axis = sec_axis(~ . / 1000, 
                          breaks = c(0,4,10,15,20),
                          name = "Non-canonical structure Count")
    ) +
    labs(
      title = paste("RLFS (", label_suffix, "), UBF, and POLR1A ChIP signal over human rDNA", sep = ""),
      x = "Human rDNA region with 100 bins"
    ) +
    theme_minimal() +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 30),
      axis.text.y.left = element_text(size = 30),
      axis.text.y.right = element_text(size = 30),
      axis.title.y.left = element_text(size = 30),
      axis.title.y.right = element_text(size = 30),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      text = element_text(size = 16),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.position = "bottom"
    )
  
}

# Create plots
zoom_plot_rlfs_entire <- zoom_plot_rlfs_with_chip(chip_with_entire, "entire")
zoom_plot_rlfs_template <- zoom_plot_rlfs_with_chip(chip_with_template, "template")
zoom_plot_rlfs_nontemplate <- zoom_plot_rlfs_with_chip(chip_with_nontemplate, "nontemplate")

# Save if needed
ggsave("RLFS_zoom_entire_with_UBF_POLR1A.tiff", zoom_plot_rlfs_entire, width = 18, height = 10, dpi = 150)
ggsave("RLFS_zoom_template_with_UBF_POLR1A.tiff", zoom_plot_rlfs_template, width = 18, height = 10, dpi = 150)
ggsave("RLFS_zoom_nontemplate_with_UBF_POLR1A.tiff", zoom_plot_rlfs_nontemplate, width = 18, height = 10, dpi = 150)




#####################################################################

#did the same for pG4CS

plot_pG4CS_with_chip <- function(data, label_suffix) {
  ggplot() +
    geom_line(data = data, aes(x = bin_midpoints, y = UBF, color = "UBF"), size = 1.2) +
    geom_line(data = data, aes(x = bin_midpoints, y = POLR1A, color = "POLR1A"), size = 1.2) +
    
    geom_line(data = data,
              aes(x = bin_midpoints, y = pG4CS_counts * 1000),
              color = "black", size = 0.8) +
    geom_point(data = data,
               aes(x = bin_midpoints, y = pG4CS_counts * 1000, shape = "pG4CS"),
               color = "black", size = 3) +
    
    scale_color_manual(name = "Signal", values = c("UBF" = "#f609b4", "POLR1A" = "#75E11E")) +
    scale_shape_manual(name = "Non-canonical structure", values = c("pG4CS" = 18)) + #change will change circle to solid diamond
    scale_x_continuous(breaks = seq(0, 50000, by = 5000)) +
    scale_y_continuous(
      name = "Signal Intensity (UBF / POLR1A)",
      breaks = c(0, 5000, 10000, 15000, 20000),
      sec.axis = sec_axis(~ . / 1000, 
                          breaks = c(0,4,10,15,20),
                          name = "Non-canonical structure Count")
    ) +
    labs(
      title = paste("pG4CS (", label_suffix, "), UBF, and POLR1A ChIP signal over human rDNA", sep = ""),
      x = "Human rDNA region with 100 bins"
    ) +
    theme_minimal() +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 30),
      axis.text.y.left = element_text(size = 30),
      axis.text.y.right = element_text(size = 30),
      axis.title.y.left = element_text(size = 30),
      axis.title.y.right = element_text(size = 30),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      text = element_text(size = 16),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.position = "bottom"
    )
  
}




# Create plots
plot_pG4CS_entire <- plot_pG4CS_with_chip(chip_with_entire, "entire")
plot_pG4CS_template <- plot_pG4CS_with_chip(chip_with_template, "template")
plot_pG4CS_nontemplate <- plot_pG4CS_with_chip(chip_with_nontemplate, "nontemplate")

# Save if needed
ggsave("pG4CS_entire_with_UBF_POLR1A.tiff", plot_pG4CS_entire, width = 18, height = 10, dpi = 150)
ggsave("pG4CS_template_with_UBF_POLR1A.tiff", plot_pG4CS_template, width = 18, height = 10, dpi = 150)
ggsave("pG4CS_nontemplate_with_UBF_POLR1A.tiff", plot_pG4CS_nontemplate, width = 18, height = 10, dpi = 150)





############ zoom in version

zoom_plot_pG4CS_with_chip <- function(zoom_data, label_suffix) {
  zoom_start <- 7000
  zoom_end <- 23000
  
  chip_zoom <- zoom_data %>% filter(bin_midpoints >= zoom_start & bin_midpoints <= zoom_end)
  
  region_marks <- c(7137, 9339, 12996, 14865, 15935, 16092, 17259, 22310,22671)
  region_labels <- c("Promoter", "5′ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3′ETS","IGS")
  
  ggplot() +
    # ChIP signal
    geom_line(data = chip_zoom, aes(x = bin_midpoints, y = UBF, color = "UBF"), size = 1.2) +
    geom_line(data = chip_zoom, aes(x = bin_midpoints, y = POLR1A, color = "POLR1A"), size = 1.2) +
    geom_line(data = chip_zoom,
              aes(x = bin_midpoints, y = pG4CS_counts * 1000),
              color = "black", size = 0.8) +
    geom_point(data = chip_zoom,
               aes(x = bin_midpoints, y = pG4CS_counts * 1000, shape = "pG4CS"),
               color = "black", size = 3) +
    
    scale_color_manual(name = "Signal", values = c("UBF" = "#f609b4", "POLR1A" = "#75E11E")) +
    scale_shape_manual(name = "Non-canonical structure", values = c("pG4CS" = 18)) +
    scale_x_continuous(breaks = region_marks, labels = region_labels) +
    scale_y_continuous(
      name = "Signal Intensity (UBF / POLR1A)",
      breaks = c(0, 5000, 10000, 15000, 20000),
      sec.axis = sec_axis(~ . / 1000, 
                          breaks = c(0,4,10,15,20),
                          name = "Non-canonical structure Count")
    ) +
    labs(
      title = paste("pG4CS (", label_suffix, "), UBF, and POLR1A ChIP signal over human rDNA", sep = ""),
      x = "Human rDNA region with 100 bins"
    ) +
    theme_minimal() +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 30),
      axis.text.y.left = element_text(size = 30),
      axis.text.y.right = element_text(size = 30),
      axis.title.y.left = element_text(size = 30),
      axis.title.y.right = element_text(size = 30),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      text = element_text(size = 16),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.position = "bottom"
    )
  
}

# Create plots
zoom_plot_pG4CS_entire <- zoom_plot_pG4CS_with_chip(chip_with_entire, "entire")
zoom_plot_pG4CS_template <- zoom_plot_pG4CS_with_chip(chip_with_template, "template")
zoom_plot_pG4CS_nontemplate <- zoom_plot_pG4CS_with_chip(chip_with_nontemplate, "nontemplate")

# Save if needed
ggsave("pG4CS_zoom_entire_with_UBF_POLR1A.tiff", zoom_plot_pG4CS_entire, width = 18, height = 10, dpi = 150)
ggsave("pG4CS_zoom_template_with_UBF_POLR1A.tiff", zoom_plot_pG4CS_template, width = 18, height = 10, dpi = 150)
ggsave("pG4CS_zoom_nontemplate_with_UBF_POLR1A.tiff", zoom_plot_pG4CS_nontemplate, width = 18, height = 10, dpi = 150)



#####################################################################
#did the same for imotif



plot_imotif_with_chip <- function(data, label_suffix) {
  ggplot() +
    geom_line(data = data, aes(x = bin_midpoints, y = UBF, color = "UBF"), size = 1.2) +
    geom_line(data = data, aes(x = bin_midpoints, y = POLR1A, color = "POLR1A"), size = 1.2) +
    
    geom_line(data = data,
              aes(x = bin_midpoints, y = imotif_counts * 1000),
              color = "black", size = 0.8) +
    geom_point(data = data,
               aes(x = bin_midpoints, y = imotif_counts * 1000, shape = "imotif"),
               color = "black", size = 3) +
    
    scale_color_manual(name = "Signal", values = c("UBF" = "#f609b4", "POLR1A" = "#75E11E")) +
    scale_shape_manual(name = "Non-canonical structure", values = c("imotif" = 15)) +
    scale_x_continuous(breaks = seq(0, 50000, by = 5000)) +
    scale_y_continuous(
      name = "Signal Intensity (UBF / POLR1A)",
      breaks = c(0, 5000, 10000, 15000, 20000),
      sec.axis = sec_axis(~ . / 1000, 
                          breaks = c(0,4,10,15,20),
                          name = "Non-canonical structure Count")
    ) +
    labs(
      title = paste("imotif (", label_suffix, "), UBF, and POLR1A ChIP signal over human rDNA", sep = ""),
      x = "Human rDNA region with 100 bins"
    ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 30),
      axis.text.y.left = element_text(size = 30),
      axis.text.y.right = element_text(size = 30),
      axis.title.y.left = element_text(size = 30),
      axis.title.y.right = element_text(size = 30),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      text = element_text(size = 16),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.position = "bottom"
    )
  
}




# Create plots
plot_imotif_entire <- plot_imotif_with_chip(chip_with_entire, "entire")
plot_imotif_template <- plot_imotif_with_chip(chip_with_template, "template")
plot_imotif_nontemplate <- plot_imotif_with_chip(chip_with_nontemplate, "nontemplate")

# Save if needed
ggsave("imotif_entire_with_UBF_POLR1A.tiff", plot_imotif_entire, width = 18, height = 10, dpi = 150)
ggsave("imotif_template_with_UBF_POLR1A.tiff", plot_imotif_template, width = 18, height = 10, dpi = 150)
ggsave("imotif_nontemplate_with_UBF_POLR1A.tiff", plot_imotif_nontemplate, width = 18, height = 10, dpi = 150)

############ zoom in version

zoom_plot_imotif_with_chip <- function(zoom_data, label_suffix) {
  zoom_start <- 7000
  zoom_end <- 23000
  
  chip_zoom <- zoom_data %>% filter(bin_midpoints >= zoom_start & bin_midpoints <= zoom_end)
  
  region_marks <- c(7137, 9339, 12996, 14865, 15935, 16092, 17259, 22310,22671)
  region_labels <- c("Promoter", "5′ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3′ETS","IGS")
  
  ggplot() +
    # ChIP signal
    geom_line(data = chip_zoom, aes(x = bin_midpoints, y = UBF, color = "UBF"), size = 1.2) +
    geom_line(data = chip_zoom, aes(x = bin_midpoints, y = POLR1A, color = "POLR1A"), size = 1.2) +
    geom_line(data = chip_zoom,
              aes(x = bin_midpoints, y = imotif_counts * 1000),
              color = "black", size = 0.8) +
    geom_point(data = chip_zoom,
               aes(x = bin_midpoints, y = imotif_counts * 1000, shape = "imotif"),
               color = "black", size = 3) +
    
    scale_color_manual(name = "Signal", values = c("UBF" = "#f609b4", "POLR1A" = "#75E11E")) +
    scale_shape_manual(name = "Non-canonical structure", values = c("imotif" = 15)) +
    scale_x_continuous(breaks = region_marks, labels = region_labels) +
    scale_y_continuous(
      name = "Signal Intensity (UBF / POLR1A)",
      breaks = c(0, 5000, 10000, 15000, 20000),
      sec.axis = sec_axis(~ . / 1000, 
                          breaks = c(0,4,10,15,20),
                          name = "Non-canonical structure Count")
    ) +
    labs(
      title = paste("imotif (", label_suffix, "), UBF, and POLR1A ChIP signal over human rDNA", sep = ""),
      x = "Human rDNA region with 100 bins"
    ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 30),
      axis.text.y.left = element_text(size = 30),
      axis.text.y.right = element_text(size = 30),
      axis.title.y.left = element_text(size = 30),
      axis.title.y.right = element_text(size = 30),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      text = element_text(size = 16),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.position = "bottom"
    )
  
}

# Create plots
zoom_plot_imotif_entire <- zoom_plot_imotif_with_chip(chip_with_entire, "entire")
zoom_plot_imotif_template <- zoom_plot_imotif_with_chip(chip_with_template, "template")
zoom_plot_imotif_nontemplate <- zoom_plot_imotif_with_chip(chip_with_nontemplate, "nontemplate")

# Save if needed
ggsave("imotif_zoom_entire_with_UBF_POLR1A.tiff", zoom_plot_imotif_entire, width = 18, height = 10, dpi = 150)
ggsave("imotif_zoom_template_with_UBF_POLR1A.tiff", zoom_plot_imotif_template, width = 18, height = 10, dpi = 150)
ggsave("imotif_zoom_nontemplate_with_UBF_POLR1A.tiff", zoom_plot_imotif_nontemplate, width = 18, height = 10, dpi = 150)



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
         x = "Human rDNA region with 100 bins",
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


geom_point(data = chip_zoom,
           aes(x = bin_midpoints, y = imotif_counts * 1000, shape = "imotif"),
           color = "black", size = 3) +


  








