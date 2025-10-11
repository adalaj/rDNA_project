# ------------------------------------------------------------------------------
# This code is part of paper: In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus.
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
#   Quantifies, compares, and visualizes the genome-wide distribution of
#   three non-canonical DNA structures — Rloop initiation Zone (RIZ), G-quadruplexes (G4FS),
#   and i-motifs (iMFS) — across the human rDNA locus.
#   Assesses their positional overlap, correlation, and co-occurrence.
#
# Major Steps:
#   1. Bin the rDNA locus (default 100 bp bins) and count RIZ, G4FS, and iMFS per bin.
#   2. Plot combined bar + density overlays showing positional distributions.
#   3. Generate correlation matrix (Pearson/Spearman) and Venn diagrams.
#   4. Produce scatter plots to assess G4FS–RIZ co-enrichment trends.
#   5. Create schematic rDNA maps (karyoploteR) for visualization reference.
#
# Input:
#   - G4FS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv
#   - RIZ_KY962518_added_3500nt_IGS_upstream_at_junctn_details_after_rule.csv
#   - imotif_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv
#
# Output:
#   CSV and figure files showing binned distributions, correlations, co-occurrence (Venn), and 
#  scatter relationships of RIZ, G4FS, and iMFS across the human rDNA locus (Fig 5 and Supplementary table 5)
# ------------------------------------------------------------------------------


#Load libraries
library(tidyverse)
library(data.table)
library(karyoploteR)
library(corrplot)
library(VennDiagram)

#read the files
entire_g4s_rdna <- fread("G4FS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", header = TRUE, sep = ",") #210
entire_g4s_rdna<- entire_g4s_rdna %>% mutate(new_start = actual_G4FS_start-1298)


riz<- fread("RIZ_KY962518_added_3500nt_IGS_upstream_at_junctn_details_after_rule.csv", header = TRUE, sep = ",") #286
riz<- riz %>% mutate(new_start = actual_RIZ_start-1298)

imotif<- fread("imotif_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", header = TRUE, sep = ",") #85
imotif<- imotif %>% mutate(new_start = actual_imotif_start-1298)

#
# mean of length, count for G4FS, RIZ, RIZ
g4_summary<- entire_g4s_rdna %>% group_by(rDNA_region, G4FS_length) %>% mean()


mean_G4FS_count<- round(mean(g4_summary$n),2)
mean_G4FS_length<- round(mean(entire_g4s_rdna$G4FS_length),2)

riz_summary<- riz %>% group_by(rDNA_region) %>% count()
mean_riz_count<- round(mean(riz_summary$n),2)
mean_riz_length<- round(mean(riz$length_RIZ),2)
#mean_riz_length<- round(mean(riz$length_riz),2) well scale is super high so i wont be taking length of riz



#thinking to do simultaneously

bin_size <- c(100)#,30, 50, 75, 100, 200, 500, 1000)


#combined raw counts for G4s and RIZ
for (i in bin_size){
  bin_size_new<- i +1
  bin_break <- seq(1299, 46137, length.out= bin_size_new) # this step is important because it is ensuring consistency across dataset 44857 is the rdna locus
  #The start and end coordinates of the bins must be identical for both datasets, so that the corresponding bins in both datasets refer to the exact same genomic region.
  
  bin_width<- diff(bin_break)[1]
  
  g4s_hist<- hist(entire_g4s_rdna[["actual_G4FS_start"]], breaks = bin_break, plot = FALSE)
  riz_hist <- hist(riz[["actual_RIZ_start"]], breaks = bin_break, plot = FALSE)
  imotif_hist <- hist(imotif[["actual_imotif_start"]], breaks = bin_break, plot = FALSE)
  
  combined_count<- data.frame(
    bin_lower_value= head(bin_break,-1),
    bin_upper_value=tail(bin_break,-1),
    bin_midpoints = (head(bin_break,-1) + tail(bin_break,-1))/2,
    G4FS_counts = g4s_hist$counts,
    RIZ_counts = riz_hist$counts,
    imotif_counts = imotif_hist$counts
  )
  
  

  fwrite(combined_count, paste("graph_input_G4FS_RIZ_imotif_count_bar_graph_", i,"bin.csv", sep=""), sep=",")
  
  
  combined_count<- fread("graph_input_G4FS_RIZ_imotif_count_bar_graph_100bin.csv", sep=",", header = TRUE)  
  bin_size_new<- i +1 #100 bins
  bin_break <- seq(1299, 46137, length.out= bin_size_new)
  bin_width<- diff(bin_break)[1]
  
  
  
  #together
  g4s_riz_imotif_count<- ggplot() + 
    geom_bar(data = combined_count, 
             aes(x = bin_midpoints, y = G4FS_counts, fill = "G4FS_counts"), 
             stat = "identity", color = "#228B22", alpha = 0.7) + 
    geom_density(data = entire_g4s_rdna,
                 aes(x = G4FS_start, 
                     y = ..density.. * length(entire_g4s_rdna$G4FS_start) * bin_width),
                 color = "#228B22", size = 2, bw = 2000, kernel = "gaussian", show.legend = FALSE)+ ##B4F609"
    geom_bar(data = combined_count, 
             aes(x = bin_midpoints, y = RIZ_counts, fill = "RIZ_counts"), 
             stat = "identity", color = "#aa2a85", alpha = 0.7) + 
    geom_density(data = riz,
                 aes(x = RIZ_start, 
                     y = ..density.. * length(riz$RIZ_start) * bin_width),
                 color = "#aa2a85", size = 2, bw = 2000, kernel = "gaussian", show.legend = FALSE)+ #"#89216B"
    geom_bar(data = combined_count, 
             aes(x = bin_midpoints, y = imotif_counts, fill = "imotif_counts"), 
             stat = "identity", color = "#32A0CD", alpha = 0.7) + 
    geom_density(data = imotif,
                 aes(x = beg, 
                     y = ..density.. * length(imotif$beg) * bin_width),
                 color =   "#32A0CD", size = 2, bw = 2000, kernel = "gaussian", show.legend = FALSE)+ ##216B89"
    
    scale_x_continuous(breaks = c(1299, 3501, 7158,9027,10097, 10254, 11421, 16472, 16833, 46137),
                       labels = c("", "5'ET", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ET", "IGS", ""))+ #changed ETS to ET so that all labels have same space
    
    labs(#title= "Non-canonical Structures", 
         x= "Human rDNA region", 
         y= "Frequency", 
         fill= NULL)+
    scale_y_continuous(breaks= seq(5, 25, by = 10), limits =c(0,25))+
    scale_fill_manual(values= c("RIZ_counts" = "#aa2a85", "G4FS_counts" = "#228B22", "imotif_counts"= "#32A0CD"),
                      breaks = c("RIZ_counts", "G4FS_counts", "imotif_counts"),
                      labels = c("RIZ_counts" = "RIZ", "G4FS_counts" = "G4FS", "imotif_counts" = "iMFS"))+
    theme_minimal()+
    theme(axis.ticks = element_line(color = "black", linewidth = 4),
          axis.ticks.length = unit(50, "pt"),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 100),
          axis.line = element_line(color = "black", linewidth = 4),
          axis.title.x = element_text(vjust = 0.9, hjust = 0.6, color = "black", margin = margin (t=20)),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, margin = margin (r=20)),  # Center Y-axis title
          axis.ticks.length.x = unit(20,"pt"),
          axis.text.x  = element_text(color = "black", size = 40),
          axis.text.y  = element_text(color = "black"),
          legend.key = element_rect(color = NA),
          legend.text = element_text(size=60),
          legend.key.size = unit(3, "cm"),
          legend.position = "top")
  
  ggsave(paste("G4FS_RIZ_imotif_count_entire_human_rdna_in_", i, "bin.png", sep = ""), 
         plot = g4s_riz_imotif_count, width = 18, height = 15, dpi = 300)
  
  

  
  
  
  
  combined_count_zoom <- combined_count %>% filter(bin_midpoints >= 1299 & bin_midpoints <= 16832)
  
  entire_g4s_rdna_zoom <- entire_g4s_rdna %>% filter(rDNA_region != "IGS")
  riz_zoom<- riz %>% filter(rDNA_region != "IGS")
  imotif_zoom<- imotif %>% filter(rDNA_region !="IGS")
  
  g4s_riz_imotif_count_zoom<- ggplot() + 
    geom_bar(data = combined_count_zoom, 
             aes(x = bin_midpoints, y = G4FS_counts, fill = "G4FS_counts"), 
             stat = "identity", color = "#228B22", alpha = 0.7) + 
    geom_density(data = entire_g4s_rdna_zoom,
                 aes(x = G4FS_start, 
                     y = ..density.. * length(entire_g4s_rdna_zoom$G4FS_start) * bin_width),
                 color = "#228B22", size = 2, bw = 2000, kernel = "gaussian", show.legend = FALSE)+ ##B4F609"
    geom_bar(data = combined_count_zoom, 
             aes(x = bin_midpoints, y = RIZ_counts, fill = "RIZ_counts"), 
             stat = "identity", color = "#aa2a85", alpha = 0.7) + 
    geom_density(data = riz_zoom,
                 aes(x = RIZ_start, 
                     y = ..density.. * length(riz_zoom$RIZ_start) * bin_width),
                 color = "#aa2a85", size = 2, bw = 2000, kernel = "gaussian", show.legend = FALSE)+ #"#89216B"
    geom_bar(data = combined_count_zoom, 
             aes(x = bin_midpoints, y = imotif_counts, fill = "imotif_counts"), 
             stat = "identity", color = "#32A0CD", alpha = 0.7) + 
    geom_density(data = imotif_zoom,
                 aes(x = beg, 
                     y = ..density.. * length(imotif_zoom$beg) * bin_width),
                 color =   "#32A0CD", size = 2, bw = 2000, kernel = "gaussian", show.legend = FALSE)+ ##216B89"
    
    
    scale_x_continuous(breaks = c(1299, 3501, 7158,9027,10097, 10254, 11421, 16472, 16832, 46137),
                       labels = c("", "5'ET", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ET", "IGS", ""))+ #changed ETS to ET so that all labels have same space
    
    labs(#title= "Non-canonical Structures",  
         x= "Human rDNA region", 
         y= "Frequency", 
         fill= NULL)+
    scale_y_continuous(breaks= seq(5, 25, by = 10), limits =c(0,25))+
    scale_fill_manual(values= c("RIZ_counts" = "#aa2a85", "G4FS_counts" = "#228B22", "imotif_counts"= "#32A0CD"),
                      breaks = c("RIZ_counts", "G4FS_counts", "imotif_counts"),
                      labels = c("RIZ_counts" = "RIZ", "G4FS_counts" = "G4FS", "imotif_counts" = "iMFS"))+
    theme_minimal()+
    theme(axis.ticks = element_line(color = "black", linewidth = 4),
          axis.ticks.length = unit(50, "pt"),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 100),
          axis.line = element_line(color = "black", linewidth = 4),
          axis.title.x = element_text(vjust = 0.9, hjust = 0.6, color = "black", margin = margin (t=20)),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5,margin = margin (r=20)),  # Center Y-axis title
          axis.ticks.length.x = unit(20,"pt"),
          axis.text.x = element_text(angle = 90, hjust=1, size = 40, colour = "black"),
          axis.text.y  = element_text(color = "black"),
          legend.key = element_rect(color = NA),
          legend.text = element_text(size=60),
          legend.key.size = unit(3, "cm"),
          legend.position = "top")
  
  ggsave(paste("G4FS_RIZ_imotif_count_entire_human_rdna_in_", i, "bin_zoom.png", sep = ""), 
         plot = g4s_riz_imotif_count_zoom, width = 18, height = 15, dpi = 300)
  
  
  
  
  #only smooth line
  g4s_riz_imotif_density<- ggplot() +
    geom_density(data = entire_g4s_rdna,
                 aes(x = G4FS_start,
                     y = ..density.. * length(entire_g4s_rdna$G4FS_start) * bin_width,
                     color = "G4FS"),
                 size = 2, bw = 2000, kernel = "gaussian") +
    
    geom_density(data = riz,
                 aes(x = RIZ_start,
                     y = ..density.. * length(riz$RIZ_start) * bin_width,
                     color = "RIZ"),
                 size = 2, bw = 2000, kernel = "gaussian") +
    
    geom_density(data = imotif,
                 aes(x = beg,
                     y = ..density.. * length(imotif$beg) * bin_width,
                     color = "imotif"),
                 size = 2, bw = 2000, kernel = "gaussian") +
    
    scale_color_manual(
      name = NULL,
      values = c("RIZ" = "#aa2a85", "G4FS" = "#228B22", "imotif" = "#32A0CD"),
      breaks = c("RIZ", "G4FS", "imotif"),
      labels = c(RIZ="RIZ", G4FS="G4FS", imotif="iMFS"),guide = guide_legend(
        override.aes = list(fill = c("#aa2a85", "#228B22", "#32A0CD"), 
                            linetype = 0)  # removes line in legend
      )
    )+
    scale_y_continuous(breaks= seq(2, 8, by = 2), limits =c(0,8))+
    scale_x_continuous(breaks = c(1299, 3501, 7158,9027,10097, 10254, 11421, 16472, 16833, 46137),
                       labels = c("", "5'ET", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ET", "IGS", ""))+ #changed ETS to ET so that all labels have same space
    
    labs(#title= "Non-canonical Structures", 
         x= "Human rDNA region", 
         y= "Bin Density (Smooth)", 
         fill= NULL) +
    
    theme_minimal() +
    theme(axis.ticks = element_line(color = "black", linewidth = 4),
          axis.ticks.length = unit(50, "pt"),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 100),
          axis.line = element_line(color = "black", linewidth = 4),
          axis.title.x = element_text(vjust = 0.9, hjust = 0.6, color = "black", margin = margin (t=20)),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5,margin = margin (r=20)),  # Center Y-axis title
          axis.ticks.length.x = unit(20,"pt"),
          axis.text.x = element_text(angle = 90, hjust=1, size = 40, colour = "black"),
          axis.text.y  = element_text(color = "black"),
          legend.key = element_rect(color = NA),
          legend.text = element_text(size=60),
          legend.key.size = unit(3, "cm"),
          legend.position = "top")
  
  

  ggsave(paste("G4FS_RIZ_imotif_density_entire_human_rdna_in_", i, "bin.png", sep = ""), 
         plot = g4s_riz_imotif_density, width = 18, height = 15, dpi = 300)
  
  
  

png("rdna_zoom_black_and_white.png", width = 15, height= 10, units= "in", res = 600)
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1299, end=46137))
kp <- plotKaryotype(genome=custom_genome, plot.type = 2)

kpRect(kp, chr = 'rDNA_locus', x0 = 1299, x1 =3500 , y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks last 1298bp from IGS representing previous rdna 

kpRect(kp, chr = 'rDNA_locus', x0 = 3501, x1 = 7157 , y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks 5'ETS (3501+(3657-1))
#3501+(3657-1) = 7157

kpRect(kp, chr = 'rDNA_locus', x0 = 7158, x1 = 9026, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 18S
#7158+ (1869-1) = 9026

kpRect(kp, chr = 'rDNA_locus', x0 = 9027, x1 = 10096, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks ITS1S
#9027+ (1070-1) = 10096 

kpRect(kp, chr = 'rDNA_locus', x0 = 10097, x1 = 10253, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 5.8S
#10097+ (157-1) = 10253

kpRect(kp, chr = 'rDNA_locus', x0 = 10254, x1 = 11420, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks ITS2
#10254+(1167-1) = 11420

kpRect(kp, chr = 'rDNA_locus', x0 = 11421, x1 = 16471, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 28S
#11421+(5051-1) = 16471

kpRect(kp, chr = 'rDNA_locus', x0 = 16472, x1 = 16832, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks 3'ETS
#16472+(361-1) = 16832

kpRect(kp, chr = 'rDNA_locus', x0 = 16833, x1 = 46137, y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks IGS
dev.off()

png("rdna_zoom_black_and_white_zoom.png", width = 15, height= 10, units= "in", res = 600)
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1299, end=16832))
kp <- plotKaryotype(genome=custom_genome, plot.type = 2)

kpRect(kp, chr = 'rDNA_locus', x0 = 1299, x1 =3500 , y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks last 1298bp from IGS representing previous rdna 

kpRect(kp, chr = 'rDNA_locus', x0 = 3501, x1 = 7157 , y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks 5'ETS (3501+(3657-1))
#3501+(3657-1) = 7157

kpRect(kp, chr = 'rDNA_locus', x0 = 7158, x1 = 9026, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 18S
#7158+ (1869-1) = 9026

kpRect(kp, chr = 'rDNA_locus', x0 = 9027, x1 = 10096, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks ITS1S
#9027+ (1070-1) = 10096 

kpRect(kp, chr = 'rDNA_locus', x0 = 10097, x1 = 10253, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 5.8S
#10097+ (157-1) = 10253

kpRect(kp, chr = 'rDNA_locus', x0 = 10254, x1 = 11420, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks ITS2
#10254+(1167-1) = 11420

kpRect(kp, chr = 'rDNA_locus', x0 = 11421, x1 = 16471, y0 = 0, y1 = 1, col = "black", data.panel = "ideogram", borders= NA) #marks 28S
#11421+(5051-1) = 16471

kpRect(kp, chr = 'rDNA_locus', x0 = 16472, x1 = 16832, y0 = 0, y1 = 1, col = "#EAEAEA", data.panel = "ideogram", borders= NA) #marks 3'ETS
#16472+(361-1) = 16832
dev.off()



cor_matrix <- cor(combined_count[, c("RIZ_counts", "G4FS_counts","imotif_counts")],
                  method = "pearson")  # or method = "spearman" for rank correlation or skewed distrubution
new_labels <- c("RIZ","G4FS", "iMFS")
dimnames(cor_matrix) <- list(new_labels, new_labels)  # row + col names
fwrite(cor_matrix, "RIZ_g4s_imfs_correlation_plot_graph_input.csv")


png("RIZ_g4s_imfs_correlation_plot.png", width = 14, height = 10, units = "in", res = 300)



corrplot(
  cor_matrix,
  method = "square", 
  type = "upper", 
  tl.col = c("#aa2a85", "#228B22","#32A0CD"),      # Text label color
  tl.cex = 3.5,          # Text label size
  col = colorRampPalette(c("white", "red"))(200),  # Red gradient
  addCoef.col = "black", # Optional: adds correlation coefficients
  number.cex = 3.5,       # Size of numbers (if added)
  cl.cex = 4.0,         # Increase scale (legend) text size
  cl.ratio = 0.3,        # Increase thickness of the color bar
  cl.align.text = "l", # c=center, r=right, l= left align text to ticks
)

dev.off()



# combined_count= 70 bins
#                   G4FS_counts RIZ_counts imotif_counts
# G4FS_counts     1.0000000  0.6937813     0.6016396
# RIZ_counts       0.6937813  1.0000000     0.4504837
# imotif_counts    0.6016396  0.4504837     1.0000000



# combined_count= 100 bins
#G4FS_counts RIZ_counts imotif_counts
#G4FS_counts     1.0000000  0.6774687     0.6513185
#RIZ_counts       0.6774687  1.0000000     0.5159703
#imotif_counts    0.6513185  0.5159703     1.0000000


#When you change the number of bins, you're changing the resolution of your data:

#More bins (e.g., 100) → finer resolution, more sensitive to local variations, possibly more noise.

#Fewer bins (e.g., 70) → smoother patterns, less noise, but may hide local peaks/dips.

#G4 and RIZ: Correlation is robust, even stronger at lower resolution → might indicate a broader co-distribution.

#G4 and iMOTIF: Correlation drops as resolution decreases → may indicate more localized co-occurrence.

#If you're interested in broad trends across the rDNA (e.g., domains of co-enrichment), use 70 bins.

#If you're interested in finer, localized overlaps, and can handle noise, use 100 bins or more.





###Venn diagram:
#Convert each count column into presence/absence (1 if count > 0, 0 if count = 0).
combined_count$RIZ_bin    <- as.integer(combined_count$RIZ_counts > 0)
combined_count$G4FS_bin     <- as.integer(combined_count$G4FS_counts > 0)
combined_count$iMFS_bin   <- as.integer(combined_count$imotif_counts > 0)
fwrite(combined_count, "graph_input_G4FS_RIZ_imotif_count_bar_graph_100bin.csv")


png("RIZ_g4s_imfs_venn_plot.png", width = 11, height = 10, units = "in", res = 300)

venn.plot <- venn.diagram(
  x = list(
    RIZ   = which(combined_count$RIZ_bin == 1),
    G4FS   = which(combined_count$G4FS_bin == 1),    
    iMFS  = which(combined_count$iMFS_bin == 1)
    
  ),
  filename = NULL,
  fill = c("#aa2a85",  "#228B22","#32A0CD"),
  alpha = 0.6,#the greater the number the higher is transparency
  cex = 6,
  cat.cex = 5,
  cat.col = c("#aa2a85", "#228B22", "#32A0CD") #will give name of structures in this color order
  #cat.pos = c(-20, 90, 20),   # angle of labels relative to circles, this didnt work
  #cat.dist = c(0.08, 0.08, 0.06)  # distance from the border# category label colors
)

grid.draw(venn.plot)

dev.off()



png("RIZ_g4s_imfs_venn_plot2.png", width = 11, height = 10, units = "in", res = 300)

venn.plot <- venn.diagram(
  x = list(
    RIZ   = which(combined_count$RIZ_bin == 1),
    G4FS = which(combined_count$G4FS_bin == 1),
    iMFS  = which(combined_count$iMFS_bin == 1)
  ),
  filename = NULL,
  fill = NA,   # no fill inside the circles
  col = c("#aa2a85", "#228B22", "#32A0CD"),  # outline colors
  lwd = 3,     # line width of circles
  alpha = 1,   # keep outline solid
  cex = 5,     # numbers inside
  cat.cex = 0  # suppress category labels
)

grid.draw(venn.plot)

dev.off()


#final
png("RIZ_g4s_imfs_venn_plot3.png", width = 11, height = 10, units = "in", res = 300)

venn.plot <- venn.diagram(
  x = list(
    RIZ   = which(combined_count$RIZ_bin == 1),
    G4FS   = which(combined_count$G4FS_bin == 1),    
    iMFS  = which(combined_count$iMFS_bin == 1)
    
  ),
  filename = NULL,
  fill = c("#aa2a85",  "#228B22","#32A0CD"),
  col=NA, #if you want black border then remove this piece of code
  alpha = 0.6,#the greater the number the higher is transparency
  cex = 6,
  cat.cex = 5,
  category.names = c("","", "") 
  #cat.col = c("#aa2a85", "#228B22", "#32A0CD") will give name of structures in this color order
  #cat.pos = c(-20, 90, 20),   # angle of labels relative to circles, this didnt work
  #cat.dist = c(0.08, 0.08, 0.06)  # distance from the border# category label colors
)


grid.draw(venn.plot)

dev.off()



#scatter plot to show how g4s and RIZ are showing similar trend 
graph_max_value <- max(max(combined_count$G4FS_counts), max(combined_count$RIZ_counts))+1

g4_riz<- ggplot(combined_count, aes(x = G4FS_counts, y = RIZ_counts)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a trend line
  labs(title = "Scatter Plot of G4FS vs RIZ Frequency in human rDNA",
       x = "G4FS Frequency",
       y = "RIZ Frequency") +
  scale_x_continuous(breaks= seq(0, graph_max_value, by = 3))+
  scale_y_continuous(breaks= seq(0, graph_max_value, by = 3))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.ticks.x = element_line(color = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 40),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))  

ggsave(paste("G4FS_vs_RIZ_scatter_plot_human_rdna_in_", i, "bin.tiff", sep = ""), 
       plot = g4_riz, width = 18, height = 10, dpi = 150)


riz_g4<- ggplot(combined_count, aes(x = RIZ_counts, y = G4FS_counts)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a trend line
  labs(title = "Scatter Plot of RIZ vs G4FS in human rDNA",
       x = "RIZ Frequency",
       y = "G4FS Frequency") +
  scale_x_continuous(breaks= seq(0, graph_max_value, by = 3))+
  scale_y_continuous(breaks= seq(0, graph_max_value, by = 3))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.ticks.x = element_line(color = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))  

ggsave(paste("RIZ_vs_G4FS_scatter_plot_human_rdna_in_", i, "bin.tiff", sep = ""), 
       plot = riz_g4, width = 18, height = 10, dpi = 150)

print(paste("RIZ_vs_G4FS_scatter_plot_correlation_in_", i, "bins", sep = ""))

print(cor.test(combined_count$RIZ_counts, combined_count$G4FS_counts, method = "spearman")) #or cor.test(combined_count$G4FS_counts,combined_count$RIZ_counts, method = "spearman") the rho value is same
}



#[1] "RIZ_vs_G4FS_scatter_plot_correlation_in_100bins in seq 1 to 45000

#Spearman's rank correlation rho

#data:  combined_count$RIZ_counts and combined_count$G4FS_counts
#S = 65242, p-value = 1.871e-11
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.6085067 


#Spearman's rank correlation rho in case of 100 bins in seq 1 to 44857.

#data:  combined_count$G4FS_counts and combined_count$RIZ_counts
#S = 60273, p-value = 9.062e-13
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
rho 
#0.638326 


##Spearman rank correlation [ρ (rho) or r]: This measures the strength and direction of the association between 2 ranked variables.




#scatter plot interpretation

#Mainly becuase on X-axis it is counting RIZ and on Y-axis it counting G4. 
#and each dot represents the bins. 
#Each dot in your scatter plot represents one genomic bin (e.g., 30 bp or 100 bp region, depending on your bin size).

#For example:

#If you see four dots at (0, y-value) on the x-axis (G4s count = 0), it means that four different bins have zero RIZs but different G4s counts on the y-axis.
#This tells you that some bins have G4s without RIZ.

#Similarly:
#If a dot appears at (x, 0) on the y-axis, it means that there are bins with RIZ but no G4s.
#If a dot is at (0,0), it means there are bins with neither G4s nor RIZ

#If many dots are clustered around (0,0), most bins have low or no G4s and RIZ.
#If dots are spread diagonally, there is a positive correlation (bins with high G4s also tend to have high RIZ).
#If dots are all over with no pattern, there’s no strong correlation between G4s and RIZ.


#On a scatter diagram, the closer the points lie to a straight line, the stronger the linear relationship between two variables. 
#To quantify the strength of the relationship, we can calculate the correlation coefficient.
#A value of the correlation coefficient close to +1 indicates a strong positive linear relationship (i.e. one variable increases with the other). 
#A value close to -1 indicates a strong negative linear relationship (i.e. one variable decreases as the other increases
#A value close to 0 indicates no linear relationship; however, there could be a nonlinear relationship (a hyperbola type) between the variables
# refer: https://pmc.ncbi.nlm.nih.gov/articles/PMC374386/


#Extra code not needed as this moment
{{

  #for individual count and length graphs
  
  g4s_length <- entire_g4s_rdna %>%
    mutate(bin = cut(actual_G4FS_start, breaks = bin_break, include.lowest = TRUE, right = FALSE)) %>%
    group_by(bin) %>%
    summarise(avg_length = round(mean(G4FS_length),2), .groups = "drop")
  
  g4s_length$bin <- gsub("\\[|\\]|\\(|\\)", "",g4s_length$bin)
  g4s_length <- g4s_length %>%
    separate(bin, c('lower_limit', 'upper_limit'), sep = ",")
  g4s_length$lower_limit<- as.numeric(g4s_length$lower_limit)
  g4s_length$upper_limit<- as.numeric(g4s_length$upper_limit)
  g4s_length<- g4s_length %>% mutate(bin_midpoints = (lower_limit + upper_limit) /2)
  g4s_length$identifier<- "G4FS"
  
  riz_length <- riz %>%
    mutate(bin = cut(actual_RIZ_start, breaks = bin_break, include.lowest = TRUE, right = FALSE)) %>%
    group_by(bin) %>%
    summarise(avg_length = round(mean(length_RIZ),2), .groups = "drop")
  
  riz_length$bin <- gsub("\\[|\\]|\\(|\\)", "", riz_length$bin)
  riz_length <- riz_length %>%
    separate(bin, c('lower_limit', 'upper_limit'), sep = ",")
  riz_length$lower_limit<- as.numeric(riz_length$lower_limit)
  riz_length$upper_limit<- as.numeric(riz_length$upper_limit)
  riz_length<- riz_length %>% mutate(bin_midpoints = (lower_limit + upper_limit) /2)
  riz_length$identifier<- "RIZ"
  
  
  combined_length<- rbind(g4s_length, riz_length)
  fwrite(combined_length, paste("graph_input_G4FS_RIZ_length_bar_graph_", i,"bin.csv", sep=""), sep=",")
  
  
  #individual G4S length
  g4s_length_entire <- ggplot() +
    geom_bar(data = g4s_length,
             aes(x = bin_midpoints, y = avg_length),
             stat = "identity", fill = "cornflowerblue", alpha = 0.7) +
    
    scale_x_continuous(breaks = c(seq(0, 45000, by = 5000), 44838), 
                       labels = c(seq(0, 45000, by = 5000), "")) +
    geom_hline(yintercept = mean_G4FS_length, color = "cornflowerblue", linetype = "dashed", size = 1)+
    
    labs(title= "G4FS length distribution in human rDNA", 
         x= paste0("Human rDNA region (",i, "100bins)"), 
         y= "G4FS average length")+
    
    scale_y_continuous(breaks= seq(0, 30, by = 10), limits =c(0,30))+
    #scale_fill_manual(values= c("G4FS_counts" = "cornflowerblue"))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.ticks.x = element_line(color = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 40),
          axis.line = element_line(color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
          axis.ticks.y = element_line(color = "black"))
  
  
  ggsave(paste("G4FS_length_entire_human_rdna_in_", i, "bin.tiff", sep = ""), 
         plot = g4s_length_entire, width = 18, height = 10, dpi = 150)
  
  
  #individual riz length
  riz_length_entire <- ggplot() +
    geom_bar(data = riz_length,
             aes(x = bin_midpoints, y = avg_length),
             stat = "identity", fill = "pink", alpha = 0.7) +
    
    scale_x_continuous(breaks = c(seq(0, 45000, by = 5000), 44838), 
                       labels = c(seq(0, 45000, by = 5000), "")) +
    geom_hline(yintercept = mean_riz_length, color = "pink", linetype = "dashed", size = 1)+
    
    labs(title= "RIZ length distribution in human rDNA", 
         x= paste0("Human rDNA region with",i," bins"), 
         y= "RIZ average length")+
    
    scale_y_continuous(breaks= seq(0, 30, by = 10), limits =c(0,30))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.ticks.x = element_line(color = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 40),
          axis.line = element_line(color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
          axis.ticks.y = element_line(color = "black"))
  
  ggsave(paste("RIZ_length_entire_human_rdna_in_", i, "bin.tiff", sep = ""), 
         plot = riz_length_entire, width = 18, height = 10, dpi = 150)
  
  
  #together
  g4s_riz_length<- ggplot(combined_length, aes(x= bin_midpoints, y = avg_length, fill= identifier)) + 
    geom_bar(stat= "identity", position ="dodge", color = "black") +
    labs(title= "G4FS vs RIZ length distribution in human rDNA", 
         x= paste0("Human rDNA region with",i," bins"), 
         y= "Average length", 
         fill = "Non-canonical structures")+
    scale_y_continuous(breaks= seq(0, 30, by = 10), limits =c(0,30))+
    scale_x_continuous(breaks = c(seq(0, 45000, by = 5000), 44838), 
                       labels = c(seq(0, 45000, by = 5000), ""))+
    geom_hline(yintercept = mean_G4FS_length, color = "cornflowerblue", linetype = "dashed", size = 1)+
    geom_hline(yintercept = mean_riz_length, color = "pink", linetype = "dashed", size = 1)+
    
    
    #geom_text(aes(label= riz_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
    scale_fill_manual(values= c(G4FS = "cornflowerblue", RIZ = "pink"), #changed the non template and template colors
                      labels = c(G4FS = "G4FS", RIZ = "RIZ"))+
    #scale_fill_manual(values = combined_colors)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
          axis.ticks.x = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
          axis.ticks.y = element_line(color = "black"))
  
  
  ggsave(paste("G4FS_vs_RIZ_length_entire_human_rdna_in_", i, "bin.tiff", sep = ""), 
         plot = g4s_riz_length, width = 18, height = 10, dpi = 150)
  
}

{g4s_riz_imotif_scale_density <- ggplot() +
  geom_density(data = entire_g4s_rdna,
               aes(x = G4FS_start, y = after_stat(scaled), color = "G4FS"),#after_stat(scaled) (formerly ..scaled..) gives a curve whose maximum is 1 for each group; perfect for visual shape comparison
               size = 1, bw = 2000, kernel = "gaussian") +
  geom_density(data = riz,
               aes(x = RIZ_start, y = after_stat(scaled), color = "RIZ"),
               size = 1, bw = 2000, kernel = "gaussian") +
  geom_density(data = imotif,
               aes(x = beg, y = after_stat(scaled), color = "imotif"),
               size = 1, bw = 2000, kernel = "gaussian") +
  scale_color_manual(
    name = "Non-canonical structures",
    values = c("RIZ" = "#aa2a85", "G4FS" = "#228B22", "imotif" = "#32A0CD"),
    breaks = c("RIZ", "G4FS", "imotif"),
    labels = c(RIZ="RIZ", G4FS="G4FS", imotif="iMFS")
  ) +
  scale_x_continuous(breaks = c(1299, 3501, 7158,9027,10097, 10254, 11421, 16472, 16833, 46137),
                     labels = c("", "5'ET", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ET", "IGS", ""))+ #changed ETS to ET so that all labels have same space
  
  labs(
    title = "G4FS vs RIZ vs imotif (scaled density, 0–1)",
    x = paste0("Human rDNA region (", i, " bins)"),
    y = "Scaled density (0–1)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.ticks.x = element_line(color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    text = element_text(size = 40),
    axis.line = element_line(color = "black"),
    axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
    axis.ticks.y = element_line(color = "black"),
    legend.key = element_rect(color = NA),
    legend.position = "top"
  )

ggsave(paste("G4FS_RIZ_imotif_scaled_density_entire_human_rdna_in_", i, "bin.png", sep = ""), 
       plot = g4s_riz_imotif_density, width = 15, height = 10, dpi = 300)


g4s_riz_imotif_density_zoom<- ggplot() +
  geom_density(data = entire_g4s_rdna_zoom,
               aes(x = G4FS_start,
                   y = ..density.. * length(entire_g4s_rdna_zoom$G4FS_start) * bin_width,
                   color = "G4FS"),
               size = 1, bw = 2000, kernel = "gaussian") +
  
  geom_density(data = riz_zoom,
               aes(x = RIZ_start,
                   y = ..density.. * length(riz_zoom$RIZ_start) * bin_width,
                   color = "RIZ"),
               size = 1, bw = 2000, kernel = "gaussian") +
  
  geom_density(data = imotif_zoom,
               aes(x = beg,
                   y = ..density.. * length(imotif_zoom$beg) * bin_width,
                   color = "imotif"),
               size = 1, bw = 2000, kernel = "gaussian") +
  
  scale_color_manual(
    name = "Non-canonical structures",
    values = c("RIZ" = "#aa2a85", "G4FS" = "#228B22", "imotif" = "#32A0CD"),
    breaks = c("RIZ", "G4FS", "imotif"),
    labels = c(RIZ="RIZ", G4FS="G4FS", imotif="iMFS"))+
  
  scale_x_continuous(breaks = c(1299, 3501, 7158,9027,10097, 10254, 11421, 16472, 16833, 46137),
                     labels = c("", "5'ET", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ET", "IGS", ""))+ #changed ETS to ET so that all labels have same space
  
  labs(title= "G4FS vs RIZ vs imotif frequency distribution in human rDNA", 
       x= paste0("Human rDNA region (",i, "bins)"), 
       y= "Bin Density (Smooth)", 
       fill= "Non-canonical structures") +
  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        axis.ticks.x = element_line(color = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 40),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"),
        legend.key = element_rect(color = NA),
        legend.position = "top")




ggsave(paste("G4FS_RIZ_imotif_density_entire_human_rdna_in_", i, "bin_zoom.png", sep = ""), 
       plot = g4s_riz_imotif_density_zoom, width = 15, height = 10, dpi = 300)

}
}



