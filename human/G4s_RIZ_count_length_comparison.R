#length doesnt make sense so will only stick to count. 
# Calculate RIZ vs G4 vs Imotif count and their coorelation



#lets make 100 bin size histogram for g4s, riz, imotif and see if they show similar trend.

#read the files
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/g4s_and_rdna/human/pG4CS_at_rdna_output/files")
entire_g4s_rdna <- fread("pG4CS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", header = TRUE, sep = ",") #210

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/output")
rlfs<- fread("RLFS_KY962518_added_3500nt_IGS_upstream_master_qmrlfs_table_after_rule.csv", header = TRUE, sep = ",") #195

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/imotif/output/files")
imotif<- fread("imotif_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", header = TRUE, sep = ",") #85

#prepare the identifiers

#1) range
#as total rdna length is 44838 (promoter to IGS)
range(entire_g4s_rdna$actual_pG4CS_start) #actual_start ensures adjust the confusion for strand orientation.
#[1]  2197 44857
range(rlfs$actual_RLFS_start)
#[1]  2267 43029

range(imotif$actual_imotif_start)
#[1]  2183 45865


#2) mean of length, count for pG4CS, RIZ, RLFS
g4_summary<- entire_g4s_rdna %>% group_by(rDNA_region) %>% count()
mean_pG4CS_count<- round(mean(g4_summary$n),2)
mean_pG4CS_length<- round(mean(entire_g4s_rdna$pG4CS_length),2)

rlfs_summary<- rlfs %>% group_by(rDNA_region) %>% count()
mean_riz_count<- round(mean(rlfs_summary$n),2)
mean_riz_length<- round(mean(rlfs$length_RIZ),2)
#mean_rlfs_length<- round(mean(rlfs$length_RLFS),2) well scale is super high so i wont be taking length of RLFS



#thinking to do simultaneously

bin_size <- c(70, 100)#,30, 50, 75, 100, 200, 500, 1000)


#combined raw counts for G4s and RIZ
for (i in bin_size){
  bin_size_new<- i +1
  bin_break <- seq(1, 46000, length.out= bin_size_new) # this step is important because it is ensuring consistency across dataset 44857 is the rdna locus
  #The start and end coordinates of the bins must be identical for both datasets, so that the corresponding bins in both datasets refer to the exact same genomic region.
  
  bin_width<- diff(bin_break)[1]
  
  g4s_hist<- hist(entire_g4s_rdna[["actual_pG4CS_start"]], breaks = bin_break, plot = FALSE)
  riz_hist <- hist(rlfs[["actual_RLFS_start"]], breaks = bin_break, plot = FALSE)
  imotif_hist <- hist(imotif[["actual_imotif_start"]], breaks = bin_break, plot = FALSE)
  
  combined_count<- data.frame(
    bin_lower_value= head(bin_break,-1),
    bin_upper_value=tail(bin_break,-1),
    bin_midpoints = (head(bin_break,-1) + tail(bin_break,-1))/2,
    pG4CS_counts = g4s_hist$counts,
    RIZ_counts = riz_hist$counts,
   imotif_counts = imotif_hist$counts
  )

  
  fwrite(combined_count, paste("graph_input_pG4CS_RIZ_imotif_count_bar_graph_", i,"bin.csv", sep=""), sep=",")
  
  
  
  #together
  g4s_riz_imotif_count<- ggplot() + 
    geom_bar(data = combined_count, 
             aes(x = bin_midpoints, y = pG4CS_counts, fill = "pG4CS_counts"), 
             stat = "identity", color = "#A0CD32", alpha = 0.7) + 
    geom_density(data = entire_g4s_rdna,
               aes(x = pG4CS_start, 
               y = ..density.. * length(entire_g4s_rdna$pG4CS_start) * bin_width),
              color = "#75E11E", size = 1, bw = 2000, kernel = "gaussian", show.legend = FALSE)+ ##B4F609"
    geom_bar(data = combined_count, 
             aes(x = bin_midpoints, y = RIZ_counts, fill = "RIZ_counts"), 
             stat = "identity", color = "#CD32A0", alpha = 0.7) + 
    geom_density(data = rlfs,
                 aes(x = start_RIZ, 
                     y = ..density.. * length(rlfs$start_RIZ) * bin_width),
                 color = "#f609b4", size = 1, bw = 2000, kernel = "gaussian", show.legend = FALSE)+ #"#89216B"
    geom_bar(data = combined_count, 
             aes(x = bin_midpoints, y = imotif_counts, fill = "imotif_counts"), 
             stat = "identity", color = "#32A0CD", alpha = 0.7) + 
    geom_density(data = imotif,
                 aes(x = beg, 
                     y = ..density.. * length(imotif$beg) * bin_width),
                 color =   "#09B4F6", size = 1, bw = 2000, kernel = "gaussian", show.legend = FALSE)+ ##216B89"
    
    
    scale_x_continuous(breaks = c(seq(0, 50000, by = 5000)), 
                       labels = c(seq(0, 50000, by = 5000))) +
    labs(title= "pG4CS vs RIZ vs imotif frequency distribution in human rDNA", 
         x= paste0("Human rDNA region with",i, " bins"), 
         y= "Frequency", 
         fill= "Non-canonical structures")+
    #scale_y_continuous(breaks= seq(0, 30, by = 10), limits =c(0,30))+
    scale_fill_manual(values= c("RIZ_counts" = "#CD32A0", "pG4CS_counts" = "#A0CD32", "imotif_counts"= "#32A0CD"),
                      labels = c("RIZ_counts" = "RIZ", "pG4CS_counts" = "pG4CS", "imotif_counts" = "imotif"))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.ticks.x = element_line(color = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
          axis.ticks.y = element_line(color = "black"),
          legend.key = element_rect(color = NA))
  
  ggsave(paste("pG4CS_RIZ_imotif_count_entire_human_rdna_in_", i, "bin.tiff", sep = ""), 
         plot = g4s_riz_imotif_count, width = 18, height = 10, dpi = 150)
  
  
  #only smooth line
  g4s_riz_imotif_denisty<- ggplot() +
    geom_density(data = entire_g4s_rdna,
                 aes(x = pG4CS_start,
                     y = ..density.. * length(entire_g4s_rdna$pG4CS_start) * bin_width,
                     color = "pG4CS"),
                 size = 1, bw = 2000, kernel = "gaussian") +
    
    geom_density(data = rlfs,
                 aes(x = start_RIZ,
                     y = ..density.. * length(rlfs$start_RIZ) * bin_width,
                     color = "RIZ"),
                 size = 1, bw = 2000, kernel = "gaussian") +
    
    geom_density(data = imotif,
                 aes(x = beg,
                     y = ..density.. * length(imotif$beg) * bin_width,
                     color = "imotif"),
                 size = 1, bw = 2000, kernel = "gaussian") +
    
    scale_color_manual(name = "Non-canonical structures",
                       values = c("pG4CS" = "#75E11E", 
                                  "RIZ" = "#f609b4", 
                                  "imotif" = "#09B4F6")) +
    
    scale_x_continuous(breaks = seq(0, 50000, by = 5000),
                       labels = seq(0, 50000, by = 5000)) +
    
    labs(title = "Non-canonical structures probability distribution in the human rDNA",
         x = paste0("Human rDNA region with ", i, " bins"),
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
  
  ggsave(paste("pG4CS_RIZ_imotif_density_entire_human_rdna_in_", i, "bin.tiff", sep = ""), 
         plot = g4s_riz_imotif_denisty, width = 18, height = 10, dpi = 150)
  
}

  
 ggplot() + 
    geom_smooth(data = combined_count, 
                aes(x = bin_midpoints, y = pG4CS_counts), 
                method = "loess", 
                color = "#A0CD32", 
                se = FALSE, 
                size = 1.5) +
    geom_smooth(data = combined_count, 
                aes(x = bin_midpoints, y = RIZ_counts), 
                method = "loess", 
                color = "#CD32A0", 
                se = FALSE, 
                size = 1.5) +
    geom_smooth(data = combined_count, 
                aes(x = bin_midpoints, y = imotif_counts), 
                method = "loess", 
                color = "#32A0CD", 
                se = FALSE, 
                size = 1.5) +
   geom_density(data= combined_count, 
                aes(x = Value, y = ..density.. * sum(combined_count$pG4CS_counts) * diff(bin_break)[1], color = "#A0CD32"), size = 1, kernel = "gaussian", bw = 0.09) +
  
  scale_x_continuous(breaks = c(seq(0, 46000, by = 4600), 44838), 
                     labels = c(seq(0, 46000, by = 4600), "")) +
    labs(title= "pG4CS vs RIZ vs imotif frequency distribution in human rDNA", 
         x= paste0("Human rDNA region with",i, " bins"), 
         y= "Frequency", 
         fill= "Non-canonical structures")+
    scale_y_continuous(breaks= seq(0, 30, by = 10), limits =c(0,30))+
    scale_fill_manual(values= c("RIZ_counts" = "#CD32A0", "pG4CS_counts" = "#A0CD32", "imotif_counts"= "#32A0CD"),
                      labels = c("RIZ_counts" = "RIZ", "pG4CS_counts" = "pG4CS", "imotif_counts" = "imotif"))+
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
 
 
 
 #density:
 ggplot() +
   geom_density(data = entire_g4s_rdna, aes(x = pG4CS_start), fill = "#A0CD32", alpha = 0.5) +
   geom_density(data = rlfs, aes(x = start_RIZ), fill = "#CD32A0", alpha = 0.5) +
   geom_density(data = imotif, aes(x = beg), fill = "#32A0CD", alpha = 0.5)
 
 
 
 #denity converted to counts and X-axis cropped
 ggplot() +
   # i-motif
   geom_density(data = imotif, 
                aes(x = beg, y = after_stat(count)), 
                fill = "#32A0CD", alpha = 0.4, adjust = 1.5) +
   
   # pG4CS
   geom_density(data = entire_g4s_rdna, 
                aes(x = pG4CS_start, y = after_stat(count)), 
                fill = "#A0CD32", alpha = 0.4, adjust = 1.5) +
   
   # RIZ
   geom_density(data = rlfs, 
                aes(x = start_RIZ, y = after_stat(count)), 
                fill = "#CD32A0", alpha = 0.4, adjust = 1.5) +
   
   labs(title = "Non-canonical Structure Enrichment (Counts)", 
        x = "Position in Human rDNA", 
        y = "Count") +
   coord_cartesian(xlim = c(0, 17000)) +
   theme_minimal(base_size = 18)
 
 library(corrplot)
 
 
 cor_matrix <- cor(combined_count[, c("pG4CS_counts", "RIZ_counts", "imotif_counts")],
                   method = "spearman")  # or method = "spearman" for rank correlation or skewed distrubution
 corrplot(
   cor_matrix,
   method = "squqre", 
   type = "upper", 
   tl.col = "black",      # Text label color
   tl.cex = 1.5,          # Text label size
   col = colorRampPalette(c("white", "red"))(200),  # Red gradient
   addCoef.col = "black", # Optional: adds correlation coefficients
   number.cex = 1.2       # Size of numbers (if added)
 )
 
 
 
# combined_count= 70 bins
#                   pG4CS_counts RIZ_counts imotif_counts
# pG4CS_counts     1.0000000  0.6937813     0.6016396
# RIZ_counts       0.6937813  1.0000000     0.4504837
# imotif_counts    0.6016396  0.4504837     1.0000000
 
 
 
 # combined_count= 100 bins
 #pG4CS_counts RIZ_counts imotif_counts
 #pG4CS_counts     1.0000000  0.6774687     0.6513185
 #RIZ_counts       0.6774687  1.0000000     0.5159703
 #imotif_counts    0.6513185  0.5159703     1.0000000
 
 
#When you change the number of bins, you're changing the resolution of your data:

#More bins (e.g., 100) → finer resolution, more sensitive to local variations, possibly more noise.

#Fewer bins (e.g., 70) → smoother patterns, less noise, but may hide local peaks/dips.

#G4 and RIZ: Correlation is robust, even stronger at lower resolution → might indicate a broader co-distribution.

#G4 and iMOTIF: Correlation drops as resolution decreases → may indicate more localized co-occurrence.

#If you're interested in broad trends across the rDNA (e.g., domains of co-enrichment), use 70 bins.
 
#If you're interested in finer, localized overlaps, and can handle noise, use 100 bins or more.


 
  #separate pg4s
  g4s_count_entire<- ggplot() + 
    geom_bar(data = combined_count[,c("bin_midpoints", "pG4CS_counts")], 
             aes(x = bin_midpoints, y = pG4CS_counts), 
             stat = "identity", fill = "cornflowerblue", alpha = 0.7) + 
  
  scale_x_continuous(breaks = c(seq(0, 45000, by = 5000), 44838), 
                     labels = c(seq(0, 45000, by = 5000), "")) +
    labs(title= "pG4CS frequency distribution in human rDNA", 
         x= paste0("Human rDNA region with",i, " bins"), 
         y= "pG4CS Frequency") +
    scale_y_continuous(breaks= seq(0, 30, by = 10), limits =c(0,30))+
    scale_fill_manual(values= c("pG4CS_counts" = "cornflowerblue"))+
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
  
  ggsave(paste("pG4CS_count_entire_human_rdna_in_", i, "bin.tiff", sep = ""), 
         plot = g4s_count_entire, width = 18, height = 10, dpi = 150)
  
  
  #just riz
  riz_count_entire<- ggplot() + 
    geom_bar(data = combined_count[,c("bin_midpoints", "RIZ_counts")], 
             aes(x = bin_midpoints, y = RIZ_counts), 
             stat = "identity", fill = "pink", alpha = 0.7) + 
    
    scale_x_continuous(breaks = c(seq(0, 45000, by = 5000), 44838), 
                       labels = c(seq(0, 45000, by = 5000), "")) +
    labs(title= "RIZ frequency distribution in human rDNA", 
         x= paste0("Human rDNA region with",i, " bins"), 
         y= "RIZ Frequency") +
    scale_y_continuous(breaks= seq(0, 30, by = 10), limits =c(0,30))+
    scale_fill_manual(values= c("pG4CS_counts" = "pink"))+
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
  
  ggsave(paste("RIZ_count_entire_human_rdna_in_", i, "bin.tiff", sep = ""), 
         plot = riz_count_entire, width = 18, height = 10, dpi = 150)
  
  
  
  #scatter plot to show how g4s and RIZ are showing similar trend 
  graph_max_value <- max(max(combined_count$pG4CS_counts), max(combined_count$RIZ_counts))+1
  
  g4_riz<- ggplot(combined_count, aes(x = pG4CS_counts, y = RIZ_counts)) +
    geom_point(color = "blue") +
    geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a trend line
    labs(title = "Scatter Plot of pG4CS vs RIZ Frequency in human rDNA",
         x = "pG4CS Frequency",
         y = "RIZ Frequency") +
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
  
  ggsave(paste("pG4CS_vs_RIZ_scatter_plot_human_rdna_in_", i, "bin.tiff", sep = ""), 
         plot = g4_riz, width = 18, height = 10, dpi = 150)
  
  
  riz_g4<- ggplot(combined_count, aes(x = RIZ_counts, y = pG4CS_counts)) +
    geom_point(color = "blue") +
    geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a trend line
    labs(title = "Scatter Plot of RIZ vs pG4CS in human rDNA",
         x = "RIZ Frequency",
         y = "pG4CS Frequency") +
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
  
  ggsave(paste("RIZ_vs_pG4CS_scatter_plot_human_rdna_in_", i, "bin.tiff", sep = ""), 
         plot = riz_g4, width = 18, height = 10, dpi = 150)
  
  print(paste("RIZ_vs_pg4CS_scatter_plot_correlation_in_", i, "bins", sep = ""))
  
  print(cor.test(combined_count$RIZ_counts, combined_count$pG4CS_counts, method = "spearman")) #or cor.test(combined_count$pG4CS_counts,combined_count$RIZ_counts, method = "spearman") the rho value is same
}



#[1] "RIZ_vs_pg4CS_scatter_plot_correlation_in_100bins in seq 1 to 45000

#Spearman's rank correlation rho

#data:  combined_count$RIZ_counts and combined_count$pG4CS_counts
#S = 65242, p-value = 1.871e-11
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.6085067 


#Spearman's rank correlation rho in case of 100 bins in seq 1 to 44857.

#data:  combined_count$pG4CS_counts and combined_count$RIZ_counts
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



{
  #moving it down as this is not needed as of now
  #for individual count and length graphs
  
  g4s_length <- entire_g4s_rdna %>%
    mutate(bin = cut(actual_pG4CS_start, breaks = bin_break, include.lowest = TRUE, right = FALSE)) %>%
    group_by(bin) %>%
    summarise(avg_length = round(mean(pG4CS_length),2), .groups = "drop")
  
  g4s_length$bin <- gsub("\\[|\\]|\\(|\\)", "",g4s_length$bin)
  g4s_length <- g4s_length %>%
    separate(bin, c('lower_limit', 'upper_limit'), sep = ",")
  g4s_length$lower_limit<- as.numeric(g4s_length$lower_limit)
  g4s_length$upper_limit<- as.numeric(g4s_length$upper_limit)
  g4s_length<- g4s_length %>% mutate(bin_midpoints = (lower_limit + upper_limit) /2)
  g4s_length$identifier<- "pG4CS"
  
  riz_length <- rlfs %>%
    mutate(bin = cut(actual_RLFS_start, breaks = bin_break, include.lowest = TRUE, right = FALSE)) %>%
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
  fwrite(combined_length, paste("graph_input_pG4CS_RIZ_length_bar_graph_", i,"bin.csv", sep=""), sep=",")
  
  
  #individual G4S length
  g4s_length_entire <- ggplot() +
    geom_bar(data = g4s_length,
             aes(x = bin_midpoints, y = avg_length),
             stat = "identity", fill = "cornflowerblue", alpha = 0.7) +
    
    scale_x_continuous(breaks = c(seq(0, 45000, by = 5000), 44838), 
                       labels = c(seq(0, 45000, by = 5000), "")) +
    geom_hline(yintercept = mean_pG4CS_length, color = "cornflowerblue", linetype = "dashed", size = 1)+
    
    labs(title= "pG4CS length distribution in human rDNA", 
         x= paste0("Human rDNA region with",i, " bins"), 
         y= "pG4CS average length")+
    
    scale_y_continuous(breaks= seq(0, 30, by = 10), limits =c(0,30))+
    #scale_fill_manual(values= c("pG4CS_counts" = "cornflowerblue"))+
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
  
  
  ggsave(paste("pG4CS_length_entire_human_rdna_in_", i, "bin.tiff", sep = ""), 
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
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
          axis.ticks.y = element_line(color = "black"))
  
  ggsave(paste("RIZ_length_entire_human_rdna_in_", i, "bin.tiff", sep = ""), 
         plot = riz_length_entire, width = 18, height = 10, dpi = 150)
  
  
  #together
  g4s_riz_length<- ggplot(combined_length, aes(x= bin_midpoints, y = avg_length, fill= identifier)) + 
    geom_bar(stat= "identity", position ="dodge", color = "black") +
    labs(title= "pG4CS vs RIZ length distribution in human rDNA", 
         x= paste0("Human rDNA region with",i," bins"), 
         y= "Average length", 
         fill = "Non-canonical structures")+
    scale_y_continuous(breaks= seq(0, 30, by = 10), limits =c(0,30))+
    scale_x_continuous(breaks = c(seq(0, 45000, by = 5000), 44838), 
                       labels = c(seq(0, 45000, by = 5000), ""))+
    geom_hline(yintercept = mean_pG4CS_length, color = "cornflowerblue", linetype = "dashed", size = 1)+
    geom_hline(yintercept = mean_riz_length, color = "pink", linetype = "dashed", size = 1)+
    
    
    #geom_text(aes(label= RLFS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
    scale_fill_manual(values= c(pG4CS = "cornflowerblue", RIZ = "pink"), #changed the non template and template colors
                      labels = c(pG4CS = "pG4CS", RIZ = "RIZ"))+
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
  
  
  ggsave(paste("pG4CS_vs_RIZ_length_entire_human_rdna_in_", i, "bin.tiff", sep = ""), 
         plot = g4s_riz_length, width = 18, height = 10, dpi = 150)
  
}





