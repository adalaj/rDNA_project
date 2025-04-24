# wanted to make new fig4: that will have 6 images 
#  i) RLFS length on Y axis present in entire rdna locus irrespective of region with two mean marks (RLFS, RIZ). 
#  ii) G4  length on Y axis present in entire rdna locus irrespective of region (RLFS,RIZ)
# iii) RIZ vs G4 coorelation
# iv, v, vi perfom the same with their count


#lets make 100 bin size histogram for g4s in riz and see if they show similar trend.

#read the files
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/g4s_and_rdna/human/pG4CS_at_rdna_output/files")
entire_g4s_rdna <- fread("pG4CS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", header = TRUE, sep = ",") #210

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/output")
rlfs<- fread("RLFS_KY962518_added_3500nt_IGS_upstream_master_qmrlfs_table_after_rule.csv", header = TRUE, sep = ",") #195


#prepare the identifiers

#1) range
#as total rdna length is 44838 (promoter to IGS)
range(entire_g4s_rdna$actual_pG4CS_start) #actual_start ensures adjust the confusion for strand orientation.
#[1]  2197 44857
range(rlfs$actual_RLFS_start)
#[1]  2267 43029


#2) mean of length, count for pG4CS, RIZ, RLFS
g4_summary<- entire_g4s_rdna %>% group_by(rDNA_region) %>% count()
mean_pG4CS_count<- round(mean(g4_summary$n),2)
mean_pG4CS_length<- round(mean(entire_g4s_rdna$pG4CS_length),2)

rlfs_summary<- rlfs %>% group_by(rDNA_region) %>% count()
mean_riz_count<- round(mean(rlfs_summary$n),2)
mean_riz_length<- round(mean(rlfs$length_RIZ),2)
#mean_rlfs_length<- round(mean(rlfs$length_RLFS),2) well scale is super high so i wont be taking length of RLFS



#thinking to do simultaneously

bin_size <- c(100)#,30, 50, 75, 100, 200, 500, 1000)


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


combined_count<- fread("graph_input_pG4CS_RIZ_count_bar_graph_100bin.csv", sep = ",", header = TRUE)

#combined raw counts for G4s and RIZ
for (i in bin_size){
  bin_size_new<- i +1
  bin_break <- seq(1, 45000, length.out= bin_size_new) # this step is important because it is ensuring consistency across dataset 44857 is the rdna locus
  #The start and end coordinates of the bins must be identical for both datasets, so that the corresponding bins in both datasets refer to the exact same genomic region.
  
  
  g4s_hist<- hist(entire_g4s_rdna[["actual_pG4CS_start"]], breaks = bin_break, plot = FALSE)
  riz_hist <- hist(rlfs[["actual_RLFS_start"]], breaks = bin_break, plot = FALSE)
  
  combined_count<- data.frame(
    bin_lower_value= head(bin_break,-1),
    bin_upper_value=tail(bin_break,-1),
    bin_midpoints = (head(bin_break,-1) + tail(bin_break,-1))/2,
    pG4CS_counts = g4s_hist$counts,
    RIZ_counts = riz_hist$counts
  )
  
  
  fwrite(combined_count, paste("graph_input_pG4CS_RIZ_count_bar_graph_", i,"bin.csv", sep=""), sep=",")
  
  #together
  g4s_riz_count<- ggplot() + 
    geom_bar(data = combined_count, 
             aes(x = bin_midpoints, y = pG4CS_counts, fill = "pG4CS_counts"), 
             stat = "identity", color = "cornflowerblue", alpha = 0.7) + 
    geom_bar(data = combined_count, 
             aes(x = bin_midpoints, y = RIZ_counts, fill = "RIZ_counts"), 
             stat = "identity", color = "pink", alpha = 0.7) + 
    #geom_point(data = combined_count, 
               #aes(x = bin_midpoints, y = pG4CS_counts), 
               #color = "cornflowerblue") + 
    #geom_smooth(data = combined_count, 
                #aes(x = bin_midpoints, y = pG4CS_counts), 
                #method = "lm", 
                #color = "cornflowerblue", 
               # se = FALSE) + 
    #geom_point(data = combined_count, 
               #aes(x = bin_midpoints, y = RIZ_counts), 
               #color = "pink") + 
    #geom_smooth(data = combined_count, 
                #aes(x = bin_midpoints, y = RIZ_counts), 
                #method = "lm", 
                #color = "pink", 
                #se = FALSE)+
    
    scale_x_continuous(breaks = c(seq(0, 45000, by = 5000), 44838), 
                       labels = c(seq(0, 45000, by = 5000), "")) +
    labs(title= "pG4CS vs RIZ frequency distribution in human rDNA", 
         x= paste0("Human rDNA region with",i, " bins"), 
         y= "Frequency", 
         fill= "Non-canonical structures")+
    scale_y_continuous(breaks= seq(0, 30, by = 10), limits =c(0,30))+
    scale_fill_manual(values= c("RIZ_counts" = "pink", "pG4CS_counts" = "cornflowerblue"),
                      labels = c("RIZ_counts" = "RIZ", "pG4CS_counts" = "pG4CS"))+
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
  
  ggsave(paste("pG4CS_RIZ_count_entire_human_rdna_in_", i, "bin.tiff", sep = ""), 
         plot = g4s_riz_count, width = 18, height = 10, dpi = 150)
  
  
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








