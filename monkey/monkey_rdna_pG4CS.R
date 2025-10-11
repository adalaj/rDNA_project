# ------------------------------------------------------------------------------
# This code is part of paper: In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus.
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
# This R script processes G-quadruplex forming sequence (G4FS) predictions in the 
# monkey rDNA locus to assign each predicted G4FS to defined 
# rDNA subregions (5′ETS, 18S, ITS1, 5.8S, ITS2, 28S, 3′ETS, and IGS).
#

# Specifically, it:
#   1.Reads G4FS predictions generated using the G4 canonical finder 
#      (run in Python: g4_canonical_finder_3.11python.py).
#   2. Separates template and non-template strand G4FS predictions.
#   3. Plots G4FS distributions on a custom rDNA ideogram using karyoploteR.
#   4. Assigns each G4FS to rDNA subregions (5′ETS, 18S, ITS1, 5.8S, ITS2, 28S, 3′ETS)
#      based on strand-specific start coordinates.
#   5. Counts and normalizes G4FS occurrences across subregions, including junctions,
#      following the rule that a G4FS is assigned to the region where it initiates.
#   6. Generates summary CSVs and high-resolution bar plots showing:
#        - Overall normalized G4FS distribution
#        - Strandwise (template vs non-template) distributions
#        - Junction-inclusive and junction-excluded comparisons
#
# Inputs:
#   - monkey rDNA FASTA sequence (GenBank: KX061890)
#   - G4FS output text file from g4_canonical_finder_3.11python.py
#
# Outputs:
#   - Annotated CSVs of G4FS counts per rDNA region (with and without junctions)
#   - Strandwise G4FS summary CSV
#   - Publication-quality TIFF figures of normalized G4FS distributions
#
# ------------------------------------------------------------------------------


#open terminal
#(python2.7) Downloads % conda activate python3.11
#(python3.11) Downloads % python g4_canonical_finder_3.11python.py nontemplate_monkey_5ets_KX061890_and_NR_146166_3ets.fasta >output_G4FS_nontemplate_monkey_5ets_KX061890_and_NR_146166_3ets.txt  
#(python3.11) Downloads % python g4_canonical_finder_3.11python.py nontemplate_monkey_5ets_KX061890_3ets.fasta >output_G4FS_nontemplate_monkey_5ets_KX061890_3ets.txt      



setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/g4s_and_rdna/monkey")
library(data.table)
library(tidyverse)
library(karyoploteR)

KX061890<- fread("output_G4FS_nontemplate_monkey_5ets_KX061890_3ets.txt", sep = "\t", header = FALSE) #114
KX061890_and_NR_146166<- fread("output_G4FS_nontemplate_monkey_5ets_KX061890_and_NR_146166_3ets.txt", sep = "\t", header = FALSE) #114

#no difference in G4FS count or annotation

datasets<- list(
  KX061890 = KX061890,
  KX061890_and_NR_146166 = KX061890_and_NR_146166
)

for (i in names(datasets)){
  filename<- i
  
  entire_g4s_rdna<- datasets[[i]]
  colnames(entire_g4s_rdna) <- c("chr", "G4FS_start", "G4FS_end", "sequence", "name", "strand")
  
  
  entire_g4s_rdna$rDNA_region <- "junction"
  
  
  #strand specificity will matter here now we want to allocate RLFS to rdna region sub components  based on their direction
  # to do that, i will be creating a column that will have actual RLFS start based on strand specificity
  
  entire_g4s_rdna<- entire_g4s_rdna %>% mutate(actual_G4FS_start = ifelse(entire_g4s_rdna$strand == "+", G4FS_start, G4FS_end))
  entire_g4s_rdna<- entire_g4s_rdna %>% mutate(actual_G4FS_end = ifelse(entire_g4s_rdna$strand == "+", G4FS_end, G4FS_start))
  entire_g4s_rdna$G4FS_length<- abs(entire_g4s_rdna$G4FS_start-entire_g4s_rdna$G4FS_end)
  
  
  entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start >= 1] <- "5'ETS"
  
  entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 3640 ] <- "5'ETS and 18S junction"
  
  entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start >= 3641] <- "18S"
  entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 5508 ] <- "18S and ITS1 junction"
  
  
  entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start >= 5509   ] <- "ITS1"
  entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 6535] <- "ITS1 and 5.8S junction"
  
  entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start >= 6536 ] <- "5.8S"
  entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 6692 ] <- "5.8S and ITS2 junction"
  
  entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start >= 6693] <- "ITS2"
  entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 7863] <- "ITS2 and 28S junction"
  
  entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start >= 7864 ] <- "28S"
  entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 12648 ] <- "28S and 3'ETS junction"
  
  entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start >= 12649] <- "3'ETS"
  
  
  fwrite(entire_g4s_rdna, paste0("G4FS_",i, "_monkey_junctn_details.csv"), sep = ",")
}

#####################################below code has not executed################


entire_g4s_rdna_summary<- entire_g4s_rdna %>% group_by(rDNA_region) %>% count()
#rows order is lost 

sum(entire_g4s_rdna_summary$n)
#[1] 65

entire_g4s_rdna_summary[,1]

#1 28S       5'ETS       ITS1       ITS2       


new_rows<- data.table(rDNA_region = c( "5'ETS and 18S junction", "18S", "18S and ITS1 junction", "ITS1 and 5.8S junction", "5.8S", 
                                       "5.8S and ITS2 junction",  "ITS2 and 28S junction", 
                                       "28S and 3'ETS junction", "3'ETS"))

new_rows$n<- 0


entire_g4s_rdna_summary<- rbind(entire_g4s_rdna_summary, new_rows)
entire_g4s_rdna_summary$rDNA_region
#[1] "28S"                    "5'ETS"                 
#[3] "ITS1"                   "ITS2"                  
#[5] "5'ETS and 18S junction" "18S"                   
#[7] "18S and ITS1 junction"  "ITS1 and 5.8S junction"
#[9] "5.8S"                   "5.8S and ITS2 junction"
#[11] "ITS2 and 28S junction"  "28S and 3'ETS junction"
#[13] "3'ETS"


row.names(entire_g4s_rdna_summary)<- entire_g4s_rdna_summary$rDNA_region
names(entire_g4s_rdna_summary)[2] <- "G4FS_count"


entire_g4s_rdna_summary<- entire_g4s_rdna_summary %>% mutate(norm_G4FS_count = G4FS_count/sum(entire_g4s_rdna_summary$G4FS_count)) %>% 
  mutate(norm_G4FS_count= round(norm_G4FS_count, 2))

fwrite(entire_g4s_rdna_summary, "G4FS_KX061890_monkey_at_junctn_graphinput.csv", sep = ",")

entire_g4s_rdna_summary$rDNA_region <- factor(entire_g4s_rdna_summary$rDNA_region, 
                                              levels = c("5'ETS", "5'ETS and 18S junction", 
                                                         "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                         "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                         "28S and 3'ETS junction", "3'ETS" ))


g4s_norm<- ggplot(entire_g4s_rdna_summary, aes(x= rDNA_region, y = norm_G4FS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized G4FS distribution in the monkey rDNA locus", 
       x= "monkey rDNA region", 
       y= "Normalized G4FS count")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.70))+
  geom_text(aes(label= G4FS_count, vjust= -0.5, size= 50))+
  scale_fill_manual(values= c("#FDCCE5", "steelblue", "#D0B6FF","darkviolet", "#EF9B20","burlywood2", "#A0322B", 
                              "pink", "#FFCC17","aquamarine", "#E5FFB6","greenyellow", "#3B8CC4"))+
  #scale_fill_manual(values = combined_colors)+
  theme_minimal()+
  theme(axis.text.x = element_blank())+ 
  theme(panel.grid = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        theme(panel.grid = element_blank()),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +  # Center Y-axis title
  theme(axis.ticks.y = element_line(color = "black"))
#coord_flip()

ggsave( "Normalized_G4FS_distribution_in_monkey_rDNA_subcomponents_incld_junctn_AR.tiff", 
        plot = g4s_norm, width=18,height=10, dpi=150)



g4s_rdna_summary<- entire_g4s_rdna_summary[!grepl("junction", entire_g4s_rdna_summary$rDNA_region),]


#When you use coord_flip(), the order of the factor levels in rDNA_region determines how the categories are displayed along the flipped y-axis.
#By default, the first factor level appears at the bottom when flipped, and the last factor level appears at the top.

g4s_rdna_summary$rDNA_region <- factor(g4s_rdna_summary$rDNA_region, 
                                       levels = rev(c("5'ETS", "18S", "ITS1", "5.8S", 
                                                      "ITS2","28S", "3'ETS")))



#To reverse the order so that "Promoter" appears at the top when flipped, modify the levels of the factor like this


g4s_norm_nojuntn<- ggplot(g4s_rdna_summary, aes(x= rDNA_region, y = norm_G4FS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized G4FS distribution in the monkey rDNA locus", 
       x= "monkey rDNA region", 
       y= "Normalized G4FS count", 
       fill = "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.70))+
  geom_text(aes(label= G4FS_count, hjust= -1.0, vjust= 0.5, size= 50))+
  scale_fill_manual(values= rev(c("#FDCCE5", "#D0B6FF", "#EF9B20", "#A0322B", 
                                  "#FFCC17", "#E5FFB6", "#3B8CC4")))+
  #guides(fill = guide_legend(reverse = TRUE))
  theme_minimal()+
  theme(axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
        axis.ticks.x = element_line(color = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),   # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))+
  coord_flip()

ggsave("Normalized_G4FS_distribution_in_monkey_rDNA_subcomponents_after_rule.tiff", 
       plot = g4s_norm_nojuntn, width=18,height=10, dpi=150)



#to make template and non-template
entire_g4s_rdna_summary2<- entire_g4s_rdna %>% group_by(rDNA_region, strand) %>% count()
names(entire_g4s_rdna_summary2)[3] <- "G4FS_count"

new_rows<- data.table(rDNA_region = c("5'ETS","3'ETS","3'ETS","18S","18S","5.8S","5.8S"),
                      strand= c("+","+", "-","+","-", "+", "-"),
                      G4FS_count = c(0,0, 0,0,0,0,0))

entire_g4s_rdna_summary2<- rbind(entire_g4s_rdna_summary2, new_rows)
entire_g4s_rdna_summary2$rDNA_region <- factor(entire_g4s_rdna_summary2$rDNA_region, 
                                               levels = c("5'ETS", "18S", "ITS1", "5.8S", 
                                                          "ITS2","28S", "3'ETS"))

entire_g4s_rdna_summary2<- entire_g4s_rdna_summary2 %>% mutate(norm_G4FS_count = G4FS_count/sum(entire_g4s_rdna_summary2$G4FS_count)) %>% 
  mutate(norm_G4FS_count= round(norm_G4FS_count, 2))

fwrite(entire_g4s_rdna_summary2, "G4FS_KX061890_monkey_no_junctn_strandwise_AR_graphinput.csv")


g4s_strandwise<- ggplot(entire_g4s_rdna_summary2, aes(x= rDNA_region, y = norm_G4FS_count, fill= strand)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized G4FS strandwise distribution in the monkey rDNA locus", 
       x= "monkey rDNA region", 
       y= "Normalized G4FS count", 
       fill= "G4FS strand")+
  scale_y_continuous(breaks= seq(0, 0.40, by = 0.1), limits =c(0,0.40))+
  geom_text(aes(label= G4FS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  scale_fill_manual(values= c("+" = "#E21515", "-" = "#1414E1"), 
                    labels = c("+" = "Non-template", "-" = "Template"))+
  #scale_fill_manual(values = combined_colors)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        #panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black")) #facet_wrap (~strand or ~rDNA region doesnt look good)



ggsave( "Normalized_strandwise_G4FS_distribution_in_monkey_rDNA_subcomponents_AR.tiff", 
        plot = g4s_strandwise, width=18,height=10, dpi=150)




nontemplate<- entire_g4s_rdna_summary2 %>% filter(strand == "+")
nontemplate$rDNA_region <- factor(nontemplate$rDNA_region, 
                                  levels = c( "5'ETS", "18S", "ITS1", "5.8S", 
                                              "ITS2","28S", "3'ETS"))


template<- entire_g4s_rdna_summary2 %>% filter(strand == "-")
template$rDNA_region <- factor(template$rDNA_region, 
                               levels = c("5'ETS", "18S", "ITS1", "5.8S", 
                                          "ITS2","28S", "3'ETS"))



g4s_nontemplate <- ggplot(nontemplate, aes(x= rDNA_region, y = norm_G4FS_count, fill= rDNA_region)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized Non-template G4FS distribution in the monkey rDNA locus", 
       x= "monkey rDNA region", 
       y= "Normalized Non-template G4FS count", 
       fill= "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.40, by = 0.1), limits =c(0,0.40))+
  geom_text(aes(label= G4FS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  
  scale_fill_manual(values= c("#FDCCE5", "#D0B6FF", "#EF9B20", "#A0322B", 
                              "#FFCC17", "#E5FFB6", "#3B8CC4"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))


ggsave( "Normalized_nontemplate_G4FS_distribution_in_monkey_rDNA_subcomponents_AR.tiff", 
        plot = g4s_nontemplate, width=18,height=10, dpi=150)

g4s_template <- ggplot(template, aes(x= rDNA_region, y = norm_G4FS_count, fill= rDNA_region)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized Template G4FS distribution in the monkey rDNA locus", 
       x= "monkey rDNA region", 
       y= "Normalized Template G4FS count", 
       fill= "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.40, by = 0.1), limits =c(0,0.40))+
  geom_text(aes(label= G4FS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  
  scale_fill_manual(values= c( "#FDCCE5", "#D0B6FF", "#EF9B20", "#A0322B", 
                               "#FFCC17", "#E5FFB6", "#3B8CC4"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))


ggsave( "Normalized_template_G4FS_distribution_in_monkey_rDNA_subcomponents_AR.tiff", 
        plot = g4s_template, width=18,height=10, dpi=150)


library(BiocManager)
BiocManager::install("karyoploteR")
library(karyoploteR)


#read the G4FS that overlapped with rdna locus
entire_rdna<- fread("output_G4FS_KX061890_monkey_rDNA_2017.txt", sep = "\t", header = FALSE) #65
entire_rdna$V1= "rDNA_locus"
colnames(entire_rdna)<- c("chr", "start", "end", "sequence", "name", "strand")

##separate as per strand
entire_rdna6_nontemplate<- entire_rdna %>% filter(strand=="+") #29
#because in NCBI keep nontemplate sequence.


entire_rdna6_template<- entire_rdna %>% filter(strand=="-") #36

##plotting begins
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=11863))
kp <- plotKaryotype(genome=custom_genome, plot.type = 2)


kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 = 3640 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
kpRect(kp, chr = 'rDNA_locus', x0 = 3641, x1 = 3659, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
kpRect(kp, chr = 'rDNA_locus', x0 = 3660, x1 = 6189, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
kpRect(kp, chr = 'rDNA_locus', x0 = 6190, x1 = 6346 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
kpRect(kp, chr = 'rDNA_locus', x0 = 6347, x1 = 7079, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
kpRect(kp, chr = 'rDNA_locus', x0 = 7080 , x1 = 11520, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
kpRect(kp, chr = 'rDNA_locus', x0 = 11521, x1 = 11863, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS


kpPlotRegions(kp, data=entire_rdna6_template, col="#1414E1", r0= -0.5, r1= -1.3)
kpPlotRegions(kp, data=entire_rdna6_nontemplate, col="#E21515", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width








