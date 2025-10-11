# ------------------------------------------------------------------------------
# This code accompanies the paper:
# "In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus"
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
# This R script processes and visualizes R-loop forming sequence (RLFS) predictions 
# across the monkey rDNA locus obtained from the QmRLFS-finder algorithm.
# It assigns each predicted RLFS to defined rDNA subregions and quantifies their 
# strand-specific and normalized distributions for both analytical and visual representation.
#
# Specifically, it:
#   1. Runs QmRLFS-finder (v1.5) on the monkey rDNA FASTA sequence 
#      (GenBank: KX061890) to identify RLFSs in BED format.
#   2. Parses and processes the RLFS output to calculate start–end coordinates, 
#      strand specificity, and RLFS lengths.
#   3. Assigns each RLFS to rDNA subcomponents (5′ETS, 18S, ITS1, 5.8S, ITS2, 28S, 3′ETS) 
#      using a rule-based system that counts RLFSs in the region where they initiate, 
#      ensuring junction-spanning RLFSs are not lost to boundary truncation.
#   4. Computes raw and normalized RLFS counts per region, including junctional regions, 
#      and generates tabular outputs suitable for statistical or graphical analyses.
#   5. Visualizes RLFS distributions using bar plots (overall, junction-excluded, 
#      and strandwise) and custom ideograms with karyoploteR for template and 
#      non-template strands.
#
# Inputs:
#   - monkey rDNA FASTA sequence (GenBank: KX061890)
#   - BED output file from QmRLFS-finder (v1.5)
#
# Outputs:
#   - Annotated CSVs summarizing RLFS counts per rDNA region (with and without junctions)
#   - Strandwise normalized RLFS count tables
#   - Publication-quality TIFF plots showing RLFS distributions
#   - KaryoploteR ideograms illustrating RLFS positions on both strands
#
# ------------------------------------------------------------------------------

library(data.table)
library(tidyverse)


#open terminal
#(base) % conda activate python2.7
#(python2.7)  Downloads % python QmRLFS-finder.py -bed -i nontemplate_monkey_5ets_KX061890_and_NR_146166_3ets.fasta -o nontemplate_monkey_5ets_KX061890_and_NR_146166_3ets_qmrlfs
#QmRLFS-finder.py (version v1.5)
#run on Wed Jun 11 2025 16:07:03 
#command line: python QmRLFS-finder.py -bed -i nontemplate_monkey_5ets_KX061890_and_NR_146166_3ets.fasta -o nontemplate_monkey_5ets_KX061890_and_NR_146166_3ets_qmrlfs

#Time used: 0.39 mins



#(python2.7)  Downloads % python QmRLFS-finder.py -bed -i nontemplate_monkey_5ets_KX061890_3ets.fasta -o nontemplate_monkey_5ets_KX061890_3ets_qmrlfs 
#QmRLFS-finder.py (version v1.5)
#run on Wed Jun 11 2025 16:14:35 
#command line: python QmRLFS-finder.py -bed -i nontemplate_monkey_5ets_KX061890_3ets.fasta -o nontemplate_monkey_5ets_KX061890_3ets_qmrlfs

#Time used: 0.40 mins

KX061890<- fread("nontemplate_monkey_5ets_KX061890_3ets_qmrlfs.out.bed", sep = "\t", header = FALSE) #131
KX061890_and_NR_146166<- fread("nontemplate_monkey_5ets_KX061890_and_NR_146166_3ets_qmrlfs.out.bed", sep = "\t", header = FALSE) #131

#no difference in RLFS count or annotation

datasets<- list(
  KX061890 = KX061890,
  KX061890_and_NR_146166 = KX061890_and_NR_146166
)

for (i in names(datasets)){
  filename<- i

entire_RLFSs_rdna<- datasets[[i]]
entire_RLFSs_rdna<- entire_RLFSs_rdna %>% select(1:6)
colnames(entire_RLFSs_rdna) <- c("GenBank_Accession", "RLFS_start", "RLFS_end", "RLFS_name","score", "strand")


entire_RLFSs_rdna$rDNA_region <- "junction"


#strand specificity will matter here now we want to allocate RLFS to rdna region sub components  based on their direction
# to do that, i will be creating a column that will have actual RLFS start based on strand specificity

entire_RLFSs_rdna<- entire_RLFSs_rdna %>% mutate(actual_rlfs_start = ifelse(entire_RLFSs_rdna$strand == "+", RLFS_start, RLFS_end))
entire_RLFSs_rdna<- entire_RLFSs_rdna %>% mutate(actual_rlfs_end = ifelse(entire_RLFSs_rdna$strand == "+", RLFS_end, RLFS_start))
entire_RLFSs_rdna$RLFS_length<- abs(entire_RLFSs_rdna$RLFS_start-entire_RLFSs_rdna$RLFS_end)


entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 1] <- "5'ETS"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 3640 ] <- "5'ETS and 18S junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 3641] <- "18S"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 5508 ] <- "18S and ITS1 junction"


entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 5509   ] <- "ITS1"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 6535] <- "ITS1 and 5.8S junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 6536 ] <- "5.8S"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 6692 ] <- "5.8S and ITS2 junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 6693] <- "ITS2"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 7863] <- "ITS2 and 28S junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 7864 ] <- "28S"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 12648 ] <- "28S and 3'ETS junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 12649] <- "3'ETS"


fwrite(entire_RLFSs_rdna, paste0("RLFS_",i, "_monkey_junctn_details.csv"), sep = ",")
}


#####################################below code has not executed################

entire_RLFSs_rdna_summary<- entire_RLFSs_rdna %>% group_by(rDNA_region) %>% count()
#rows order is lost 

sum(entire_RLFSs_rdna_summary$n)
#[1] 91

entire_RLFSs_rdna_summary[,1]

#1 18S         28S       3'ETS       
#4 5'ETS       ITS1       ITS2       


new_rows<- data.table(rDNA_region = c("5'ETS and 18S junction", "18S and ITS1 junction", "ITS1 and 5.8S junction", "5.8S", 
                                      "5.8S and ITS2 junction",  "ITS2 and 28S junction",
                                      "28S and 3'ETS junction" ))

new_rows$n<- 0


entire_RLFSs_rdna_summary<- rbind(entire_RLFSs_rdna_summary, new_rows)
entire_RLFSs_rdna_summary$rDNA_region
#[1] "18S"                    "28S"                   
#[3] "3'ETS"                  "5'ETS"                 
#[5] "ITS1"                   "ITS2"                  
#[7] "5'ETS and 18S junction" "18S and ITS1 junction" 
#[9] "ITS1 and 5.8S junction" "5.8S"                  
#[11] "5.8S and ITS2 junction" "ITS2 and 28S junction" 
#[13] "28S and 3'ETS junction"


row.names(entire_RLFSs_rdna_summary)<- entire_RLFSs_rdna_summary$rDNA_region
names(entire_RLFSs_rdna_summary)[2] <- "RLFS_count"


entire_RLFSs_rdna_summary<- entire_RLFSs_rdna_summary %>% mutate(norm_RLFS_count = RLFS_count/sum(entire_RLFSs_rdna_summary$RLFS_count)) %>% 
  mutate(norm_RLFS_count= round(norm_RLFS_count, 2))

fwrite(entire_RLFSs_rdna_summary, "RLFS_KX061890_monkey_at_junctn_graphinput.csv", sep = ",")

entire_RLFSs_rdna_summary$rDNA_region <- factor(entire_RLFSs_rdna_summary$rDNA_region, 
                                                levels = c("5'ETS", "5'ETS and 18S junction", 
                                                           "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                           "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                           "28S and 3'ETS junction", "3'ETS" ))


RLFS_norm<- ggplot(entire_RLFSs_rdna_summary, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized RLFS distribution in the monkey rDNA locus", 
       x= "monkey rDNA region", 
       y= "Normalized RLFS count")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.60))+
  geom_text(aes(label= RLFS_count, vjust= -0.5, size= 50))+
  scale_fill_manual(values= c( "#FDCCE5", "steelblue", "#D0B6FF","darkviolet", "#EF9B20","burlywood2", "#A0322B", 
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

ggsave( "Normalized_RLFS_distribution_in_monkey_rDNA_subcomponents_incld_junctn_AR.tiff", 
        plot = RLFS_norm, width=18,height=10, dpi=150)






RLFSs_rdna_summary<- entire_RLFSs_rdna_summary[!grepl("junction", entire_RLFSs_rdna_summary$rDNA_region),]


#When you use coord_flip(), the order of the factor levels in rDNA_region determines how the categories are displayed along the flipped y-axis.
#By default, the first factor level appears at the bottom when flipped, and the last factor level appears at the top.

RLFSs_rdna_summary$rDNA_region <- factor(RLFSs_rdna_summary$rDNA_region, 
                                         levels = rev(c("5'ETS", "18S", "ITS1", "5.8S", 
                                                        "ITS2","28S", "3'ETS")))



#To reverse the order so that "Promoter" appears at the top when flipped, modify the levels of the factor like this



RLFS_norm_nojuntn<- ggplot(RLFSs_rdna_summary, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized RLFS distribution in the monkey rDNA locus", 
       x= "monkey rDNA region", 
       y= "Normalized RLFS count", 
       fill = "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.60))+
  geom_text(aes(label= RLFS_count, hjust= -1.0, vjust= 0.5, size= 50))+
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

ggsave("Normalized_RLFS_distribution_in_monkey_rDNA_subcomponents_after_rule.tiff", 
       plot = RLFS_norm_nojuntn, width=18,height=10, dpi=150)



#to make template and non-template
entire_RLFSs_rdna_summary2<- entire_RLFSs_rdna %>% group_by(rDNA_region, strand) %>% count()
names(entire_RLFSs_rdna_summary2)[3] <- "RLFS_count"

new_rows<- data.table(rDNA_region = c("18S","5.8S","5.8S"),
                      strand= c("-","-","+"),
                      RLFS_count = c(0, 0,0))

entire_RLFSs_rdna_summary2<- rbind(entire_RLFSs_rdna_summary2, new_rows)
entire_RLFSs_rdna_summary2$rDNA_region <- factor(entire_RLFSs_rdna_summary2$rDNA_region, 
                                                 levels = c("5'ETS", "18S", "ITS1", "5.8S", 
                                                            "ITS2","28S", "3'ETS" ))

entire_RLFSs_rdna_summary2<- entire_RLFSs_rdna_summary2 %>% mutate(norm_RLFS_count = RLFS_count/sum(entire_RLFSs_rdna_summary2$RLFS_count)) %>% 
  mutate(norm_RLFS_count= round(norm_RLFS_count, 2))

fwrite(entire_RLFSs_rdna_summary2, "RLFS_KY962518_monkey_no_junctn_strandwise_AR_graphinput.csv")


rlfs_strandwise<- ggplot(entire_RLFSs_rdna_summary2, aes(x= rDNA_region, y = norm_RLFS_count, fill= strand)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized RLFS strandwise distribution in the monkey rDNA locus", 
       x= "monkey rDNA region", 
       y= "Normalized RLFS count", 
       fill= "RLFS strand")+
  scale_y_continuous(breaks= seq(0, 0.40, by = 0.1), limits =c(0,0.40))+
  geom_text(aes(label= RLFS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
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



ggsave( "Normalized_strandwise_RLFS_distribution_in_monkey_rDNA_subcomponents_AR.tiff", 
        plot = rlfs_strandwise, width=18,height=10, dpi=150)




nontemplate<- entire_RLFSs_rdna_summary2 %>% filter(strand == "+")
nontemplate$rDNA_region <- factor(nontemplate$rDNA_region, 
                                  levels = c("5'ETS", "18S", "ITS1", "5.8S", 
                                             "ITS2","28S", "3'ETS"))


template<- entire_RLFSs_rdna_summary2 %>% filter(strand == "-")
template$rDNA_region <- factor(template$rDNA_region, 
                               levels = c("5'ETS", "18S", "ITS1", "5.8S", 
                                          "ITS2","28S", "3'ETS"))



rlfs_nontemplate <- ggplot(nontemplate, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized Non-template RLFS distribution in the monkey rDNA locus", 
       x= "monkey rDNA region", 
       y= "Normalized Non-template RLFS count", 
       fill= "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.40, by = 0.1), limits =c(0,0.40))+
  geom_text(aes(label= RLFS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  
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


ggsave( "Normalized_nontemplate_RLFS_distribution_in_monkey_rDNA_subcomponents_AR.tiff", 
        plot = rlfs_nontemplate, width=18,height=10, dpi=150)

rlfs_template <- ggplot(template, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized Template RLFS distribution in the monkey rDNA locus", 
       x= "monkey rDNA region", 
       y= "Normalized Template RLFS count", 
       fill= "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.40, by = 0.1), limits =c(0,0.40))+
  geom_text(aes(label= RLFS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  
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


ggsave( "Normalized_template_RLFS_distribution_in_monkey_rDNA_subcomponents_AR.tiff", 
        plot = rlfs_template, width=18,height=10, dpi=150)




