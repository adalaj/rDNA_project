# ------------------------------------------------------------------------------
# This code is part of paper: In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus.
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
# This R script processes G-quadruplex forming sequence (G4FS) predictions in the 
# Human rDNA locus to assign each predicted G4FS to defined 
# rDNA subregions (Promoter, 5′ETS, 18S, ITS1, 5.8S, ITS2, 28S, 3′ETS, and IGS).
#

# Specifically, it:
#   1.Reads G4FS predictions generated using the G4 canonical finder 
#      (run in Python: g4_canonical_finder_3.11python.py).
#   2. Separates template and non-template strand G4FS predictions.
#   4. Assigns each G4FS to rDNA subregions (Promoter, 5′ETS, 18S, ITS1, 5.8S, ITS2, 28S, 3′ETS)
#      based on strand-specific start coordinates.
#   5. Counts and normalizes G4FS occurrences across subregions, including junctions,
#      following the rule that a G4FS is assigned to the region where it initiates.
#   6. Generates summary CSVs and high-resolution bar plots showing:
#        - Overall normalized G4FS distribution
#        - Strandwise (template vs non-template) distributions
#
# Inputs:
#   - Human rDNA FASTA sequence (GenBank: KY962518)
#   - G4FS output text file from g4_canonical_finder_3.11python.py
#
# Outputs:
#   - Annotated CSV file and plots of G4FS counts, density and proportion per rDNA region.
#
# ------------------------------------------------------------------------------





library(tidyverse)
library(data.table)


##i want plot bar graph after the rule for G4s that has 3500 bp added to 5ETS and also has promoter at the end. 
##basically input file will be output_G4FS_KY962518_added_3500nt_IGS_upstream_humanrDNA.txt after passing through G4FS algorithm.

##the rule: Counting the presence of G4s where it is first detected. For example, if G4s start towards the
# end of 5'ETS but stretches to 18S then it will counted under 5'ETS. 
# if we would have defined boundary in the first place then we would have lost this junction in our counting. 


entire_g4s_rdna<- fread("output_G4FS_KY962518_added_3500nt_IGS_upstream_humanrDNA.txt", sep = "\t", header = FALSE) #222
#222 is total number of G4FS that are found when both  the promoter (2202bp) and some  part of IGS (1297) are present twice
# this means after doing junction analysis we should have less than 222 
# to avoid double count i will filter the output_G4FS_KY962518_added_3500nt_IGS_upstream_humanrDNA.txt right after the the actual IGS length ends that is 46136.

# as from the bar graph visualization it is clear that we have 12 extra G4FS from promoter. 
colnames(entire_g4s_rdna) <- c("GenBank_Accession", "G4FS_start", "G4FS_end", "G4FS_sequence", "G4FS_name", "strand")

## strand specificty doesnt matter here because we are cropping region of our interest
entire_g4s_rdna<- entire_g4s_rdna[entire_g4s_rdna$G4FS_start >1299 & entire_g4s_rdna$G4FS_start <46137]#210 .. total length is 48338
entire_g4s_rdna$rDNA_region <- "junction"

entire_g4s_rdna$pGCS_length<- abs(entire_g4s_rdna$G4FS_end-entire_g4s_rdna$G4FS_start)


#strand specificity will matter here now we want to allocate G4s to rdna region sub components  based on their direction
# to do that, i will be creating a column that will have actual G4s start based on strand specificity

entire_g4s_rdna<- entire_g4s_rdna %>% mutate(actual_G4FS_start = ifelse(entire_g4s_rdna$strand == "+", G4FS_start, G4FS_end))
entire_g4s_rdna<- entire_g4s_rdna %>% mutate(actual_G4FS_end = ifelse(entire_g4s_rdna$strand=="+", G4FS_end, G4FS_start))


entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 1299 & entire_g4s_rdna$actual_G4FS_start < 3500] <- "Promoter"
sum(entire_g4s_rdna$rDNA_region=="Promoter")
#12

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 3500] <- "Promoter and 5'ETS junction"


entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start >= 3501] <- "5'ETS"

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 7157 ] <- "5'ETS and 18S junction"

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start >= 7158 ] <- "18S"
entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 9026 ] <- "18S and ITS1 junction"


entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start >= 9027  ] <- "ITS1"
entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 10096] <- "ITS1 and 5.8S junction"

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start >= 10097 ] <- "5.8S"
entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 10253 ] <- "5.8S and ITS2 junction"

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start >= 10254] <- "ITS2"
entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 11420 ] <- "ITS2 and 28S junction"

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start >= 11421 ] <- "28S"
entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 16471 ] <- "28S and 3'ETS junction"

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start >= 16472 ] <- "3'ETS"
entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 16832 ] <- "3'ETS and IGS junction"

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_G4FS_start > 16833 & entire_g4s_rdna$actual_G4FS_start < 46137] <- "IGS" 

#rDNA_region_length is needed for downstream normalization. 
#we will doing nromalization over rDNA length 


entire_g4s_rdna <- entire_g4s_rdna %>%
  mutate(rDNA_region_length = case_when(
    rDNA_region == "Promoter" ~ 2202,
    rDNA_region == "5'ETS"    ~ 3657,
    rDNA_region == "18S"      ~ 1869,
    rDNA_region == "ITS1"     ~ 1070,
    rDNA_region == "5.8S"     ~ 157,
    rDNA_region =="ITS2"      ~ 1167,
    rDNA_region =="28S"       ~ 5051,
    rDNA_region == "3'ETS"    ~ 361,
    rDNA_region == "IGS"      ~ 29305, 
    TRUE ~0 )) #unmatched make it zero


fwrite(entire_g4s_rdna, "G4FS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", sep = ",")

entire_g4s_rdna_summary<- entire_g4s_rdna %>% group_by(rDNA_region,rDNA_region_length) %>% count()
#rows order is lost 
sum(entire_g4s_rdna_summary$n)
#210, here we see that 12 extra G4FS from promoter are removed.

entire_g4s_rdna_summary[,1]

#1 28S         3'ETS      5'ETS      
#4 IGS          ITS1       ITS2       
#7 Promoter



new_rows<- data.table(rDNA_region = c("Promoter and 5'ETS junction", "5'ETS and 18S junction", "18S", 
                                      "18S and ITS1 junction","ITS1 and 5.8S junction",
                                      "5.8S", "5.8S and ITS2 junction","ITS2 and 28S junction",
                                      "28S and 3'ETS junction", "3'ETS and IGS junction"),
                      n=0,
                      rDNA_region_length= 0)




entire_g4s_rdna_summary<- rbind(entire_g4s_rdna_summary, new_rows)
entire_g4s_rdna_summary$rDNA_region
#[1] "28S"                         "3'ETS"                      
#[3] "5'ETS"                       "IGS"                        
#[5] "ITS1"                        "ITS2"                       
#[7] "Promoter"                    "Promoter and 5'ETS junction"
#[9] "5'ETS and 18S junction"      "18S"                        
#[11] "18S and ITS1 junction"       "ITS1 and 5.8S junction"     
#[13] "5.8S"                        "5.8S and ITS2 junction"     
#[15] "ITS2 and 28S junction"       "28S and 3'ETS junction"     
#[17] "3'ETS and IGS junction"     


row.names(entire_g4s_rdna_summary)<- entire_g4s_rdna_summary$rDNA_region


entire_g4s_rdna_summary <- entire_g4s_rdna_summary[c("Promoter","Promoter and 5'ETS junction", "5'ETS", "5'ETS and 18S junction", 
                                                     "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                     "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                     "28S and 3'ETS junction", "3'ETS", "3'ETS and IGS junction", "IGS" ),]

names(entire_g4s_rdna_summary)[2] <- "rDNA_region_length"
names(entire_g4s_rdna_summary)[3] <- "G4FS_count"
entire_g4s_rdna_summary <- entire_g4s_rdna_summary %>%
  mutate(rDNA_region_length = case_when(
    rDNA_region == "Promoter" ~ 2202,
    rDNA_region == "5'ETS"    ~ 3657,
    rDNA_region == "18S"      ~ 1869,
    rDNA_region == "ITS1"     ~ 1070,
    rDNA_region == "5.8S"     ~ 157,
    rDNA_region =="ITS2"      ~ 1167,
    rDNA_region =="28S"       ~ 5051,
    rDNA_region == "3'ETS"    ~ 361,
    rDNA_region == "IGS"      ~ 29305, 
    TRUE ~ 0)) #unmatched make it zero


#Which regions are enriched for G4FS relative to their size.
#Removes the bias of region length — highlights hotspots where G4FS are concentrated

entire_g4s_rdna_summary<- entire_g4s_rdna_summary %>% mutate(G4FS_density = G4FS_count/rDNA_region_length) %>% 
  mutate(G4FS_density= round(G4FS_density, 3))

fwrite(entire_g4s_rdna_summary, "G4FS_KY962518_added_3500nt_IGS_upstream_at_junctn_after_rule_graphinput.csv", sep = ",")


entire_g4s_rdna_summary<- fread("G4FS_KY962518_added_3500nt_IGS_upstream_at_junctn_after_rule_graphinput.csv", sep = ",", header = TRUE)
#Out of all G4FS across the rDNA locus, what fraction comes from each region?
entire_g4s_rdna_summary<- entire_g4s_rdna_summary %>% mutate(G4FS_proportion_perc = round((G4FS_count/sum(G4FS_count)*100),2))
#Which region contributes the largest share of G4FS overall.
#Biased toward longer regions (they naturally accumulate more G4FS simply because they have more bases)
fwrite(entire_g4s_rdna_summary, "G4FS_KY962518_added_3500nt_IGS_upstream_at_junctn_after_rule_graphinput.csv", sep = ",")



#Fig 3C to 3E
entire_g4s_rdna_summary<- fread("G4FS_KY962518_added_3500nt_IGS_upstream_at_junctn_after_rule_graphinput.csv", sep = ",", header = TRUE)
g4s_rdna_summary<- entire_g4s_rdna_summary[!grepl("junction", entire_g4s_rdna_summary$rDNA_region),]


#When you use coord_flip(), the order of the factor levels in rDNA_region determines how the categories are displayed along the flipped y-axis.
#By default, the first factor level appears at the bottom when flipped, and the last factor level appears at the top.

g4s_rdna_summary$rDNA_region <- factor(g4s_rdna_summary$rDNA_region, 
                                         levels = rev(c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                        "ITS2","28S", "3'ETS", "IGS" )))



#To reverse the order so that "Promoter" appears at the top when flipped, modify the levels of the factor like this

max_value<- max(entire_g4s_rdna_summary$G4FS_density, na.rm = TRUE) #i added 0.001 because round was making 0.015 to  0.01 instead of 0.02. 

G4FS_norm_3500igs_nojuntn<- ggplot(g4s_rdna_summary, aes(x= rDNA_region, y = G4FS_density, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(#title= "Normalized G4FS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "G4FS density", 
       fill = "rDNA")+
  scale_y_continuous(breaks= seq(0, max_value, by = 0.005), limits =c(0,max_value), expand = expansion(mult = c(0, 0.08)))+
  #geom_text(aes(label= G4FS_count, hjust= -0.2, vjust= 0.5), size= 30)+
  scale_fill_manual(values= rev(c("#B6FFF4", "#FDCCE5","#D0B6FF", "#EF9B20", "#A0322B", 
                                  "#FFCC17", "#E5FFB6", "#3B8CC4", "#A4A2A8")))+
  #guides(fill = guide_legend(reverse = TRUE))
  theme_minimal()+
  theme(axis.ticks = element_line(color = "black", linewidth = 4),
        axis.ticks.length = unit(50, "pt"),
       panel.grid = element_blank(),
       plot.title = element_text(hjust = 0.5), #face = "bold"),
       plot.subtitle = element_text(hjust = 0.5),
       text = element_text(size = 100),
       axis.line = element_line(color = "black", linewidth = 4),
       axis.title.x = element_text(vjust = 0.7, hjust = 0.5, color = "black"),
       axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, margin = margin (r=20)),  # Center Y-axis title
       axis.text.x  = element_text(color = "black"),
       axis.text.y  = element_text(color = "black"), 
       legend.position = "none")+
  coord_flip()

ggsave( "Normalized_G4FS_distribution_in_human_rDNA_subcomponents_after_rule.png", 
        plot = G4FS_norm_3500igs_nojuntn, width=30,height=18, dpi=300)



max_value<- max(g4s_rdna_summary$G4FS_proportion_perc, na.rm = TRUE)

G4FS_prop_3500igs_nojuntn<- ggplot(g4s_rdna_summary, aes(x= rDNA_region, y = G4FS_proportion_perc, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(#title= "Normalized G4FS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "G4FS proportion (%)", 
       fill = "rDNA")+
  scale_y_continuous(breaks= seq(0, max_value, by = 10), limits =c(0,max_value), expand = expansion(mult = c(0, 0.08)))+
  #geom_text(aes(label= G4FS_count, hjust= -0.2, vjust= 0.5), size= 30)+
  scale_fill_manual(values= rev(c("#B6FFF4", "#FDCCE5","#D0B6FF", "#EF9B20", "#A0322B", 
                                  "#FFCC17", "#E5FFB6", "#3B8CC4", "#A4A2A8")))+
  #guides(fill = guide_legend(reverse = TRUE))
  theme_minimal()+
  theme(axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.ticks.length = unit(50, "pt"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 100),
        axis.line = element_line(color = "black", linewidth = 4),
        axis.title.x = element_text(vjust = 0.7, hjust = 0.5, colour = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, margin = margin(r=20)),   # Center Y-axis title
        axis.text.x  = element_text(color = "black"),
        axis.text.y  = element_text(color = "black"),
        legend.position = "none")+
  coord_flip()

ggsave("G4FS_proportion_distribution_in_human_rDNA_subcomponents_after_rule.png", 
       plot = G4FS_prop_3500igs_nojuntn, width=30,height=18, dpi=600)



entire_g4s_rdna<- fread("G4FS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", sep = ",", header = TRUE)

#to make template and non-template
entire_g4s_rdna_summary2<- entire_g4s_rdna %>% group_by(rDNA_region, rDNA_region_length, strand) %>% count()
names(entire_g4s_rdna_summary2)[4] <- "G4FS_count"

new_rows<- data.table(rDNA_region = c( "Promoter", "ITS1", "18S","18S", "5.8S","5.8S"),
                      rDNA_region_length = c(0,0,0,0,0,0),
                      strand= c("+", "-", "+", "-", "+", "-"),
                      G4FS_count = c(0, 0, 0,0,0,0))

entire_g4s_rdna_summary2<- rbind(entire_g4s_rdna_summary2, new_rows)

entire_g4s_rdna_summary2 <- entire_g4s_rdna_summary2 %>%
  mutate(rDNA_region_length = case_when(
    rDNA_region == "Promoter" ~ 2202,
    rDNA_region == "5'ETS"    ~ 3657,
    rDNA_region == "18S"      ~ 1869,
    rDNA_region == "ITS1"     ~ 1070,
    rDNA_region == "5.8S"     ~ 157,
    rDNA_region =="ITS2"      ~ 1167,
    rDNA_region =="28S"       ~ 5051,
    rDNA_region == "3'ETS"    ~ 361,
    rDNA_region == "IGS"      ~ 29305, 
    TRUE ~ 0)) #unmatched make it zero



entire_g4s_rdna_summary2$rDNA_region <- factor(entire_g4s_rdna_summary2$rDNA_region, 
                                                 levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                            "ITS2","28S", "3'ETS", "IGS" ))

entire_g4s_rdna_summary2<- entire_g4s_rdna_summary2 %>% mutate(G4FS_density = G4FS_count/rDNA_region_length) %>% 
  mutate(G4FS_density= round(G4FS_density, 4))

fwrite(entire_g4s_rdna_summary2, "G4FS_KY962518_added_3500nt_IGS_upstream_no_junctn_strandwise_graphinput.csv")

entire_g4s_rdna_summary2<- fread("G4FS_KY962518_added_3500nt_IGS_upstream_no_junctn_strandwise_graphinput.csv", sep = ",", header = TRUE)

entire_g4s_rdna_summary2$strand <- factor(
  entire_g4s_rdna_summary2$strand,
  levels = rev(c("+", "-"))  # "+" = Non-template first, "-" = Template second
)

entire_g4s_rdna_summary2$rDNA_region <- factor(entire_g4s_rdna_summary2$rDNA_region, 
                                                 levels = rev(c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                                "ITS2","28S", "3'ETS", "IGS" )))

max_value<- max(entire_g4s_rdna_summary2$G4FS_density, na.rm = TRUE)#i added 0.001 because round was making 0.015 to  0.01 instead of 0.02. 

G4FS_strandwise_flip<- ggplot(entire_g4s_rdna_summary2, aes(x= rDNA_region, y = G4FS_density, fill= strand)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(#title= "Normalized G4FS strandwise distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "G4FS density", 
       fill= NULL)+
  scale_y_continuous(breaks= seq(0, max_value, by = 0.003), limits =c(0,max_value),expand = expansion(mult = c(0, 0.05)))+
  #geom_text(aes(label= G4FS_count, hjust=-0.2, vjust=0.5), size=20, position = position_dodge(width = 0.9))+
  scale_fill_manual(values= c("+" = "#E21515", "-" = "#1414E1"), 
                    labels = c("+" = "Non-template strand", "-" = "Template strand"),
                    breaks = c("+", "-"))+
  #scale_fill_manual(values = combined_colors)+
  theme_minimal()+
  theme(axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.ticks.length = unit(50, "pt"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5), #face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 100),
        axis.line = element_line(color = "black", linewidth = 4),
        axis.title.x = element_text(vjust = 0.6, hjust = 0.5, colour = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),   # Center Y-axis title
        axis.text.x  = element_text(color = "black"),
        axis.text.y  = element_text(color = "black"),
        legend.position = "top", 
        legend.text = element_text(size=80),
        legend.key.size = unit(3, "cm"))+
  coord_flip()


ggsave( "Normalized_strandwise_G4FS_flipped_distribution_in_human_rDNA_subcomponents_after_rule.png", 
        plot = G4FS_strandwise_flip, width=31,height=18, dpi=600)







