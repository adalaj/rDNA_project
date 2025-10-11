# ------------------------------------------------------------------------------
# This code is part of paper: In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus.
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
#   Quantifies and visualizes R-loop initiation zone (RIZ) across the human
#   rDNA locus (KY962518 + 3.5 kb upstream IGS) based on their origin.
#   RIZ are assigned to the region where they first appear (e.g., an RIZ
#   spanning 5′ETS–18S is counted under 5′ETS). Outputs regional counts,
#   normalized densities, and strand-specific RIZ distributions.
#
# Major Steps:
#   1. Import QMRLFS output BED file and filter to promoter–IGS range (1299–46136 bp).
#   2. Assign RIZ to rDNA subregions based on start coordinate and strand.
#   3. Compute region-specific counts, lengths, normalized densities, and proportions.
#   4. Generate bar plots:
#        - RIZ proportion per rDNA region
#        - Length-normalized RIZ density
#        - Strandwise RIZ density (template vs non-template)
#
# Input:
#   - Human rDNA FASTA sequence (GenBank: KY962518)
#   - RIZ output text file from QmRLFS_RIZ_finder.py
#
#
# Output:
#  - Annotated CSV file "RIZ_KY962518_added_3500nt_IGS_upstream_at_junctn_details_after_rule.csv"and 
#  plots of RIZ counts, density and proportion per rDNA region (not deposited in git).
# ------------------------------------------------------------------------------


library(tidyverse)
library(data.table)

{#(python2.7)  Downloads % python QmRLFS_RIZ_finder.py -bed -i test.fasta -o test_rlfs
#(python2.7)  %python QmRLFS-finder.py -bed -i test.fasta -o test_riz


rlfs<- fread("test_rlfs.out.bed", sep = "\t", header = FALSE) #50
rlfs_plus<- rlfs %>% filter(V6=="+") %>% 
  select(1,2)

rlfs_minus<- rlfs %>% filter(V6=="-") %>% 
  select(1,3) # this column because I am only interested in where RLFS or RIZ start, just to see if all my initiation part in RLFS code  
#captured by updated RIZ code. 


riz<- fread("test_riz.out.bed", sep = "\t", header = FALSE) #64

riz_plus<- riz %>% filter(V6=="+") %>% 
  select(1,2)

riz_minus<- riz %>% filter(V6=="-") %>% 
  select(1,3)



setdiff(rlfs_plus, riz_plus)
#All plus ones are captured
#Empty data.table (0 rows and 2 cols): V1,V2



setdiff(rlfs_minus, riz_minus)
#Empty data.table (0 rows and 2 cols): V1,V3

#this shows that all riz code captures more of RIZ than RLFS, simpley because alot of REZ can also form RIZ if they satisfy the RIZ condition.

}

#Now, as RIZ code is up and running i would like to run on KY962518_added_3500nt_IGS_upstream_nontemplate.fasta 

#(python2.7) Downloads % python QmRLFS_RIZ_finder.py -bed -i KY962518_added_3500nt_IGS_upstream_nontemplate.fasta -o KY962518_added_3500nt_IGS_upstream_RIZ
#QmRLFS_RIZ_finder.py (version v1.5)
#run on Thu Aug 14 2025 18:06:22 
#command line: python QmRLFS_RIZ_finder.py -bed -i KY962518_added_3500nt_IGS_upstream_nontemplate.fasta -o KY962518_added_3500nt_IGS_upstream_RIZ
#Time used: 0.00 mins

#(python2.7) Downloads % python QmRLFS_RIZ_finder.py -i KY962518_added_3500nt_IGS_upstream_nontemplate.fasta -o KY962518_added_3500nt_IGS_upstream_RIZ 
#QmRLFS_RIZ_finder.py (version v1.5)
#run on Thu Aug 14 2025 18:06:51 
#command line: python QmRLFS_RIZ_finder.py -i KY962518_added_3500nt_IGS_upstream_nontemplate.fasta -o KY962518_added_3500nt_IGS_upstream_RIZ
#Time used: 0.00 mins

entire_RIZs_rdna<- fread("KY962518_added_3500nt_IGS_upstream_RIZ.out.bed", sep = "\t", header = FALSE) #304
#208 was RLFS which is now changed to 304 just RIZ when both  the promoter (2202bp) and some  part of IGS (1297) are present twice

# this means after doing junction analysis we should have less than 304
# to avoid double count i will filter the KY962518_added_3500nt_IGS_upstream_RIZ.out.bed right after the the actual IGS length ends that is 46136.

# as from the bar graph visualization it is clear that we have 13 extra RIZ from promoter. 

entire_RIZs_rdna<- entire_RIZs_rdna %>% select(1:6)
colnames(entire_RIZs_rdna) <- c("GenBank_Accession", "RIZ_start", "RIZ_end", "RIZ_name","score", "strand")

# strand specificity doesnt matter here because we are cropping region of our interest
entire_RIZs_rdna<- entire_RIZs_rdna[entire_RIZs_rdna$RIZ_start >1299 & entire_RIZs_rdna$RIZ_start < 46137] #286 RIZ and 195 RLFS
#after 1298 is where the promoter begins and 46137 is where IGS end (48338- promoter length which is 2202=46136)

entire_RIZs_rdna$rDNA_region <- "junction"



#strand specificity will matter here now we want to allocate RIZ to rdna region sub components  based on their direction
# to do that, i will be creating a column that will have actual RIZ start based on strand specificity

entire_RIZs_rdna<- entire_RIZs_rdna %>% mutate(actual_RIZ_start = ifelse(entire_RIZs_rdna$strand == "+", RIZ_start, RIZ_end))
entire_RIZs_rdna<- entire_RIZs_rdna %>% mutate(actual_RIZ_end = ifelse(entire_RIZs_rdna$strand == "+", RIZ_end, RIZ_start))

entire_RIZs_rdna$RIZ_length<- abs(entire_RIZs_rdna$RIZ_start-entire_RIZs_rdna$RIZ_end)

entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start > 1299 & entire_RIZs_rdna$actual_RIZ_start < 3501] <- "Promoter"
sum(entire_RIZs_rdna$rDNA_region=="Promoter")#increased from 13 RLFS  to 18RIZ

entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start > 3500] <- "Promoter and 5'ETS junction"

entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start >= 3501] <- "5'ETS"
entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start > 7157] <- "5'ETS and 18S junction"

entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start >= 7158] <- "18S"
entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start > 9026] <- "18S and ITS1 junction"

entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start >= 9027] <- "ITS1"
entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start > 10096] <- "ITS1 and 5.8S junction"

entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start >= 10097] <- "5.8S"
entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start > 10253] <- "5.8S and ITS2 junction"

entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start >= 10254] <- "ITS2"
entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start > 11420] <- "ITS2 and 28S junction"

entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start >= 11421] <- "28S"
entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start > 16471] <- "28S and 3'ETS junction"

entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start >= 16472] <- "3'ETS"
entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start > 16832] <- "3'ETS and IGS junction"

entire_RIZs_rdna$rDNA_region[entire_RIZs_rdna$actual_RIZ_start > 16833 & entire_RIZs_rdna$actual_RIZ_start < 46137] <- "IGS"


#rDNA_region_length is needed for downstream normalization. 
#we will doing normalization over rDNA length 


entire_RIZs_rdna <- entire_RIZs_rdna %>%
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



fwrite(entire_RIZs_rdna, "RIZ_KY962518_added_3500nt_IGS_upstream_at_junctn_details_after_rule.csv", sep = ",")


#Extra and not needed for manuscript
{entire_RIZs_rdna_summary<- entire_RIZs_rdna %>% group_by(rDNA_region, rDNA_region_length) %>% count()
#rows order is lost 

sum(entire_RIZs_rdna_summary$n)
#[1] 286 

entire_RIZs_rdna_summary[,1]

#1 18S      28S         3'ETS            
#4 5'ETS    IGS          ITS1       
#7 ITS2     Promoter


new_rows<- data.table(rDNA_region = c("Promoter and 5'ETS junction", "5'ETS and 18S junction", 
                                      "18S and ITS1 junction","ITS1 and 5.8S junction",
                                      "5.8S", "5.8S and ITS2 junction","ITS2 and 28S junction",
                                      "28S and 3'ETS junction", "3'ETS and IGS junction"), 
                      n= 0,
                      rDNA_region_length=0)

entire_RIZs_rdna_summary<- rbind(entire_RIZs_rdna_summary, new_rows)
entire_RIZs_rdna_summary$rDNA_region
#[1] "18S"                         "28S"                        
#[3] "3'ETS"                       "5'ETS"                      
#[5] "IGS"                         "ITS1"                       
#[7] "ITS2"                        "Promoter"                   
#[9] "Promoter and 5'ETS junction" "5'ETS and 18S junction"     
#[11] "18S and ITS1 junction"       "ITS1 and 5.8S junction"     
#[13] "5.8S"                        "5.8S and ITS2 junction"     
#[15] "ITS2 and 28S junction"       "28S and 3'ETS junction"     
#[17] "3'ETS and IGS junction" 


row.names(entire_RIZs_rdna_summary)<- entire_RIZs_rdna_summary$rDNA_region

entire_RIZs_rdna_summary <- entire_RIZs_rdna_summary[c("Promoter","Promoter and 5'ETS junction", "5'ETS", "5'ETS and 18S junction", 
                                                       "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                       "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                       "28S and 3'ETS junction", "3'ETS", "3'ETS and IGS junction", "IGS" ),]
names(entire_RIZs_rdna_summary)[2] <- "rDNA_region_length"
names(entire_RIZs_rdna_summary)[3] <- "RIZ_count"
entire_RIZs_rdna_summary <- entire_RIZs_rdna_summary %>%
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


entire_RIZs_rdna_summary<- entire_RIZs_rdna_summary %>% mutate(RIZ_density = RIZ_count/rDNA_region_length)%>% 
  mutate(RIZ_density= round(RIZ_density, 3))

fwrite(entire_RIZs_rdna_summary, "RIZ_KY962518_added_3500nt_IGS_upstream_at_junctn_after_rule_graphinput.csv", sep = ",")

entire_RIZs_rdna_summary$rDNA_region <- factor(entire_RIZs_rdna_summary$rDNA_region, 
                                               levels = c("Promoter","Promoter and 5'ETS junction", "5'ETS", "5'ETS and 18S junction", 
                                                          "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                          "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                          "28S and 3'ETS junction", "3'ETS", "3'ETS and IGS junction", "IGS" ))


max_value<- round(max(entire_RIZs_rdna_summary$RIZ_density, na.rm = TRUE),2)

RIZ_norm_3500igs<- ggplot(entire_RIZs_rdna_summary, aes(x= rDNA_region, y = RIZ_density, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized RIZ distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "RIZ density", 
       fill = "rDNA region")+
  scale_y_continuous(breaks= seq(0, max_value, by = 0.01), limits =c(0,max_value))+
  geom_text(aes(label= RIZ_count), vjust= -0.5, size= 5)+
  scale_fill_manual(values= c( "#FFB6C1","maroon", "#D0B6FF", "steelblue", "#E5FFB6","darkviolet", "#FFE0C2","burlywood2", "#B6FFF4", 
                               "pink4", "#FFFFE0","aquamarine", "#E8E8FB","greenyellow", "#B6E5FF","turquoise2", "#DCDCDC"))+
  #scale_fill_manual(values = combined_colors)+
  theme_minimal()+
  theme(axis.text.x = element_blank(), 
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),   # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))
#coord_flip()

ggsave( "Normalized_RIZ_distribution_in_human_rDNA_subcomponents_incld_junctn_AR.tiff", 
        plot = RIZ_norm_3500igs, width=18,height=10, dpi=150) #AR is after rule, bcoz powerpoint is not accepting too long image name



RIZs_rdna_summary<- entire_RIZs_rdna_summary[!grepl("junction", entire_RIZs_rdna_summary$rDNA_region),]


#When you use coord_flip(), the order of the factor levels in rDNA_region determines how the categories are displayed along the flipped y-axis.
#By default, the first factor level appears at the bottom when flipped, and the last factor level appears at the top.

RIZs_rdna_summary$rDNA_region <- factor(RIZs_rdna_summary$rDNA_region, 
                                        levels = rev(c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                       "ITS2","28S", "3'ETS", "IGS" )))


#To reverse the order so that "Promoter" appears at the top when flipped, modify the levels of the factor like this

max_value<- round(max(RIZs_rdna_summary$RIZ_density, na.rm = TRUE),2)


RIZ_norm_3500igs_nojuntn<- ggplot(RIZs_rdna_summary, aes(x= rDNA_region, y = RIZ_density, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized RIZ distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "RIZ density", 
       fill = "rDNA")+
  scale_y_continuous(breaks= seq(0, max_value, by = 0.01), limits =c(0,max_value))+
  geom_text(aes(label= RIZ_count, hjust= -1.0, vjust= 0.5, size= 50))+
  scale_fill_manual(values= rev(c("#B6FFF4", "#FDCCE5","#D0B6FF", "#EF9B20", "#A0322B", 
                                  "#FFCC17", "#E5FFB6", "#3B8CC4", "#A4A2A8")))+
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

ggsave("Normalized_RIZ_distribution_in_human_rDNA_subcomponents_after_rule.tiff", 
       plot = RIZ_norm_3500igs_nojuntn, width=18,height=10, dpi=150)



#to make template and non-template

entire_RIZs_rdna<- fread("RIZ_KY962518_added_3500nt_IGS_upstream_at_junctn_details_after_rule.csv", sep = ",", header = TRUE)


entire_RIZs_rdna_summary2<- entire_RIZs_rdna %>% group_by(rDNA_region, rDNA_region_length, strand) %>% count()
names(entire_RIZs_rdna_summary2)[4] <- "RIZ_count"

new_rows<- data.table(rDNA_region = c("5.8S","5.8S"),
                      rDNA_region_length = c(0, 0),
                      strand= c( "+", "-"),
                      RIZ_count = c(0,0))

entire_RIZs_rdna_summary2<- rbind(entire_RIZs_rdna_summary2, new_rows)


entire_RIZs_rdna_summary2 <- entire_RIZs_rdna_summary2 %>%
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


entire_RIZs_rdna_summary2$rDNA_region <- factor(entire_RIZs_rdna_summary2$rDNA_region, 
                                                levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                           "ITS2","28S", "3'ETS", "IGS" ))

entire_RIZs_rdna_summary2<- entire_RIZs_rdna_summary2 %>% mutate(RIZ_density = RIZ_count/rDNA_region_length) %>% 
  mutate(RIZ_density= round(RIZ_density, 4))

fwrite(entire_RIZs_rdna_summary2, "RIZ_KY962518_added_3500nt_IGS_upstream_no_junctn_strandwise_AR_graphinput.csv")


max_value<- round(max(entire_RIZs_rdna_summary2$RIZ_density, na.rm = TRUE),4)

RIZ_strandwise<- ggplot(entire_RIZs_rdna_summary2, aes(x= rDNA_region, y = RIZ_density, fill= strand)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized RIZ strandwise distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "RIZ density", 
       fill= "RIZ strand")+
  scale_y_continuous(breaks= seq(0, max_value, by = 0.002), limits =c(0,max_value))+
  geom_text(aes(label= RIZ_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  scale_fill_manual(values= c("+" = "#E21515", "-" = "#1414E1"), #changed the non template and template colors
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
        axis.ticks.y = element_line(color = "black"))#facet_wrap (~strand or ~rDNA region doesnt look good)



ggsave( "Normalized_strandwise_RIZ_distribution_in_human_rDNA_subcomponents_AR.tiff", 
        plot = RIZ_strandwise, width=18,height=10, dpi=150)


entire_RIZs_rdna_summary2$rDNA_region <- factor(entire_RIZs_rdna_summary2$rDNA_region, 
                                                levels = rev(c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                               "ITS2","28S", "3'ETS", "IGS" )))

max_value<- round(max(entire_RIZs_rdna_summary2$RIZ_density, na.rm = TRUE),4)

RIZ_strandwise_flip<- ggplot(entire_RIZs_rdna_summary2, aes(x= rDNA_region, y = RIZ_density, fill= strand)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized RIZ strandwise distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "RIZ density", 
       fill= "RIZ strand")+
  scale_y_continuous(breaks= seq(0, max_value, by = 0.004), limits =c(0,max_value))+
  geom_text(aes(label= RIZ_count, hjust= -1.0, vjust= 0.5, size= 50), position = position_dodge(width = 0.9))+
  scale_fill_manual(values= c("+" = "#E21515", "-" = "#1414E1"), #changed the non template and template colors
                    labels = c("+" = "Non-template", "-" = "Template"))+
  #scale_fill_manual(values = combined_colors)+
  theme_minimal()+
  theme(axis.title.x = element_text(vjust = -0.5, hjust = 0.5),
        axis.ticks.x = element_line(color = "black"), 
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),   # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))+
  coord_flip()


ggsave( "Normalized_strandwise_RIZ_flipped_distribution_in_human_rDNA_subcomponents_AR.tiff", 
        plot = RIZ_strandwise_flip, width=18,height=10, dpi=150)



nontemplate<- entire_RIZs_rdna_summary2 %>% filter(strand == "+")
nontemplate$rDNA_region <- factor(nontemplate$rDNA_region, 
                                  levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                             "ITS2","28S", "3'ETS", "IGS" ))
max_value<- round(max(nontemplate$RIZ_density),4)




RIZ_nontemplate <- ggplot(nontemplate, aes(x= rDNA_region, y = RIZ_density, fill= rDNA_region)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Non-template RIZ density distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Non-template RIZ density", 
       fill= "rDNA")+
  scale_y_continuous(breaks= seq(0, max_value, by = 0.002), limits =c(0,max_value))+
  geom_text(aes(label= RIZ_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  
  scale_fill_manual(values= c("#FFB6C1", "#D0B6FF", "#E5FFB6", "#FFE0C2", "#B6FFF4", 
                              "#FFFFE0", "#E8E8FB", "#B6E5FF", "#DCDCDC"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))


ggsave( "Normalized_nontemplate_RIZ_distribution_in_human_rDNA_subcomponents_AR.tiff", 
        plot = RIZ_nontemplate, width=18,height=10, dpi=150)



template<- entire_RIZs_rdna_summary2 %>% filter(strand == "-")
template$rDNA_region <- factor(template$rDNA_region, 
                               levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                          "ITS2","28S", "3'ETS", "IGS" ))

max_value<- round(max(template$RIZ_density),4)


RIZ_template <- ggplot(template, aes(x= rDNA_region, y = RIZ_density, fill= rDNA_region)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Template RIZ density distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Template RIZ density", 
       fill= "rDNA")+
  scale_y_continuous(breaks= seq(0, max_value, by = 0.002), limits =c(0,max_value))+
  geom_text(aes(label= RIZ_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  
  scale_fill_manual(values= c("#FFB6C1", "#D0B6FF", "#E5FFB6", "#FFE0C2", "#B6FFF4", 
                              "#FFFFE0", "#E8E8FB", "#B6E5FF", "#DCDCDC"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))


ggsave( "Normalized_template_RIZ_distribution_in_human_rDNA_subcomponents_AR.tiff", 
        plot = RIZ_template, width=18,height=10, dpi=150)

}
