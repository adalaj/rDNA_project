
#To make plot for G4S distribution in human rDNA locus 


library(tidyverse)
library(data.table)


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/g4s_and_rdna/human")



##the best approach that got approved by Bruce!
##i want plot bar graph after the rule for G4s that has 3500 bp added to 5ETS and also has promoter at the end. 
##basically input file will be output_pG4CS_KY962518_added_3500nt_IGS_upstream_humanrDNA.txt after passing through pG4CS algorithm.

##the rule: Counting the presence of G4s where it is first detected. For example, if G4s start towards the
# end of 5'ETS but stretches to 18S then it will counted under 5'ETS. 
# if we would have defined boundary in the first place then we would have lost this junction in our counting. 


entire_g4s_rdna<- fread("output_pG4CS_KY962518_added_3500nt_IGS_upstream_humanrDNA.txt", sep = "\t", header = FALSE) #222
#222 is total number of pG4CS that are found when both  the promoter (2202bp) and some  part of IGS (1297) are present twice
# this means after doing junction analysis we should have less than 222 
# to avoid double count i will filter the output_pG4CS_KY962518_added_3500nt_IGS_upstream_humanrDNA.txt right after the the actual IGS length ends that is 46136.

# as from the bar graph visualization it is clear that we have 12 extra pG4CS from promoter. 
colnames(entire_g4s_rdna) <- c("GenBank_Accession", "pG4CS_start", "pG4CS_end", "pG4CS_sequence", "pG4CS_name", "strand")

## strand specificty doesnt matter here because we are cropping region of our interest
entire_g4s_rdna<- entire_g4s_rdna[entire_g4s_rdna$pG4CS_start >1299 & entire_g4s_rdna$pG4CS_start <46137]#210 .. total length is 48338
entire_g4s_rdna$rDNA_region <- "junction"

entire_g4s_rdna$pGCS_length<- abs(entire_g4s_rdna$pG4CS_end-entire_g4s_rdna$pG4CS_start)


#strand specificity will matter here now we want to allocate G4s to rdna region sub components  based on their direction
# to do that, i will be creating a column that will have actual G4s start based on strand specificity

entire_g4s_rdna<- entire_g4s_rdna %>% mutate(actual_pG4CS_start = ifelse(entire_g4s_rdna$strand == "+", pG4CS_start, pG4CS_end))


entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start > 1299 & entire_g4s_rdna$actual_pG4CS_start < 3500] <- "Promoter"
sum(entire_g4s_rdna$rDNA_region=="Promoter")
#12

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start > 3500] <- "Promoter and 5'ETS junction"


entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start >= 3501] <- "5'ETS"

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start > 7157 ] <- "5'ETS and 18S junction"

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start >= 7158 ] <- "18S"
entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start > 9026 ] <- "18S and ITS1 junction"


entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start >= 9027  ] <- "ITS1"
entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start > 10096] <- "ITS1 and 5.8S junction"

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start >= 10097 ] <- "5.8S"
entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start > 10253 ] <- "5.8S and ITS2 junction"

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start >= 10254] <- "ITS2"
entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start > 11420 ] <- "ITS2 and 28S junction"

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start >= 11421 ] <- "28S"
entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start > 16471 ] <- "28S and 3'ETS junction"

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start >= 16472 ] <- "3'ETS"
entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start > 16832 ] <- "3'ETS and IGS junction"

entire_g4s_rdna$rDNA_region[entire_g4s_rdna$actual_pG4CS_start > 16833 & entire_g4s_rdna$actual_pG4CS_start < 46137] <- "IGS" 



fwrite(entire_g4s_rdna, "pG4CS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", sep = ",")

entire_g4s_rdna_summary<- entire_g4s_rdna %>% group_by(rDNA_region) %>% count()
#rows order is lost 
sum(entire_g4s_rdna_summary$n)
#210, here we see that 12 extra pG4CS from promoter are removed.

entire_g4s_rdna_summary[,1]

#1 28S         3'ETS      5'ETS      
#4 IGS          ITS1       ITS2       
#7 Promoter



new_rows<- data.table(rDNA_region = c("Promoter and 5'ETS junction", "5'ETS and 18S junction", "18S", 
                                      "18S and ITS1 junction","ITS1 and 5.8S junction",
                                      "5.8S", "5.8S and ITS2 junction","ITS2 and 28S junction",
                                      "28S and 3'ETS junction", "3'ETS and IGS junction"))

new_rows$n<- 0


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


names(entire_g4s_rdna_summary)[2] <- "pG4CS_count"

entire_g4s_rdna_summary<- entire_g4s_rdna_summary %>% mutate(norm_pG4CS_count = pG4CS_count/sum(entire_g4s_rdna_summary$pG4CS_count)) %>% 
  mutate(norm_pG4CS_count= round(norm_pG4CS_count, 2))

fwrite(entire_g4s_rdna_summary, "pG4CS_KY962518_added_3500nt_IGS_upstream_at_junctn_graphinput.csv", sep = ",")

entire_g4s_rdna_summary$rDNA_region <- factor(entire_g4s_rdna_summary$rDNA_region, 
                                                levels = c("Promoter","Promoter and 5'ETS junction", "5'ETS", "5'ETS and 18S junction", 
                                                           "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                           "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                           "28S and 3'ETS junction", "3'ETS", "3'ETS and IGS junction", "IGS" ))


pG4CS_norm_3500igs<- ggplot(entire_g4s_rdna_summary, aes(x= rDNA_region, y = norm_pG4CS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized pG4CS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalized pG4CS count", 
       fill = "rDNA region")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.60))+
  geom_text(aes(label= pG4CS_count), vjust= -0.5, size= 5)+
  scale_fill_manual(values= c( "#F5FEFB","maroon", "#E21515", "steelblue", "#5AAA46","darkviolet", "#F36017","burlywood2", "#6B1519", 
                               "pink", "#818689","aquamarine", "#ECE612","greenyellow", "#E07F80","turquoise2", "#DE9A22"))+
  #scale_fill_manual(values = combined_colors)+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5), +  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))
        #coord_flip()

ggsave( "Normalized_pG4CS_distribution_in_human_rDNA_subcomponents_incld_junctn_after_rule.tiff", 
        plot = pG4CS_norm_3500igs, width=18,height=10, dpi=150)


g4s_rdna_summary<- entire_g4s_rdna_summary[!grepl("junction", entire_g4s_rdna_summary$rDNA_region),]


#When you use coord_flip(), the order of the factor levels in rDNA_region determines how the categories are displayed along the flipped y-axis.
#By default, the first factor level appears at the bottom when flipped, and the last factor level appears at the top.

g4s_rdna_summary$rDNA_region <- factor(g4s_rdna_summary$rDNA_region, 
                                         levels = rev(c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                        "ITS2","28S", "3'ETS", "IGS" )))



#To reverse the order so that "Promoter" appears at the top when flipped, modify the levels of the factor like this



pG4CS_norm_3500igs_nojuntn<- ggplot(g4s_rdna_summary, aes(x= rDNA_region, y = norm_pG4CS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized pG4CS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalized pG4CS count", 
       fill = "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.60))+
  geom_text(aes(label= pG4CS_count, hjust= -1.0, vjust= 0.5, size= 50))+
  scale_fill_manual(values= rev(c("#F5FEFB", "#E21515", "#5AAA46", "#F36017", "#6B1519", 
                                  "#818689", "#ECE612", "#E07F80", "#DE9A22")))+
  #guides(fill = guide_legend(reverse = TRUE))
  theme_minimal()+
  theme(axis.text.x = element_blank(), 
       panel.grid = element_blank(),
       plot.title = element_text(hjust = 0.5, face = "bold"),
       plot.subtitle = element_text(hjust = 0.5),
       text = element_text(size = 30),
       axis.line = element_line(color = "black"),
       axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
       axis.ticks.y = element_line(color = "black"))+
  coord_flip()

ggsave( "Normalized_pG4CS_distribution_in_human_rDNA_subcomponents_after_rule.tiff", 
        plot = pG4CS_norm_3500igs_nojuntn, width=18,height=10, dpi=150)



#to make template and non-template
entire_g4s_rdna_summary2<- entire_g4s_rdna %>% group_by(rDNA_region, strand) %>% count()
names(entire_g4s_rdna_summary2)[3] <- "pG4CS_count"

new_rows<- data.table(rDNA_region = c( "Promoter", "ITS1", "18S","18S", "5.8S","5.8S"),
                      strand= c("+", "-", "+", "-", "+", "-"),
                      pG4CS_count = c(0, 0, 0,0,0,0))

entire_g4s_rdna_summary2<- rbind(entire_g4s_rdna_summary2, new_rows)
entire_g4s_rdna_summary2$rDNA_region <- factor(entire_g4s_rdna_summary2$rDNA_region, 
                                                 levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                            "ITS2","28S", "3'ETS", "IGS" ))

entire_g4s_rdna_summary2<- entire_g4s_rdna_summary2 %>% mutate(norm_pG4CS_count = pG4CS_count/sum(entire_g4s_rdna_summary2$pG4CS_count)) %>% 
  mutate(norm_pG4CS_count= round(norm_pG4CS_count, 2))

fwrite(entire_g4s_rdna_summary2, "pG4CS_KY962518_added_3500nt_IGS_upstream_no_junctn_strandwise_graphinput.csv")


pG4CS_strandwise<- ggplot(entire_g4s_rdna_summary2, aes(x= rDNA_region, y = norm_pG4CS_count, fill= strand)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized pG4CS strandwise distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalized pG4CS count", 
       fill= "pG4CS strand")+
  scale_y_continuous(breaks= seq(0, 0.30, by = 0.1), limits =c(0,0.30))+
  geom_text(aes(label= pG4CS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  scale_fill_manual(values= c("+" = "royalblue", "-" = "maroon3"), 
                    labels = c("+" = "Non-template", "-" = "Template"))+
  #scale_fill_manual(values = combined_colors)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))
#coord_flip()

ggsave( "Normalized_strandwise_pG4CS_distribution_in_human_rDNA_subcomponents_after_rule.tiff", 
        plot = pG4CS_strandwise, width=18,height=10, dpi=150)


nontemplate<- entire_g4s_rdna_summary2 %>% filter(strand == "+")
nontemplate$rDNA_region <- factor(nontemplate$rDNA_region, 
                                  levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                             "ITS2","28S", "3'ETS", "IGS" ))


template<- entire_g4s_rdna_summary2 %>% filter(strand == "-")
template$rDNA_region <- factor(template$rDNA_region, 
                               levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                          "ITS2","28S", "3'ETS", "IGS" ))



g4s_nontemplate <- ggplot(nontemplate, aes(x= rDNA_region, y = norm_pG4CS_count, fill= rDNA_region)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized Non-template pG4CS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalized Non-template pG4CS count", 
       fill= "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.30, by = 0.1), limits =c(0,0.30))+
  geom_text(aes(label= pG4CS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  
  scale_fill_manual(values= c("#F5FEFB", "#E21515", "#5AAA46", "#F36017", "#6B1519", 
                              "#818689", "#ECE612", "#E07F80", "#DE9A22"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))


ggsave( "Normalized_nontemplate_pG4CS_distribution_in_human_rDNA_subcomponents_AR.tiff", 
        plot = g4s_nontemplate, width=18,height=10, dpi=150)

g4s_template <- ggplot(template, aes(x= rDNA_region, y = norm_pG4CS_count, fill= rDNA_region)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized Template RLFS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalized Template RLFS count", 
       fill= "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.30, by = 0.1), limits =c(0,0.30))+
  geom_text(aes(label= pG4CS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  
  scale_fill_manual(values= c("#F5FEFB", "#E21515", "#5AAA46", "#F36017", "#6B1519", 
                              "#818689", "#ECE612", "#E07F80", "#DE9A22"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))


ggsave( "Normalized_template_pG4CS_distribution_in_human_rDNA_subcomponents_AR.tiff", 
        plot = g4s_template, width=18,height=10, dpi=150)



