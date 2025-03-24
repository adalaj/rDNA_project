
#To make plot for G4S distribution in human rDNA locus including and excluding promoter


library(tidyverse)
library(data.table)


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/g4s_and_rdna/human")

g4s_rdna<- fread("output_pG4CS_KY962518_1_humanrDNA_2018_nontemplate.txt", header = FALSE, sep = "\t")

g4s_rdna$v7 <- "junction"
g4s_rdna$v8<- g4s_rdna$V3-g4s_rdna$V2


g4s_rdna$v7[g4s_rdna$V3 < 3657] <- "5'ETS"
sum(g4s_rdna$v7=="5'ETS")
#23 ##here, we see identified 23 pG4CS in 5'ETS, but in powerpoint we can only see 18, this is because few G4S are in close vicinity and appear darker. 
#Darker lines indicate more than one pG4CS.
g4s_rdna$v7[g4s_rdna$V2 < 3657 & g4s_rdna$V3 > 3657 ] <- "5'ETS and 18S junction"


g4s_rdna$v7[g4s_rdna$V2 > 3657 & g4s_rdna$V3 < 5526 ] <- "18S"
g4s_rdna$v7[g4s_rdna$V2 < 5526 & g4s_rdna$V3 > 5526 ] <- "18S and ITS1 junction"


g4s_rdna$v7[g4s_rdna$V2 > 5526 & g4s_rdna$V3 < 6596 ] <- "ITS1"
g4s_rdna$v7[g4s_rdna$V2 < 6596 & g4s_rdna$V3 > 6596] <- "ITS1 and 5.8S junction"

g4s_rdna$v7[g4s_rdna$V2 > 6596 & g4s_rdna$V3 < 6753 ] <- "5.8S"
g4s_rdna$v7[g4s_rdna$V2 < 6753 & g4s_rdna$V3 > 6753] <- "5.8S and ITS2 junction"

g4s_rdna$v7[g4s_rdna$V2 > 6753 & g4s_rdna$V3 < 7920 ] <- "ITS2"
g4s_rdna$v7[g4s_rdna$V2 < 7920 & g4s_rdna$V3 > 7920 ] <- "ITS2 and 28S junction"

g4s_rdna$v7[g4s_rdna$V2 > 7920 & g4s_rdna$V3 < 12971 ] <- "28S"
g4s_rdna$v7[g4s_rdna$V2 < 12971 & g4s_rdna$V3 > 12971 ] <- "28S and 3'ETS junction"

g4s_rdna$v7[g4s_rdna$V2 > 12971 & g4s_rdna$V3 < 13332 ] <- "3'ETS"
g4s_rdna$v7[g4s_rdna$V2 < 13332 & g4s_rdna$V3 > 13332  ] <- "3'ETS and IGS junction"

g4s_rdna$v7[g4s_rdna$V2 > 13332 & g4s_rdna$V3 < 44838 ] <- "IGS"

g4s_rdna<- g4s_rdna %>% mutate(v9 = V3-V2) # calculate R-loop lengthcolnames(g4s_rdna)<- c("GenBank_Accession", "pG4CS_start", "pG4CS_end", "details", "score", "strand", "rDNA_region", "pG4CS_length")
colnames(g4s_rdna) <- c("GenBank_Accession", "pG4CS_start", "pG4CS_end", "pG4CS_sequence", "pG4CS_name", "strand", "rDNA_region", "pG4CS_length")

fwrite(g4s_rdna, "rDNA_hg38_chr21_2018_pG4CS_at_junctn.csv", sep = ",")

g4s_rdna_summary<- g4s_rdna %>% group_by(rDNA_region) %>% dplyr::count()
#rows order is lost 

new_rows<- data.table(rDNA_region = c("5'ETS and 18S junction", "18S", "18S and ITS1 junction","ITS1 and 5.8S junction",
                                      "5.8S", "5.8S and ITS2 junction","ITS2 and 28S junction",
                                      "28S and 3'ETS junction","3'ETS and IGS junction"))

new_rows$n<- 0


g4s_rdna_summary<- rbind(g4s_rdna_summary, new_rows)

g4s_rdna_summary[,1]

#1 28S, 3'ETS,             5'ETS,               IGS,             ITS1                  
#6 ITS2, 5'ETS and 18S junction, 18S,        18S and ITS1 junction, ITS1 and 5.8S junction
#11 5.8S                  
#12 5.8S and ITS2 junction
#13 ITS2 and 28S junction 
#14 28S and 3'ETS junction
#15 3'ETS and IGS junction

g4s_rdna_summary1 <- g4s_rdna_summary [c(3,7,8,9,5,10,11,12, 6, 13, 1,14, 2, 15, 4),]

names(g4s_rdna_summary1)[2] <- "pG4CS_count"


g4s_rdna_summary1<- g4s_rdna_summary1 %>% mutate(norm_pG4CS_count= pG4CS_count/sum(g4s_rdna_summary1$pG4CS_count)) %>% 
    mutate(norm_pG4CS_count = round(norm_pG4CS_count,2))

fwrite(g4s_rdna_summary1, "rDNA_hg38_chr21_2018_pG4CS_at_junctn_graphinput.csv", sep = ",")
fwrite(g4s_rdna, "rDNA_hg38_chr21_2018_pG4CS_at_junctn.csv", sep = ",")




g4s_rdna_summary1$rDNA_region <- factor(g4s_rdna_summary1$rDNA_region, 
                                   levels = c("5'ETS", "5'ETS and 18S junction", "18S", "18S and ITS1 junction", 
                                              "ITS1", "ITS1 and 5.8S junction", "5.8S", "5.8S and ITS2 junction", 
                                              "ITS2", "ITS2 and 28S junction","28S", "28S and 3'ETS junction", 
                                              "3'ETS", "3'ETS and IGS junction", "IGS" ))



library(RColorBrewer)

# Extract colors from both Set3 and Set1 palettes
colors_set3 <- brewer.pal(12, "Set3") # Set3 has 12 colors
colors_set1 <- brewer.pal(9, "Set1")  # Set1 has 9 colors

combined_colors <- c(colors_set1,colors_set3)

rDNA_pG4CS_graph<- ggplot(g4s_rdna_summary1, aes(x= rDNA_region, y = pG4CS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "pG4CS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "pG4CS count")+
  scale_y_continuous(breaks= seq(0, 90, by = 10), limits =c(0,90))+
  geom_text(aes(label= pG4CS_count, vjust= -0.5))+
  scale_fill_manual(values = combined_colors)+
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


ggsave( "pG4CS distribution in rDNA subcomponents.tiff", 
        plot = rDNA_pG4CS_graph, width=15,height=10, dpi=300)

rDNA_npG4CS_graph<- ggplot(g4s_rdna_summary1, aes(x= rDNA_region, y = norm_pG4CS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized pG4CS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalized pG4CS count")+
  scale_y_continuous(breaks= seq(0, 0.50, by = 0.1), limits =c(0,0.50))+
  geom_text(aes(label= pG4CS_count, vjust= -0.5))+
  scale_fill_manual(values = combined_colors)+
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

ggsave( "Normalized pG4CS distribution in rDNA subcomponents.tiff", 
        plot = rDNA_npG4CS_graph, width=20,height=10, dpi=300)





#repeat the same steps for promoter as well
promoter_g4s_rdna<- fread("output_pG4CS_KY962518_inclu_2kb_promoter_nontemplate.txt", header = FALSE, sep = "\t")
promoter_g4s_rdna$v7 <- "junction"

promoter_g4s_rdna$v8<- promoter_g4s_rdna$V3-promoter_g4s_rdna$V2

promoter_g4s_rdna$v7[promoter_g4s_rdna$V3 < 2000] <- "Promoter"
sum(promoter_g4s_rdna$v7=="Promoter")
#12

promoter_g4s_rdna$v7[promoter_g4s_rdna$V3 > 2000 & promoter_g4s_rdna$V2< 2000] <- "Promoter and 5'ETS junction"


promoter_g4s_rdna$v7[promoter_g4s_rdna$V2 > 2000 & promoter_g4s_rdna$V3< 5657] <- "5'ETS"

promoter_g4s_rdna$v7[promoter_g4s_rdna$V3 > 5657 & promoter_g4s_rdna$V2 < 5657 ] <- "5'ETS and 18S junction"

promoter_g4s_rdna$v7[promoter_g4s_rdna$V2 > 5658 & promoter_g4s_rdna$V3 < 7526 ] <- "18S"
promoter_g4s_rdna$v7[promoter_g4s_rdna$V3 > 7526 & promoter_g4s_rdna$V2 < 7526 ] <- "18S and ITS1 junction"


promoter_g4s_rdna$v7[promoter_g4s_rdna$V2 > 7527 & promoter_g4s_rdna$V3 < 8596 ] <- "ITS1"
promoter_g4s_rdna$v7[promoter_g4s_rdna$V3 > 6596 & promoter_g4s_rdna$V2 < 6596 ] <- "ITS1 and 5.8S junction"

promoter_g4s_rdna$v7[promoter_g4s_rdna$V2 > 8597 & promoter_g4s_rdna$V3 < 8753 ] <- "5.8S"
promoter_g4s_rdna$v7[promoter_g4s_rdna$V3 > 8753 & promoter_g4s_rdna$V2 < 8753 ] <- "5.8S and ITS2 junction"

promoter_g4s_rdna$v7[promoter_g4s_rdna$V2 > 8754 & promoter_g4s_rdna$V3 < 9920 ] <- "ITS2"
promoter_g4s_rdna$v7[promoter_g4s_rdna$V3 > 9920 & promoter_g4s_rdna$V2 < 9920 ] <- "ITS2 and 28S junction"

promoter_g4s_rdna$v7[promoter_g4s_rdna$V2 > 9921 & promoter_g4s_rdna$V3 < 14971 ] <- "28S"
promoter_g4s_rdna$v7[promoter_g4s_rdna$V3 > 14971 & promoter_g4s_rdna$V2 < 14971 ] <- "28S and 3'ETS junction"

promoter_g4s_rdna$v7[promoter_g4s_rdna$V2 > 14972 & promoter_g4s_rdna$V3 < 15332 ] <- "3'ETS"
promoter_g4s_rdna$v7[promoter_g4s_rdna$V3 > 15332 & promoter_g4s_rdna$V2 < 15332 ] <- "3'ETS and IGS junction"

promoter_g4s_rdna$v7[promoter_g4s_rdna$V2 > 15333 & promoter_g4s_rdna$V3 < 46838 ] <- "IGS"


colnames(promoter_g4s_rdna) <- c("GenBank_Accession", "pG4CS_start", "pG4CS_end", "pG4CS_sequence", "pG4CS_name", "strand", "rDNA_region", "pG4CS_length")

fwrite(promoter_g4s_rdna, "rDNA_hg38_chr21_2018_incld_promoter_pG4CS_at_junctn.csv", sep = ",")

promoter_g4s_rdna_summary<- promoter_g4s_rdna %>% group_by(rDNA_region) %>% count()
#rows order is lost 
promoter_g4s_rdna_summary[,1]

#1 28S         3'ETS       5'ETS      
#4 IGS          ITS1       ITS2       
#7 Promoter


new_rows<- data.table(rDNA_region = c("Promoter and 5'ETS junction", "5'ETS and 18S junction", "18S", 
                                      "18S and ITS1 junction","ITS1 and 5.8S junction",
                                      "5.8S", "5.8S and ITS2 junction","ITS2 and 28S junction",
                                      "28S and 3'ETS junction","3'ETS and IGS junction"))

new_rows$n<- 0


promoter_g4s_rdna_summary<- rbind(promoter_g4s_rdna_summary, new_rows)

promoter_g4s_rdna_summary <- promoter_g4s_rdna_summary[c(7,8,3,9,10,11,5,12,13,14,6,15,1,16,2,17,4),]

names(promoter_g4s_rdna_summary)[2] <- "pG4CS_count"


promoter_g4s_rdna_summary<- promoter_g4s_rdna_summary %>% mutate(norm_pG4CS_count = pG4CS_count/sum(promoter_g4s_rdna_summary$pG4CS_count)) %>% 
  mutate(norm_pG4CS_count= round(norm_pG4CS_count, 2))

fwrite(promoter_g4s_rdna_summary, "rDNA_hg38_chr21_2018_incld_promoter_pG4CS_at_junctn_graphinput.csv", sep = ",")

promoter_g4s_rdna_summary$rDNA_region <- factor(promoter_g4s_rdna_summary$rDNA_region, 
                                            levels = c("Promoter","Promoter and 5'ETS junction", "5'ETS", "5'ETS and 18S junction", 
                                                       "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                       "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                       "28S and 3'ETS junction", "3'ETS", "3'ETS and IGS junction", "IGS" ))


colors_set3 <- brewer.pal(12, "Set3") # Set3 has 12 colors
colors_set1 <- brewer.pal(9, "Set1")  # Set1 has 9 colors
newcolors<- c("mintcream", "violetred4")
combined_colors <- c(newcolors,colors_set1,colors_set3)

promoter_rDNA_pG4CS_graph<- ggplot(promoter_g4s_rdna_summary, aes(x= rDNA_region, y = pG4CS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "pG4CS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "pG4CS count")+
  scale_y_continuous(breaks= seq(0, 90, by = 10), limits =c(0,90))+
  geom_text(aes(label= pG4CS_count, vjust= -0.5))+
  scale_fill_manual(values = combined_colors)+
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

ggsave( "pG4CS distribution in rDNA subcomponents incld promoters.tiff", 
        plot = promoter_rDNA_pG4CS_graph, width=15,height=10, dpi=300)


promoter_rDNA_npG4CS_graph<- ggplot(promoter_g4s_rdna_summary, aes(x= rDNA_region, y = norm_pG4CS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized pG4CS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalized pG4CS count")+
  scale_y_continuous(breaks= seq(0, 0.50, by = 0.1), limits =c(0,0.50))+
  geom_text(aes(label= pG4CS_count, vjust= -0.5))+
  scale_fill_manual(values = combined_colors)+
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

ggsave( "Normalized pG4CS distribution in rDNA subcomponents incld promoters.tiff", 
        plot = promoter_rDNA_npG4CS_graph, width=15,height=10, dpi=300)


##As there is no pG4CS in juncion, i decided to remove all teh junctions rows 

filtered_data <- promoter_g4s_rdna_summary %>%
       filter(!grepl("junction", rDNA_region))

fwrite(filtered_data, "rDNA_hg38_chr21_2018_incld_promoter_pG4CS_no_junctn_graphinput.csv", sep = ",")


filtered_data$rDNA_region <- factor(filtered_data$rDNA_region, 
                                                levels = c("Promoter","5'ETS",  
                                                           "18S", "ITS1",  "5.8S",  "ITS2", 
                                                           "28S",  "3'ETS", "IGS" ))


no_junction_promoter_rDNA_npG4CS_graph<- ggplot(filtered_data, aes(x= rDNA_region, y = norm_pG4CS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized pG4CS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalized pG4CS count")+
  scale_y_continuous(breaks= seq(0, 0.50, by = 0.1), limits =c(0,0.50))+
  geom_text(aes(label= pG4CS_count, vjust= -0.5))+
  scale_fill_manual(values= c("mintcream", "red", "green", "orange", "brown", "pink", "lightyellow", "salmon", "peachpuff"))+
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

ggsave( "Normalized pG4CS distribution in rDNA subcomponents incld promoters but no junctn.tiff", 
        plot = no_junction_promoter_rDNA_npG4CS_graph, width=15,height=10, dpi=300)


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

entire_g4s_rdna<- entire_g4s_rdna[entire_g4s_rdna$V2 >1299 & entire_g4s_rdna$V2 <46137]#210 .. total length is 48338

entire_g4s_rdna$v7 <- "junction"

entire_g4s_rdna$v8<- entire_g4s_rdna$V3-entire_g4s_rdna$V2

entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 1299 & entire_g4s_rdna$V2 < 3500] <- "Promoter"
sum(entire_g4s_rdna$v7=="Promoter")
#12

entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 3500] <- "Promoter and 5'ETS junction"


entire_g4s_rdna$v7[entire_g4s_rdna$V2 >= 3501] <- "5'ETS"

entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 7157 ] <- "5'ETS and 18S junction"

entire_g4s_rdna$v7[entire_g4s_rdna$V2 >= 7158 ] <- "18S"
entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 9026 ] <- "18S and ITS1 junction"


entire_g4s_rdna$v7[entire_g4s_rdna$V2 >= 9027  ] <- "ITS1"
entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 10096] <- "ITS1 and 5.8S junction"

entire_g4s_rdna$v7[entire_g4s_rdna$V2 >= 10097 ] <- "5.8S"
entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 10253 ] <- "5.8S and ITS2 junction"

entire_g4s_rdna$v7[entire_g4s_rdna$V2 >= 10254] <- "ITS2"
entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 11420 ] <- "ITS2 and 28S junction"

entire_g4s_rdna$v7[entire_g4s_rdna$V2 >= 11421 ] <- "28S"
entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 16471 ] <- "28S and 3'ETS junction"

entire_g4s_rdna$v7[entire_g4s_rdna$V2 >= 16472 ] <- "3'ETS"
entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 16832 ] <- "3'ETS and IGS junction"

entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 16833 & entire_g4s_rdna$V2 < 46137] <- "IGS" 


colnames(entire_g4s_rdna) <- c("GenBank_Accession", "pG4CS_start", "pG4CS_end", "pG4CS_sequence", "pG4CS_name", "strand", "rDNA_region", "pG4CS_length")

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

entire_g4s_rdna_summary <- entire_g4s_rdna_summary[c(7,8,3,9,10,11,5,12,13,14,6,15, 1,16,2,17,4),]

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
  theme(axis.text.x = element_blank())+ 
  theme(panel.grid = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        theme(panel.grid = element_blank()),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +  # Center Y-axis title
  theme(axis.ticks.y = element_line(color = "black"))+
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
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ 
  theme(panel.grid = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        theme(panel.grid = element_blank()),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +  # Center Y-axis title
  theme(axis.ticks.y = element_line(color = "black"))
#coord_flip()

ggsave( "Normalized_strandwise_pG4CS_distribution_in_human_rDNA_subcomponents_after_rule.tiff", 
        plot = pG4CS_strandwise, width=18,height=10, dpi=150)



