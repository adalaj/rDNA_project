setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rloop_and_rdna/one_rDNA_seq/input")
library(Biostrings)
library(stringr)
library(data.table)
library(tidyverse)


##I have created my own excel sheet manually filling nums and normalised
without_rule<- fread("rdna_2018_without_rule_and_junctn.csv", sep = ",", header = TRUE)

without_rule$rDNA_region <- factor(without_rule$rDNA_region, 
                                   levels = c("5'ETS", "18S", "ITS1", "5.8S", "ITS2","28S", "3'ETS", "IGS" ))


rDNA_nrlfs_graph<- ggplot(without_rule, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalised RLFS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalised RLFS count")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.60))+
  geom_text(aes(label= RLFS_count, vjust= -0.5))+
  scale_fill_manual(values= c("red", "green", "orange", "brown", "pink", "lightyellow", "salmon", "peachpuff"))+
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

ggsave( "Normalised RLFS distribution in rDNA subcomponents without rule.tiff", 
        plot = rDNA_nrlfs_graph, width=20,height=10, dpi=150)



##Bruce suggested that take entire 2018 chr21 fasta file, look for RLFS in the entire one rDNA region
## then as you already know the boundary of each components, separate those RLFS according to their boundary 
## and you will also retain the rlfs determined at the boundary of two components for eg btw 5'ETS and 18S
## show the results in terms of bar graph and/ or pie chart

##steps
##open KY962518_1_humanrDNA_2018.fasta file in sublime text, 
##Only retain >KY962518_1_humanrDNA_2018.fasta and removed .1 Homo sapiens external transcribed spacer 18S ribosomal RNA gene, internal transcribed spacer 1, 5.8S ribosomal RNA gene, internal transcribed spacer 2, 28S ribosomal RNA gene, external transcribed spacer, and CDC27-like protein pseudogene, complete sequence
##saved this file as KY962518.txt file in downloads, 
## open terminal, cd to downloads and run this code
#python QmRLFS-finder.py -bed -i KY962518_1_humanrDNA_2018.fasta -o KY962518_1_humanrDNA_2018_qmrlfs

#KY962518_1_humanrDNA_2018_qmrlfs.out.bed is the output file saved under, one_rdna_seq>output

##the bed gives start and end of RLFS, and in name provide the total length of RLFS along with m1,m2 model details


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rloop_and_rdna/human/one_rDNA_seq/output/Qmrlfs_results_2018")
rdna_2018<- fread("KY962518_1_humanrDNA_2018_qmrlfs.out.bed", header = FALSE, sep = "\t")
rdna_2018_clean<- rdna_2018 %>% select(1:6)

rdna_2018_clean$v7 <- "junction"
rdna_2018_clean$v8<- rdna_2018_clean$V3-rdna_2018_clean$V2

rdna_2018_clean_v1<- rdna_2018_clean

rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V3 < 3657] <- "5'ETS"
sum(rdna_2018_clean_v1$v7=="5'ETS")
#38 ##here, we see identified 38 RLFS in 5'ETS, compared to 41 RLFS in 5'ETS when boundary were defined
## suggesting that few RLFS are formed at junction 
rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V2 < 3657 & rdna_2018_clean_v1$V3 > 3657 ] <- "5'ETS and 18S junction"


#3657+1869=5526
rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V2 > 3657 & rdna_2018_clean_v1$V3 < 5526 ] <- "18S"
rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V2 < 5526 & rdna_2018_clean_v1$V3 > 5526 ] <- "18S and ITS1 junction"

#5526+1070=6596
rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V2 > 5526 & rdna_2018_clean_v1$V3 < 6596 ] <- "ITS1"
rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V2 < 6596 & rdna_2018_clean_v1$V3 > 6596] <- "ITS1 and 5.8S junction"

rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V2 > 6596 & rdna_2018_clean_v1$V3 < 6753 ] <- "5.8S"
rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V2 < 6753 & rdna_2018_clean_v1$V3 > 6753] <- "5.8S and ITS2 junction"

rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V2 > 6753 & rdna_2018_clean_v1$V3 < 7920 ] <- "ITS2"
rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V2 < 7920 & rdna_2018_clean_v1$V3 > 7920 ] <- "ITS2 and 28S junction"

rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V2 > 7920 & rdna_2018_clean_v1$V3 < 12971 ] <- "28S"
rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V2 < 12971 & rdna_2018_clean_v1$V3 > 12971 ] <- "28S and 3'ETS junction"

rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V2 > 12971 & rdna_2018_clean_v1$V3 < 13332 ] <- "3'ETS"
rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V2 < 13332 & rdna_2018_clean_v1$V3 > 13332  ] <- "3'ETS and IGS junction"

rdna_2018_clean_v1$v7[rdna_2018_clean_v1$V2 > 13332 & rdna_2018_clean_v1$V3 < 44838 ] <- "IGS"

rdna_2018_clean_v1<- rdna_2018_clean_v1 %>% mutate(v9 = V3-V2) # calculate R-loop lengthcolnames(rdna_2018_clean_v1)<- c("GenBank_Accession", "RLFS_start", "RLFS_end", "details", "score", "strand", "rDNA_region", "RLFS_length")
colnames(rdna_2018_clean_v1) <- c("GenBank_Accession", "RLFS_start", "RLFS_end", "details", "score", "strand", "rDNA_region", "RLFS_length")

fwrite(rdna_2018_clean_v1, "rDNA_hg38_chr21_2018_RLFS_at_junctn.csv", sep = ",")

rdna_summary<- rdna_2018_clean_v1 %>% group_by(rDNA_region) %>% dplyr::count()
#rows order is lost 
rdna_summary[,1]

#1= 28S, 2= 28S and 3'ETS junction, 3= 3'ETS and IGS junction, 4= 5'ETS 
#5=5'ETS and 18S junction, 6=5.8S and ITS2 junction, 7=IGS
#8= ITS1, 9= ITS2, 10= ITS2 and 28S junction

rdna_summary <- rdna_summary [c(4,5,8,6,9,10,1,2,3,7),]

names(rdna_summary)[2] <- "RLFS_count"

rdna_summary<- rdna_summary %>% mutate(norm_RLFS_count = RLFS_count/sum(RLFS_count)) %>% 
  mutate(norm_RLFS_count= round(norm_RLFS_count, 2))

fwrite(rdna_summary, "rDNA_hg38_chr21_2018_RLFS_at_junctn_graphinput.csv", sep = ",")


##add 18s and junction as zero 
#Manually added zero in junction_graphinput


rdna_summary<- fread("rDNA_hg38_chr21_2018_RLFS_at_junctn_graphinput.csv", sep = ",", header = TRUE)
rdna_summary$rDNA_region <- factor(rdna_summary$rDNA_region, 
                                   levels = c("5'ETS", "5'ETS and 18S junction", "18S", "18S and ITS1 junction", 
                                              "ITS1", "ITS1 and 5.8S junction", "5.8S", "5.8S and ITS2 junction", 
                                              "ITS2", "ITS2 and 28S junction","28S", "28S and 3'ETS junction", 
                                              "3'ETS", "3'ETS and IGS junction", "IGS" ))

library(RColorBrewer)

# Extract colors from both Set3 and Set1 palettes
colors_set3 <- brewer.pal(12, "Set3") # Set3 has 12 colors
colors_set1 <- brewer.pal(9, "Set1")  # Set1 has 9 colors

combined_colors <- c(colors_set1,colors_set3)

rDNA_rlfs_graph<- ggplot(rdna_summary, aes(x= rDNA_region, y = RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "RLFS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "RLFS count")+
  scale_y_continuous(breaks= seq(0, 60, by = 10), limits =c(0,60))+
  geom_text(aes(label= RLFS_count, vjust= -0.5))+
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


ggsave( "RLFS distribution in rDNA subcomponents.tiff", 
        plot = rDNA_rlfs_graph, width=15,height=10, dpi=600)

rDNA_nrlfs_graph<- ggplot(rdna_summary, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalised RLFS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalised RLFS count")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.60))+
  geom_text(aes(label= RLFS_count, vjust= -0.5))+
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

ggsave( "Normalised RLFS distribution in rDNA subcomponents.tiff", 
        plot = rDNA_nrlfs_graph, width=20,height=10, dpi=600)




##i need to update the steps:
## to include the promoter which is 2kb upstream of 5'ets
##I first blast the 5'ETS in the ucsc genome browser 
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rloop_and_rdna/one_rDNA_seq/output/Qmrlfs_results_2018")
rdna_2018_promoter<- fread("rDNA_promoter_hg38_qmrlfs.out.bed", header = FALSE, sep = "\t")
rdna_2018_promoter<- rdna_2018_promoter %>% select(1:6)

##steps
##open KY962518_1_humanrDNA_2018.fasta file in sublime text, 
##Only retain >KY962518_1_humanrDNA_2018.fasta and removed .1 Homo sapiens external transcribed spacer 18S ribosomal RNA gene, internal transcribed spacer 1, 5.8S ribosomal RNA gene, internal transcribed spacer 2, 28S ribosomal RNA gene, external transcribed spacer, and CDC27-like protein pseudogene, complete sequence
##saved this file as KY962518.txt file in downloads, 



## open terminal, cd to downloads and run this code for bed files
# conda activate python2.7
#python QmRLFS-finder.py -bed -i KY962518_inclu_2kb_promoter.txt -o KY962518_inclu_2kb_promoter_qmrlfs

#and this to get table version
#python QmRLFS-finder.py -i KY962518_inclu_2kb_promoter.txt -o KY962518_inclu_2kb_promoter_qmrlfs

#KY962518_inclu_2kb_promoter_qmrlfs.out.bed is the output file saved under, one_rdna_seq>output>Qmrlfs_results_2018

##the bed gives start and end of RLFS, and in name provide the total length of RLFS along with m1,m2 model details

promoter_rdna_2018<- fread("KY962518_inclu_2kb_promoter_qmrlfs.out.bed", header = FALSE, sep = "\t")
promoter_rdna_2018_clean<- promoter_rdna_2018 %>% select(1:6)


promoter_rdna_2018_clean$v8 <- "junction"

promoter_rdna_2018_clean_v1<- promoter_rdna_2018_clean


promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V3 < 2000] <- "Promoter"
sum(promoter_rdna_2018_clean_v1$v8=="Promoter")
#8 ##  ##here, we see identified 8 RLFS in promoter, compared to 13 RLFS in promoter when boundary were defined look rDNA_promoter_hg38_qmrlfs.out.bed
## suggesting that few RLFS are formed at junction

promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V3 > 2000 & promoter_rdna_2018_clean_v1$V2< 2000] <- "Promoter and 5'ETS junction"
sum(promoter_rdna_2018_clean_v1$v8=="Promoter and 5'ETS junction")
#5  adding from 8 plus 5 is equal to 13 which we find in rDNA_promoter_hg38_qmrlfs.out.bed


promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V2 > 2000 & promoter_rdna_2018_clean_v1$V3< 5657] <- "5'ETS"
sum(promoter_rdna_2018_clean_v1$v8=="5'ETS")

#38 ##here, we see identified 38 RLFS in 5'ETS, compared to 41 RLFS in 5'ETS when boundary were defined
## suggesting that few RLFS are formed at junction 
promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V3 > 5657 & promoter_rdna_2018_clean_v1$V2 < 5657 ] <- "5'ETS and 18S junction"

promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V2 > 5658 & promoter_rdna_2018_clean_v1$V3 < 7526 ] <- "18S"
promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V3 > 7526 & promoter_rdna_2018_clean_v1$V2 < 7526 ] <- "18S and ITS1 junction"


promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V2 > 7527 & promoter_rdna_2018_clean_v1$V3 < 8596 ] <- "ITS1"
promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V3 > 6596 & promoter_rdna_2018_clean_v1$V2 < 6596 ] <- "ITS1 and 5.8S junction"

promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V2 > 8597 & promoter_rdna_2018_clean_v1$V3 < 8753 ] <- "5.8S"
promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V3 > 8753 & promoter_rdna_2018_clean_v1$V2 < 8753 ] <- "5.8S and ITS2 junction"

promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V2 > 8754 & promoter_rdna_2018_clean_v1$V3 < 9920 ] <- "ITS2"
promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V3 > 9920 & promoter_rdna_2018_clean_v1$V2 < 9920 ] <- "ITS2 and 28S junction"

promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V2 > 9921 & promoter_rdna_2018_clean_v1$V3 < 14971 ] <- "28S"
promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V3 > 14971 & promoter_rdna_2018_clean_v1$V2 < 14971 ] <- "28S and 3'ETS junction"

promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V2 > 14972 & promoter_rdna_2018_clean_v1$V3 < 15332 ] <- "3'ETS"
promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V3 > 15332 & promoter_rdna_2018_clean_v1$V2 < 15332 ] <- "3'ETS and IGS junction"

promoter_rdna_2018_clean_v1$v8[promoter_rdna_2018_clean_v1$V2 > 15333 & promoter_rdna_2018_clean_v1$V3 < 46838 ] <- "IGS"

promoter_rdna_2018_clean_v1<- promoter_rdna_2018_clean_v1 %>% mutate(v9 = V3-V2) # calculate R-loop lengthcolnames(rdna_2018_clean_v1)<- c("GenBank_Accession", "RLFS_start", "RLFS_end", "details", "score", "strand", "rDNA_region", "RLFS_length")

colnames(promoter_rdna_2018_clean_v1) <- c("GenBank_Accession", "RLFS_start", "RLFS_end", "details", "score", "strand", "rDNA_region", "RLFS_length")

fwrite(promoter_rdna_2018_clean_v1, "rDNA_hg38_chr21_2018_incld_promoter_RLFS_at_junctn.csv", sep = ",")

promoter_rdna_summary<- promoter_rdna_2018_clean_v1 %>% group_by(rDNA_region) %>% count()
#rows order is lost 
promoter_rdna_summary[,1]

#1= 28S, 2= 28S and 3'ETS junction, 3= 3'ETS and IGS junction, 4= 5'ETS 
#5=5'ETS and 18S junction, 6=5.8S and ITS2 junction, 7=IGS
#8= ITS1, 9= ITS2, 10= ITS2 and 28S junction, 11= Promoter, 12= Promoter and 5'ETS junction

promoter_rdna_summary <- promoter_rdna_summary [c(11, 12, 4,5,8,6,9,10,1,2,3,7),]

names(promoter_rdna_summary)[2] <- "RLFS_count"
fwrite(promoter_rdna_summary, "rDNA_hg38_chr21_2018_incld_promoter_RLFS_at_junctn_summary.csv", sep = ",")

promoter_rdna_summary<- promoter_rdna_summary %>% mutate(norm_RLFS_count = RLFS_count/sum(RLFS_count)) %>% 
                       mutate(norm_RLFS_count= round(norm_RLFS_count, 2))

fwrite(promoter_rdna_summary, "rDNA_hg38_chr21_2018_incld_promoter_RLFS_at_junctn_graphinput.csv", sep = ",")

promoter_rdna_summary$rDNA_region <- factor(promoter_rdna_summary$rDNA_region, 
                                   levels = c("Promoter", "Promoter and 5'ETS junction","5'ETS", "5'ETS and 18S junction", "18S", 
                                              "18S and ITS1 junction","ITS1", "ITS1 and 5.8S junction","5.8S" , "5.8S and ITS2 junction", "ITS2", 
                                              "ITS2 and 28S junction","28S", "28S and 3'ETS junction", "3'ETS",
                                              "3'ETS and IGS junction", "IGS" ))


colors_set3 <- brewer.pal(12, "Set3") # Set3 has 12 colors
colors_set1 <- brewer.pal(9, "Set1")  # Set1 has 9 colors
newcolors<- c("mintcream", "violetred4")
combined_colors <- c(newcolors,colors_set1,colors_set3)

promoter_rDNA_rlfs_graph<- ggplot(promoter_rdna_summary, aes(x= rDNA_region, y = RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "RLFS distribution in rDNA subcomponents", 
       x= "Human rDNA region", 
       y= "RLFS count")+
  scale_y_continuous(breaks= seq(0, 60, by = 10), limits =c(0,60))+
  geom_text(aes(label= RLFS_count, vjust= -0.5))+
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


ggsave( "RLFS distribution in rDNA subcomponents incld promoters.tiff", 
        plot = promoter_rDNA_rlfs_graph, width=15,height=10, dpi=300)


promoter_rDNA_nrlfs_graph<- ggplot(promoter_rdna_summary, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalised RLFS distribution in rDNA subcomponents", 
       x= "Human rDNA region", 
       y= "Normalized RLFS count")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.60))+
  geom_text(aes(label= RLFS_count, vjust= -0.5))+
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

ggsave( "Normalised RLFS distribution in rDNA subcomponents incld promoters.tiff", 
        plot = promoter_rDNA_nrlfs_graph, width=15,height=10, dpi=300)




##i want plot bar graph after the rule for RLFSs that has 3500 bp added to 5ETS.  
##basically input file will be output_RLFS_KY962518_added_3500nt_IGS_upstream_humanrDNA.txt after passing through QmRLFS algorithm.

##the rule: Counting the presence of RLFSs where it is first detected. For example, if RLFSs start towards the
# end of 5'ETS but stretches to 18S then it will counted under 5'ETS. 
# if we would have defined boundary in the first place then we would have lost this junction in our counting. 


entire_RLFSs_rdna<- fread("KY962518_added_3500nt_IGS_upstream_qmrlfs.out.bed", sep = "\t", header = FALSE) #208
#208 is total number of RLFS are found when both  the promoter (2202bp) and some  part of IGS (1297) are present twice
# this means after doing junction analysis we should have less than 208 
# to avoid double count i will filter the KY962518_added_3500nt_IGS_upstream_qmrlfs.out.bed right after the the actual IGS length ends that is 46136.

# as from the bar graph visualization it is clear that we have 13 extra RLFS from promoter. 

entire_RLFSs_rdna<- entire_RLFSs_rdna %>% select(1:6)
colnames(entire_RLFSs_rdna) <- c("GenBank_Accession", "RLFS_start", "RLFS_end", "RLFS_name","score", "strand")

# strand specificty doesnt matter here because we are cropping region of our interest
entire_RLFSs_rdna<- entire_RLFSs_rdna[entire_RLFSs_rdna$RLFS_start >1299 & entire_RLFSs_rdna$RLFS_start < 46137] #195 
entire_RLFSs_rdna$rDNA_region <- "junction"



#strand specificity will matter here now we want to allocate RLFS to rdna region sub components  based on their direction
# to do that, i will be creating a column that will have actual RLFS start based on strand specificity

entire_RLFSs_rdna<- entire_RLFSs_rdna %>% mutate(actual_rlfs_start = ifelse(entire_RLFSs_rdna$strand == "+", RLFS_start, RLFS_end))
entire_RLFSs_rdna$RLFS_length<- abs(entire_RLFSs_rdna$RLFS_start-entire_RLFSs_rdna$RLFS_end)

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 1299 & entire_RLFSs_rdna$actual_rlfs_start < 3500] <- "Promoter"
sum(entire_RLFSs_rdna$rDNA_region=="Promoter")#13

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 3500] <- "Promoter and 5'ETS junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 3501] <- "5'ETS"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 7157] <- "5'ETS and 18S junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 7158] <- "18S"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 9026] <- "18S and ITS1 junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 9027] <- "ITS1"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 10096] <- "ITS1 and 5.8S junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 10097] <- "5.8S"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 10253] <- "5.8S and ITS2 junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 10254] <- "ITS2"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 11420] <- "ITS2 and 28S junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 11421] <- "28S"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 16471] <- "28S and 3'ETS junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 16472] <- "3'ETS"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 16832] <- "3'ETS and IGS junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 16833 & entire_RLFSs_rdna$actual_rlfs_start < 46137] <- "IGS"


fwrite(entire_RLFSs_rdna, "RLFS_KY962518_added_3500nt_IGS_upstream_at_junctn_details_after_rule.csv", sep = ",")

entire_RLFSs_rdna_summary<- entire_RLFSs_rdna %>% group_by(rDNA_region) %>% count()
#rows order is lost 

sum(entire_RLFSs_rdna_summary$n)
#[1] 195 , here we see that 13 extra RLFS from promoter are removed.

entire_RLFSs_rdna_summary[,1]

#1 18S      28S         3'ETS            
#4 5'ETS    IGS          ITS1       
#7 ITS2     Promoter


new_rows<- data.table(rDNA_region = c("Promoter and 5'ETS junction", "5'ETS and 18S junction", 
                                      "18S and ITS1 junction","ITS1 and 5.8S junction",
                                      "5.8S", "5.8S and ITS2 junction","ITS2 and 28S junction",
                                      "28S and 3'ETS junction", "3'ETS and IGS junction"))

new_rows$n<- 0


entire_RLFSs_rdna_summary<- rbind(entire_RLFSs_rdna_summary, new_rows)
entire_RLFSs_rdna_summary$rDNA_region
#[1] "18S"                         "28S"                        
#[3] "3'ETS"                       "5'ETS"                      
#[5] "IGS"                         "ITS1"                       
#[7] "ITS2"                        "Promoter"                   
#[9] "Promoter and 5'ETS junction" "5'ETS and 18S junction"     
#[11] "18S and ITS1 junction"       "ITS1 and 5.8S junction"     
#[13] "5.8S"                        "5.8S and ITS2 junction"     
#[15] "ITS2 and 28S junction"       "28S and 3'ETS junction"     
#[17] "3'ETS and IGS junction" 


row.names(entire_RLFSs_rdna_summary)<- entire_RLFSs_rdna_summary$rDNA_region

entire_RLFSs_rdna_summary <- entire_RLFSs_rdna_summary[c("Promoter","Promoter and 5'ETS junction", "5'ETS", "5'ETS and 18S junction", 
                                                           "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                           "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                           "28S and 3'ETS junction", "3'ETS", "3'ETS and IGS junction", "IGS" ),]

names(entire_RLFSs_rdna_summary)[2] <- "RLFS_count"


entire_RLFSs_rdna_summary<- entire_RLFSs_rdna_summary %>% mutate(norm_RLFS_count = RLFS_count/sum(entire_RLFSs_rdna_summary$RLFS_count)) %>% 
  mutate(norm_RLFS_count= round(norm_RLFS_count, 2))

fwrite(entire_RLFSs_rdna_summary, "RLFS_KY962518_added_3500nt_IGS_upstream_at_junctn_after_rule_graphinput.csv", sep = ",")

entire_RLFSs_rdna_summary$rDNA_region <- factor(entire_RLFSs_rdna_summary$rDNA_region, 
                                              levels = c("Promoter","Promoter and 5'ETS junction", "5'ETS", "5'ETS and 18S junction", 
                                                         "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                         "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                         "28S and 3'ETS junction", "3'ETS", "3'ETS and IGS junction", "IGS" ))


RLFS_norm_3500igs<- ggplot(entire_RLFSs_rdna_summary, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized RLFS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalized RLFS count", 
       fill = "rDNA region")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.60))+
  geom_text(aes(label= RLFS_count), vjust= -0.5, size= 5)+
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
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),   # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))
#coord_flip()

ggsave( "Normalized_RLFS_distribution_in_human_rDNA_subcomponents_incld_junctn_AR.tiff", 
        plot = RLFS_norm_3500igs, width=18,height=10, dpi=150) #AR is after rule, bcoz powerpoint is not accepting too long image name


RLFSs_rdna_summary<- entire_RLFSs_rdna_summary[!grepl("junction", entire_RLFSs_rdna_summary$rDNA_region),]


#When you use coord_flip(), the order of the factor levels in rDNA_region determines how the categories are displayed along the flipped y-axis.
#By default, the first factor level appears at the bottom when flipped, and the last factor level appears at the top.

RLFSs_rdna_summary$rDNA_region <- factor(RLFSs_rdna_summary$rDNA_region, 
                                                 levels = rev(c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                            "ITS2","28S", "3'ETS", "IGS" )))



#To reverse the order so that "Promoter" appears at the top when flipped, modify the levels of the factor like this



RLFS_norm_3500igs_nojuntn<- ggplot(RLFSs_rdna_summary, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized RLFS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalized RLFS count", 
       fill = "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.60))+
  geom_text(aes(label= RLFS_count, hjust= -1.0, vjust= 0.5, size= 50))+
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
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),   # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))+
  coord_flip()

ggsave("Normalized_RLFS_distribution_in_human_rDNA_subcomponents_after_rule.tiff", 
        plot = RLFS_norm_3500igs_nojuntn, width=18,height=10, dpi=150)



#to make template and non-template
entire_RLFSs_rdna_summary2<- entire_RLFSs_rdna %>% group_by(rDNA_region, strand) %>% count()
names(entire_RLFSs_rdna_summary2)[3] <- "RLFS_count"

new_rows<- data.table(rDNA_region = c("18S", "5.8S","5.8S"),
                      strand= c("+", "+", "-"),
                      RLFS_count = c(0, 0,0))

entire_RLFSs_rdna_summary2<- rbind(entire_RLFSs_rdna_summary2, new_rows)
entire_RLFSs_rdna_summary2$rDNA_region <- factor(entire_RLFSs_rdna_summary2$rDNA_region, 
                                                levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                        "ITS2","28S", "3'ETS", "IGS" ))
                                                            
entire_RLFSs_rdna_summary2<- entire_RLFSs_rdna_summary2 %>% mutate(norm_RLFS_count = RLFS_count/sum(entire_RLFSs_rdna_summary2$RLFS_count)) %>% 
  mutate(norm_RLFS_count= round(norm_RLFS_count, 2))

fwrite(entire_RLFSs_rdna_summary2, "RLFS_KY962518_added_3500nt_IGS_upstream_no_junctn_strandwise_AR_graphinput.csv")


rlfs_strandwise<- ggplot(entire_RLFSs_rdna_summary2, aes(x= rDNA_region, y = norm_RLFS_count, fill= strand)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized RLFS strandwise distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalized RLFS count", 
       fill= "RLFS strand")+
  scale_y_continuous(breaks= seq(0, 0.30, by = 0.1), limits =c(0,0.30))+
  geom_text(aes(label= RLFS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  scale_fill_manual(values= c("+" = "royalblue", "-" = "maroon3"), 
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

        

ggsave( "Normalized_strandwise_RLFS_distribution_in_human_rDNA_subcomponents_AR.tiff", 
        plot = rlfs_strandwise, width=18,height=10, dpi=150)




nontemplate<- entire_RLFSs_rdna_summary2 %>% filter(strand == "+")
nontemplate$rDNA_region <- factor(nontemplate$rDNA_region, 
                                                 levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                            "ITS2","28S", "3'ETS", "IGS" ))


template<- entire_RLFSs_rdna_summary2 %>% filter(strand == "-")
template$rDNA_region <- factor(template$rDNA_region, 
                                                 levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                            "ITS2","28S", "3'ETS", "IGS" ))



rlfs_nontemplate <- ggplot(nontemplate, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized Non-template RLFS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalized Non-template RLFS count", 
       fill= "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.30, by = 0.1), limits =c(0,0.30))+
  geom_text(aes(label= RLFS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  
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


ggsave( "Normalized_nontemplate_RLFS_distribution_in_human_rDNA_subcomponents_AR.tiff", 
        plot = rlfs_nontemplate, width=18,height=10, dpi=150)

rlfs_template <- ggplot(template, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized Template RLFS distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalized Template RLFS count", 
       fill= "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.30, by = 0.1), limits =c(0,0.30))+
  geom_text(aes(label= RLFS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  
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


ggsave( "Normalized_template_RLFS_distribution_in_human_rDNA_subcomponents_AR.tiff", 
        plot = rlfs_template, width=18,height=10, dpi=150)



