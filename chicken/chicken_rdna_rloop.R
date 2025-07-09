##i want plot bar graph after to chicken rDNA sequence  
##basically input file will be KT445934_chicken_rDNA_2017.fasta for QmRLFS algorithm.

##the rule: Counting the presence of RLFSs where it is first detected. For example, if RLFSs start towards the
# end of 5'ETS but stretches to 18S then it will counted under 5'ETS. 
# if we would have defined boundary in the first place then we would have lost this junction during counting. 



setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/chicken")

library(data.table)
library(tidyverse)


#open terminal
#(base) jyotiadala@Jyotis-MacBook-Pro ~ % conda activate python2.7
#(python2.7) jyotiadala@Jyotis-MacBook-Pro ~ % cd Downloads 
#(python2.7) jyotiadala@Jyotis-MacBook-Pro Downloads % python QmRLFS-finder.py -bed -i KT445934_chicken_rDNA_2017.fasta -o KT445934_chicken_rDNA_2017_qmrlfs
#QmRLFS-finder.py (version v1.5)
#run on Fri Apr 18 2025 14:55:03 
#command line: python QmRLFS-finder.py -bed -i KT445934_chicken_rDNA_2017.fasta -o KT445934_chicken_rDNA_2017_qmrlfs
#Time used: 0.27 mins
#(python2.7) jyotiadala@Jyotis-MacBook-Pro Downloads % 


#(python2.7) jyotiadala@Jyotis-MacBook-Pro Downloads % python QmRLFS-finder.py -i KT445934_chicken_rDNA_2017.fasta -o KT445934_chicken_rDNA_2017_qmrlfs 
#QmRLFS-finder.py (version v1.5)
#run on Fri Apr 18 2025 14:56:29 
#command line: python QmRLFS-finder.py -i KT445934_chicken_rDNA_2017.fasta -o KT445934_chicken_rDNA_2017_qmrlfs
#Time used: 0.25 mins




entire_RLFSs_rdna<- fread("KT445934_chicken_rDNA_2017_qmrlfs.out.bed", sep = "\t", header = FALSE) #91
entire_RLFSs_rdna<- entire_RLFSs_rdna %>% select(1:6)
colnames(entire_RLFSs_rdna) <- c("GenBank_Accession", "RLFS_start", "RLFS_end", "RLFS_name","score", "strand")


entire_RLFSs_rdna$rDNA_region <- "junction"



#strand specificity will matter here now we want to allocate RLFS to rdna region sub components  based on their direction
# to do that, i will be creating a column that will have actual RLFS start based on strand specificity

entire_RLFSs_rdna<- entire_RLFSs_rdna %>% mutate(actual_rlfs_start = ifelse(entire_RLFSs_rdna$strand == "+", RLFS_start, RLFS_end))
entire_RLFSs_rdna<- entire_RLFSs_rdna %>% mutate(actual_rlfs_end = ifelse(entire_RLFSs_rdna$strand == "+", RLFS_end, RLFS_start))
entire_RLFSs_rdna$RLFS_length<- abs(entire_RLFSs_rdna$RLFS_start-entire_RLFSs_rdna$RLFS_end)


entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 1] <- "5'ETS"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 1836 ] <- "5'ETS and 18S junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 1837] <- "18S"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 3659 ] <- "18S and ITS1 junction"


entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 3660   ] <- "ITS1"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 6189] <- "ITS1 and 5.8S junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 6190 ] <- "5.8S"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 6346 ] <- "5.8S and ITS2 junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 6347] <- "ITS2"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 7079] <- "ITS2 and 28S junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 7080 ] <- "28S"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 11520 ] <- "28S and 3'ETS junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 11521] <- "3'ETS"


fwrite(entire_RLFSs_rdna, "RLFS_KT445934_chicken_junctn_details.csv", sep = ",")

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

fwrite(entire_RLFSs_rdna_summary, "RLFS_KT445934_chicken_at_junctn_graphinput.csv", sep = ",")

entire_RLFSs_rdna_summary$rDNA_region <- factor(entire_RLFSs_rdna_summary$rDNA_region, 
                                                levels = c("5'ETS", "5'ETS and 18S junction", 
                                                           "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                           "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                           "28S and 3'ETS junction", "3'ETS" ))


RLFS_norm<- ggplot(entire_RLFSs_rdna_summary, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized RLFS distribution in the chicken rDNA locus", 
       x= "Chicken rDNA region", 
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

ggsave( "Normalized_RLFS_distribution_in_chicken_rDNA_subcomponents_incld_junctn_AR.tiff", 
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
  labs(title= "Normalized RLFS distribution in the chicken rDNA locus", 
       x= "Chicken rDNA region", 
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

ggsave("Normalized_RLFS_distribution_in_chicken_rDNA_subcomponents_after_rule.tiff", 
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

fwrite(entire_RLFSs_rdna_summary2, "RLFS_KY962518_chicken_no_junctn_strandwise_AR_graphinput.csv")


rlfs_strandwise<- ggplot(entire_RLFSs_rdna_summary2, aes(x= rDNA_region, y = norm_RLFS_count, fill= strand)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized RLFS strandwise distribution in the chicken rDNA locus", 
       x= "Chicken rDNA region", 
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



ggsave( "Normalized_strandwise_RLFS_distribution_in_chicken_rDNA_subcomponents_AR.tiff", 
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
  labs(title= "Normalized Non-template RLFS distribution in the chicken rDNA locus", 
       x= "Chicken rDNA region", 
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


ggsave( "Normalized_nontemplate_RLFS_distribution_in_chicken_rDNA_subcomponents_AR.tiff", 
        plot = rlfs_nontemplate, width=18,height=10, dpi=150)

rlfs_template <- ggplot(template, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized Template RLFS distribution in the chicken rDNA locus", 
       x= "Chicken rDNA region", 
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


ggsave( "Normalized_template_RLFS_distribution_in_chicken_rDNA_subcomponents_AR.tiff", 
        plot = rlfs_template, width=18,height=10, dpi=150)




##R-loop visualization:
library(karyoploteR)

entire_rdna<- fread("KT445934_chicken_rDNA_2017_qmrlfs.out.bed", sep = "\t", header = FALSE) #91
entire_rdna$V1= "rDNA_locus"
entire_rdna6<- entire_rdna %>% select(1:6)
colnames(entire_rdna6)<- c("chr", "start", "end", "name", "score", "strand")

##separate as per strand
entire_rdna6_nontemplate<- entire_rdna6 %>% filter(strand=="+") #51
#because in NCBI keep nontemplate sequence.


entire_rdna6_template<- entire_rdna6 %>% filter(strand=="-")#40


##plotting begins
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=11863))
kp <- plotKaryotype(genome=custom_genome, plot.type = 2)


kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 = 1836 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram", borders= NA) #marks 5'ETS
kpRect(kp, chr = 'rDNA_locus', x0 = 1837, x1 = 3659, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram", borders= NA) #marks 18S
kpRect(kp, chr = 'rDNA_locus', x0 = 3660, x1 = 6189, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram", borders= NA) #marks ITS1
kpRect(kp, chr = 'rDNA_locus', x0 = 6190, x1 = 6346 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram", borders= NA) #marks 5.8S
kpRect(kp, chr = 'rDNA_locus', x0 = 6347, x1 = 7079, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram", borders= NA) #marks ITS2
kpRect(kp, chr = 'rDNA_locus', x0 = 7080 , x1 = 11520, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram", borders= NA) #marks 28S
kpRect(kp, chr = 'rDNA_locus', x0 = 11521, x1 = 11863, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram", borders= NA) #marks 3'ETS


kpPlotRegions(kp, data=entire_rdna6_template, col="#1414E1", r0= -0.5, r1= -1.3)
kpPlotRegions(kp, data=entire_rdna6_nontemplate, col="#E21515", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width


#Playing with colors use zoom option, took screenshot and edited in powerpoint

kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 = 1836 , y0 = 0, y1 = 1, col = "#F5F5DC", data.panel = "ideogram", borders= NA) #marks 5'ETS
kpRect(kp, chr = 'rDNA_locus', x0 = 1837, x1 = 3659, y0 = 0, y1 = 1, col = "#90EE90", data.panel = "ideogram", borders= NA) #marks 18S
kpRect(kp, chr = 'rDNA_locus', x0 = 3660, x1 = 6189, y0 = 0, y1 = 1, col = "#E6E6FA", data.panel = "ideogram", borders= NA) #marks ITS1
kpRect(kp, chr = 'rDNA_locus', x0 = 6190, x1 = 6346 , y0 = 0, y1 = 1, col = "#ADD8E6", data.panel = "ideogram", borders= NA) #marks 5.8S
kpRect(kp, chr = 'rDNA_locus', x0 = 6347, x1 = 7079, y0 = 0, y1 = 1, col = "#D3D3D3", data.panel = "ideogram", borders= NA) #marks ITS2
kpRect(kp, chr = 'rDNA_locus', x0 = 7080 , x1 = 11520, y0 = 0, y1 = 1, col = "#FFFFE0", data.panel = "ideogram", borders= NA) #marks 28S
kpRect(kp, chr = 'rDNA_locus', x0 = 11521, x1 = 11863, y0 = 0, y1 = 1, col = "#FFB6C1", data.panel = "ideogram", borders= NA) #marks 3'ETS


#Playing with colors use zoom option, took screenshot and edited in powerpoint
kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 = 1836 , y0 = 0, y1 = 1, col = "#FFB6C1", data.panel = "ideogram", borders= NA) #marks 5'ETS
kpRect(kp, chr = 'rDNA_locus', x0 = 1837, x1 = 3659, y0 = 0, y1 = 1, col = "#ADD8E6", data.panel = "ideogram", borders= NA) #marks 18S
kpRect(kp, chr = 'rDNA_locus', x0 = 3660, x1 = 6189, y0 = 0, y1 = 1, col = "#FFDAB9", data.panel = "ideogram", borders= NA) #marks ITS1
kpRect(kp, chr = 'rDNA_locus', x0 = 6190, x1 = 6346 , y0 = 0, y1 = 1, col = "#E6E6FA", data.panel = "ideogram", borders= NA) #marks 5.8S
kpRect(kp, chr = 'rDNA_locus', x0 = 6347, x1 = 7079, y0 = 0, y1 = 1, col = "#FFFFE0", data.panel = "ideogram", borders= NA) #marks ITS2
kpRect(kp, chr = 'rDNA_locus', x0 = 7080 , x1 = 11520, y0 = 0, y1 = 1, col = "#AAF0D1", data.panel = "ideogram", borders= NA) #marks 28S
kpRect(kp, chr = 'rDNA_locus', x0 = 11521, x1 = 11863, y0 = 0, y1 = 1, col = "#F7E7CE", data.panel = "ideogram", borders= NA) #marks 3'ETS



kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 = 1836 , y0 = 0, y1 = 1, col = "#FFB6C1", data.panel = "ideogram", borders= NA) #marks 5'ETS
kpRect(kp, chr = 'rDNA_locus', x0 = 1837, x1 = 3659, y0 = 0, y1 = 1, col = "#AAF0D1", data.panel = "ideogram", borders= NA) #marks 18S
kpRect(kp, chr = 'rDNA_locus', x0 = 3660, x1 = 6189, y0 = 0, y1 = 1, col = "#FFDAB9", data.panel = "ideogram", borders= NA) #marks ITS1
kpRect(kp, chr = 'rDNA_locus', x0 = 6190, x1 = 6346 , y0 = 0, y1 = 1, col = "#924045", data.panel = "ideogram", borders= NA) #marks 5.8S
kpRect(kp, chr = 'rDNA_locus', x0 = 6347, x1 = 7079, y0 = 0, y1 = 1, col = "#FFFFE0", data.panel = "ideogram", borders= NA) #marks ITS2
kpRect(kp, chr = 'rDNA_locus', x0 = 7080 , x1 = 11520, y0 = 0, y1 = 1, col = "#E6E6FA", data.panel = "ideogram", borders= NA) #marks 28S
kpRect(kp, chr = 'rDNA_locus', x0 = 11521, x1 = 11863, y0 = 0, y1 = 1, col = "#ADD8E6", data.panel = "ideogram", borders= NA) #marks 3'ETS

kpPlotRegions(kp, data=entire_rdna6_template, col="#1414E1", r0= -0.5, r1= -1.3)
kpPlotRegions(kp, data=entire_rdna6_nontemplate, col="#FDCCE5", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width



##i like these but bruce needs to confirm.

##template only
kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 = 1836 , y0 = 0, y1 = 1, col = "#FFB6C1", data.panel = "ideogram", borders= NA) #marks 5'ETS
kpRect(kp, chr = 'rDNA_locus', x0 = 1837, x1 = 3659, y0 = 0, y1 = 1, col = "#B5F2DC", data.panel = "ideogram", borders= NA) #marks 18S
kpRect(kp, chr = 'rDNA_locus', x0 = 3660, x1 = 6189, y0 = 0, y1 = 1, col = "#FFE0C2", data.panel = "ideogram", borders= NA) #marks ITS1
kpRect(kp, chr = 'rDNA_locus', x0 = 6190, x1 = 6346 , y0 = 0, y1 = 1, col = "#C98686", data.panel = "ideogram", borders= NA) #marks 5.8S
kpRect(kp, chr = 'rDNA_locus', x0 = 6347, x1 = 7079, y0 = 0, y1 = 1, col = "#FFFFE0", data.panel = "ideogram", borders= NA) #marks ITS2
kpRect(kp, chr = 'rDNA_locus', x0 = 7080 , x1 = 11520, y0 = 0, y1 = 1, col = "#E8E8FB", data.panel = "ideogram", borders= NA) #marks 28S
kpRect(kp, chr = 'rDNA_locus', x0 = 11521, x1 = 11863, y0 = 0, y1 = 1, col = "#A2CCE9", data.panel = "ideogram", borders= NA) #marks 3'ETS

kpPlotRegions(kp, data=entire_rdna6_template, col="#1414E1", r0= -0.5, r1= -1.3)



##non template only
kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 = 1836 , y0 = 0, y1 = 1, col = "#FFB6C1", data.panel = "ideogram", borders= NA) #marks 5'ETS
kpRect(kp, chr = 'rDNA_locus', x0 = 1837, x1 = 3659, y0 = 0, y1 = 1, col = "#B5F2DC", data.panel = "ideogram", borders= NA) #marks 18S
kpRect(kp, chr = 'rDNA_locus', x0 = 3660, x1 = 6189, y0 = 0, y1 = 1, col = "#FFE0C2", data.panel = "ideogram", borders= NA) #marks ITS1
kpRect(kp, chr = 'rDNA_locus', x0 = 6190, x1 = 6346 , y0 = 0, y1 = 1, col = "#C98686", data.panel = "ideogram", borders= NA) #marks 5.8S
kpRect(kp, chr = 'rDNA_locus', x0 = 6347, x1 = 7079, y0 = 0, y1 = 1, col = "#FFFFE0", data.panel = "ideogram", borders= NA) #marks ITS2
kpRect(kp, chr = 'rDNA_locus', x0 = 7080 , x1 = 11520, y0 = 0, y1 = 1, col = "#E8E8FB", data.panel = "ideogram", borders= NA) #marks 28S
kpRect(kp, chr = 'rDNA_locus', x0 = 11521, x1 = 11863, y0 = 0, y1 = 1, col = "#A2CCE9", data.panel = "ideogram", borders= NA) #marks 3'ETS

kpPlotRegions(kp, data=entire_rdna6_nontemplate, col="#FDCCE5", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width

"#DCDCDC","#FFB6C1","#D0B6FF","#E5FFB6","#FFE0C2","#B6FFF4","#FFFFE0", "#E8E8FB","#B6E5FF","#E9BFA2"

