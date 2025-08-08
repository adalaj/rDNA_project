##i want plot bar graph after the rule for RLFSs that has 5000 bp added to 5ETS.  
##basically input file will be output_RLFS_KY962518_added_3500nt_IGS_upstream_humanrDNA.txt after passing through QmRLFS algorithm.

##the rule: Counting the presence of RLFSs where it is first detected. For example, if RLFSs start towards the
# end of 5'ETS but stretches to 18S then it will counted under 5'ETS. 
# if we would have defined boundary in the first place then we would have lost this junction during counting. 



setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/mouse/output/QmRLFS")

library(data.table)
library(tidyverse)


#(python2.7) jyotiadala@Jyotis-MacBook-Pro Downloads % python QmRLFS-finder.py -bed -i BK000964_added_5000nt_IGS_upstream_nontemplate.fasta -o BK000964_added_5000nt_IGS_upstream_nontemplate_qmrlfs
#QmRLFS-finder.py (version v1.5)
#run on Fri Apr 04 2025 13:14:27 
#command line: python QmRLFS-finder.py -bed -i BK000964_added_5000nt_IGS_upstream_nontemplate.fasta -o BK000964_added_5000nt_IGS_upstream_nontemplate_qmrlfs

#Time used: 0.36 mins



#(python2.7) jyotiadala@Jyotis-MacBook-Pro Downloads % python QmRLFS-finder.py -i BK000964_added_5000nt_IGS_upstream_nontemplate.fasta -o BK000964_added_5000nt_IGS_upstream_nontemplate_qmrlfs 
#QmRLFS-finder.py (version v1.5)
#run on Fri Apr 04 2025 13:15:12 
#command line: python QmRLFS-finder.py -i BK000964_added_5000nt_IGS_upstream_nontemplate.fasta -o BK000964_added_5000nt_IGS_upstream_nontemplate_qmrlfs


entire_RLFSs_rdna<- fread("BK000964_added_5000nt_IGS_upstream_nontemplate_qmrlfs.out.bed", sep = "\t", header = FALSE) #71
entire_RLFSs_rdna<- entire_RLFSs_rdna %>% select(1:6)
colnames(entire_RLFSs_rdna) <- c("GenBank_Accession", "RLFS_start", "RLFS_end", "RLFS_name","score", "strand")

# strand specificity doesnt matter here because we are cropping region of our interest
entire_RLFSs_rdna<- entire_RLFSs_rdna[entire_RLFSs_rdna$RLFS_start >2001 & entire_RLFSs_rdna$RLFS_start < 47307] #68
#After 2000 is where the promoter begins and 45338 is where IGS end (50306-promoter length which is 3000 = 47306)

entire_RLFSs_rdna$rDNA_region <- "junction"



#strand specificity will matter here now we want to allocate RLFS to rdna region sub components  based on their direction
# to do that, i will be creating a column that will have actual RLFS start based on strand specificity

entire_RLFSs_rdna<- entire_RLFSs_rdna %>% mutate(actual_rlfs_start = ifelse(entire_RLFSs_rdna$strand == "+", RLFS_start, RLFS_end))
entire_RLFSs_rdna$RLFS_length<- abs(entire_RLFSs_rdna$RLFS_start-entire_RLFSs_rdna$RLFS_end)


entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 2001 & entire_RLFSs_rdna$actual_rlfs_start < 5001] <- "Promoter"
sum(entire_RLFSs_rdna$rDNA_region=="Promoter")
#3

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 5000] <- "Promoter and 5'ETS junction"


entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 5001] <- "5'ETS"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 9007 ] <- "5'ETS and 18S junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 9008] <- "18S"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 10877 ] <- "18S and ITS1 junction"


entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 10878   ] <- "ITS1"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 11877] <- "ITS1 and 5.8S junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 11878 ] <- "5.8S"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 12034 ] <- "5.8S and ITS2 junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 12035] <- "ITS2"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 13122] <- "ITS2 and 28S junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 13123 ] <- "28S"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 17852 ] <- "28S and 3'ETS junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start >= 17853] <- "3'ETS"
entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 18403 ] <- "3'ETS and IGS junction"

entire_RLFSs_rdna$rDNA_region[entire_RLFSs_rdna$actual_rlfs_start > 18404 & entire_RLFSs_rdna$actual_rlfs_start < 45339] <- "IGS"


fwrite(entire_RLFSs_rdna, "RLFS_BK000964_added_5000nt_IGS_upstream_at_junctn_details.csv", sep = ",")

entire_RLFSs_rdna_summary<- entire_RLFSs_rdna %>% group_by(rDNA_region) %>% count()
#rows order is lost 

sum(entire_RLFSs_rdna_summary$n)
#[1] 68

entire_RLFSs_rdna_summary[,1]

#1 28S         3'ETS       5'ETS       
#4 IGS          ITS1       ITS2       promoter


new_rows<- data.table(rDNA_region = c("Promoter and 5'ETS junction", "5'ETS and 18S junction", "18S", 
                                      "18S and ITS1 junction","ITS1 and 5.8S junction","5.8S",
                                      "5.8S and ITS2 junction","ITS2 and 28S junction",
                                      "28S and 3'ETS junction", "3'ETS and IGS junction"))

new_rows$n<- 0


entire_RLFSs_rdna_summary<- rbind(entire_RLFSs_rdna_summary, new_rows)
entire_RLFSs_rdna_summary$rDNA_region
#[1] "28S"                         "3'ETS"                      
#[3] "5'ETS"                       "IGS"                        
#[5] "ITS1"                        "ITS2"                       
#[7] "Promoter"                    "Promoter and 5'ETS junction"
#[9] "5'ETS and 18S junction"      "18S"                        
#[11] "18S and ITS1 junction"       "ITS1 and 5.8S junction"     
#[13] "5.8S"                         "5.8S and ITS2 junction"     
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

fwrite(entire_RLFSs_rdna_summary, "RLFS_BK000964_added_5000nt_IGS_upstream_at_junctn_graphinput.csv", sep = ",")

entire_RLFSs_rdna_summary$rDNA_region <- factor(entire_RLFSs_rdna_summary$rDNA_region, 
                                                levels = c("Promoter","Promoter and 5'ETS junction", "5'ETS", "5'ETS and 18S junction", 
                                                           "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                           "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                           "28S and 3'ETS junction", "3'ETS", "3'ETS and IGS junction", "IGS" ))


RLFS_norm_5000igs<- ggplot(entire_RLFSs_rdna_summary, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized RLFS distribution in the Mouse rDNA locus", 
       x= "Mouse rDNA region", 
       y= "Normalized RLFS count")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.60))+
  geom_text(aes(label= RLFS_count, vjust= -0.5, size= 50))+
  scale_fill_manual(values= c( "#B6FFF4","maroon", "#FDCCE5", "steelblue", "#D0B6FF","darkviolet", "#EF9B20","burlywood2", "#A0322B", 
                               "pink", "#FFCC17","aquamarine", "#E5FFB6","greenyellow", "#3B8CC4","turquoise2", "#A4A2A8"))+
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

ggsave( "Normalized_RLFS_distribution_in_mouse_rDNA_subcomponents_incld_junctn_AR.tiff", 
        plot = RLFS_norm_3500igs, width=18,height=10, dpi=150)






RLFSs_rdna_summary<- entire_RLFSs_rdna_summary[!grepl("junction", entire_RLFSs_rdna_summary$rDNA_region),]


#When you use coord_flip(), the order of the factor levels in rDNA_region determines how the categories are displayed along the flipped y-axis.
#By default, the first factor level appears at the bottom when flipped, and the last factor level appears at the top.

RLFSs_rdna_summary$rDNA_region <- factor(RLFSs_rdna_summary$rDNA_region, 
                                         levels = rev(c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                        "ITS2","28S", "3'ETS", "IGS" )))



#To reverse the order so that "Promoter" appears at the top when flipped, modify the levels of the factor like this



RLFS_norm_5000igs_nojuntn<- ggplot(RLFSs_rdna_summary, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized RLFS distribution in the Mouse rDNA locus", 
       x= "Mouse rDNA region", 
       y= "Normalized RLFS count", 
       fill = "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.60))+
  geom_text(aes(label= RLFS_count, hjust= -1.0, vjust= 0.5, size= 50))+
  scale_fill_manual(values= rev(c("#B6FFF4", "#FDCCE5", "#D0B6FF", "#EF9B20", "#A0322B", 
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

ggsave("Normalized_RLFS_distribution_in_mouse_rDNA_subcomponents_after_rule.tiff", 
       plot = RLFS_norm_5000igs_nojuntn, width=18,height=10, dpi=300)



#to make template and non-template
entire_RLFSs_rdna_summary2<- entire_RLFSs_rdna %>% group_by(rDNA_region, strand) %>% count()
names(entire_RLFSs_rdna_summary2)[3] <- "RLFS_count"

new_rows<- data.table(rDNA_region = c("Promoter","3'ETS","18S","18S","5.8S","5.8S"),
                      strand= c("-","-","+","-", "+", "-"),
                      RLFS_count = c(0, 0,0,0,0,0))

entire_RLFSs_rdna_summary2<- rbind(entire_RLFSs_rdna_summary2, new_rows)
entire_RLFSs_rdna_summary2$rDNA_region <- factor(entire_RLFSs_rdna_summary2$rDNA_region, 
                                                 levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                            "ITS2","28S", "3'ETS", "IGS" ))

entire_RLFSs_rdna_summary2<- entire_RLFSs_rdna_summary2 %>% mutate(norm_RLFS_count = RLFS_count/sum(entire_RLFSs_rdna_summary2$RLFS_count)) %>% 
  mutate(norm_RLFS_count= round(norm_RLFS_count, 2))

fwrite(entire_RLFSs_rdna_summary2, "RLFS_KY962518_added_5000nt_IGS_upstream_no_junctn_strandwise_AR_graphinput.csv")


rlfs_strandwise<- ggplot(entire_RLFSs_rdna_summary2, aes(x= rDNA_region, y = norm_RLFS_count, fill= strand)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized RLFS strandwise distribution in the mouse rDNA locus", 
       x= "Mouse rDNA region", 
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



ggsave( "Normalized_strandwise_RLFS_distribution_in_mouse_rDNA_subcomponents_AR.tiff", 
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
  labs(title= "Normalized Non-template RLFS distribution in the Mouse rDNA locus", 
       x= "Mouse rDNA region", 
       y= "Normalized Non-template RLFS count", 
       fill= "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.40, by = 0.1), limits =c(0,0.40))+
  geom_text(aes(label= RLFS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  
  scale_fill_manual(values= c("#B6FFF4", "#FDCCE5", "#D0B6FF", "#EF9B20", "#A0322B", 
                              "#FFCC17", "#E5FFB6", "#3B8CC4", "#A4A2A8"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))


ggsave( "Normalized_nontemplate_RLFS_distribution_in_mouse_rDNA_subcomponents_AR.tiff", 
        plot = rlfs_nontemplate, width=18,height=10, dpi=150)

rlfs_template <- ggplot(template, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized Template RLFS distribution in the Mouse rDNA locus", 
       x= "Mouse rDNA region", 
       y= "Normalized Template RLFS count", 
       fill= "rDNA")+
  scale_y_continuous(breaks= seq(0, 0.40, by = 0.1), limits =c(0,0.40))+
  geom_text(aes(label= RLFS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  
  scale_fill_manual(values= c("#B6FFF4", "#FDCCE5", "#D0B6FF", "#EF9B20", "#A0322B", 
                              "#FFCC17", "#E5FFB6", "#3B8CC4", "#A4A2A8"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))


ggsave( "Normalized_template_RLFS_distribution_in_mouse_rDNA_subcomponents_AR.tiff", 
        plot = rlfs_template, width=18,height=10, dpi=150)




##R-loop visualization:
library(karyoploteR)

entire_rdna<- fread("BK000964_added_5000nt_IGS_upstream_nontemplate_qmrlfs.out.bed", sep = "\t", header = FALSE) #71
entire_rdna$V1= "rDNA_locus"
entire_rdna6<- entire_rdna %>% select(1:6)
colnames(entire_rdna6)<- c("chr", "start", "end", "name", "score", "strand")

##separate as per strand
entire_rdna6_nontemplate<- entire_rdna6 %>% filter(strand=="+") #46
#because in NCBI keep nontemplate sequence.


entire_rdna6_template<- entire_rdna6 %>% filter(strand=="-")#25


##plotting begins
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=50306))
#end is 50306 because I added 5000 to 45306.

##make separate data.table for each rdna components

#"yellowgreen", "red", "green", "orange", "brown", "pink", "lightyellow", "salmon", "peachpuff"
png("mouse_rdna_both_strand_rlfs.png", width = 50, height= 30, units= "in", res = 150)


kp <- plotKaryotype(genome=custom_genome, plot.type = 2)

kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =2000 , y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks last 2000bp from IGS representing previous rdna 
kpRect(kp, chr = 'rDNA_locus', x0 = 2001, x1 =5000 , y0 = 0, y1 = 1, col = "#B6FFF4", data.panel = "ideogram", borders= NA) #marks 3000 bp of  promoter

kpRect(kp, chr = 'rDNA_locus', x0 = 5001, x1 = 9007 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram", borders= NA) #marks 5'ETS
#5001+(4007-1) = 9007

kpRect(kp, chr = 'rDNA_locus', x0 = 9008, x1 = 10877, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram", borders= NA) #marks 18S
#9008+(1870-1) = 10877

kpRect(kp, chr = 'rDNA_locus', x0 = 10877, x1 = 11877, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram", borders= NA) #marks ITS1
#10877+(1000-1) = 11877

kpRect(kp, chr = 'rDNA_locus', x0 = 11878, x1 = 12034 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram", borders= NA) #marks 5.8S
#11878+(157-1) = 12034

kpRect(kp, chr = 'rDNA_locus', x0 = 12035, x1 = 13122, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram", borders= NA) #marks ITS2
#12035+(1088-1) = 13122

kpRect(kp, chr = 'rDNA_locus', x0 = 13123 , x1 = 17852, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram", borders= NA) #marks 28S
#13123+(4730-1)= 17852

kpRect(kp, chr = 'rDNA_locus', x0 = 17853, x1 = 18403, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram", borders= NA) #marks 3'ETS
#17853+(551-1) = 18403

kpRect(kp, chr = 'rDNA_locus', x0 = 18404, x1 = 50306, y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks IGS
#18404+(31903-1)= 50306


kpRect(kp, chr = 'rDNA_locus', x0 = 47306, x1 = 50306, y0 = 0, y1 = 1, col = "#B6FFF4", data.panel = "ideogram", borders= NA) #marks 2202bp of promoter
#50306-3000 = 47306

kpPlotRegions(kp, data=entire_rdna6_template, col="#1414E1", r0= -0.5, r1= -1.3)
kpPlotRegions(kp, data=entire_rdna6_nontemplate, col="#E21515", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width
dev.off()

#use zoom option, took screenshot and edited in powerpoint

##wanted to plot only till 3'ets and nontemplate
##plotting begins

png("mouse_rdna_nontemplate_rlfs_coverage.png", width = 30, height= 30, units= "in", res = 150)
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=20403)) 

##make separate data.table for each rdna components

#"yellowgreen", "red", "green", "orange", "brown", "pink", "lightyellow", "salmon", "peachpuff"

kp <- plotKaryotype(genome=custom_genome, plot.type = 2)

kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =2000 , y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks last 2000bp from IGS representing previous rdna 
kpRect(kp, chr = 'rDNA_locus', x0 = 2001, x1 =5000 , y0 = 0, y1 = 1, col = "#B6FFF4", data.panel = "ideogram", borders= NA) #marks 3000 bp of  promoter

kpRect(kp, chr = 'rDNA_locus', x0 = 5001, x1 = 9007 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram", borders= NA) #marks 5'ETS
#5001+(4007-1) = 9007

kpRect(kp, chr = 'rDNA_locus', x0 = 9008, x1 = 10877, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram", borders= NA) #marks 18S
#9008+(1870-1) = 10877

kpRect(kp, chr = 'rDNA_locus', x0 = 10877, x1 = 11877, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram", borders= NA) #marks ITS1
#10877+(1000-1) = 11877

kpRect(kp, chr = 'rDNA_locus', x0 = 11878, x1 = 12034 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram", borders= NA) #marks 5.8S
#11878+(157-1) = 12034

kpRect(kp, chr = 'rDNA_locus', x0 = 12035, x1 = 13122, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram", borders= NA) #marks ITS2
#12035+(1088-1) = 13122

kpRect(kp, chr = 'rDNA_locus', x0 = 13123 , x1 = 17852, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram", borders= NA) #marks 28S
#13123+(4730-1)= 17852

kpRect(kp, chr = 'rDNA_locus', x0 = 17853, x1 = 18403, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram", borders= NA) #marks 3'ETS
#17853+(551-1) = 18403

kpRect(kp, chr = 'rDNA_locus', x0 = 18404, x1 = 20403, y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks 3'ETS




kpPlotRegions(kp, data=entire_rdna6_nontemplate, col="#E21515", r0= -0.5, r1= -1.0) #-1.5 to make blue with more width

kpPlotCoverage(kp, data=entire_rdna6_nontemplate, col = "#E21515")
kpPlotRegions(kp, data=entire_rdna6_nontemplate, data.panel=2, col = "#E21515")
dev.off()


png("mouse_rdna_template_rlfs_coverage.png", width = 30, height= 30, units= "in", res = 150)
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=20403))
#end is 50306 because I added 5000 to 45306.

##make separate data.table for each rdna components

#"yellowgreen", "red", "green", "orange", "brown", "pink", "lightyellow", "salmon", "peachpuff"

kp <- plotKaryotype(genome=custom_genome, plot.type = 2)

kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =2000 , y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks last 2000bp from IGS representing previous rdna 
kpRect(kp, chr = 'rDNA_locus', x0 = 2001, x1 =5000 , y0 = 0, y1 = 1, col = "#B6FFF4", data.panel = "ideogram", borders= NA) #marks 3000 bp of  promoter

kpRect(kp, chr = 'rDNA_locus', x0 = 5001, x1 = 9007 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram", borders= NA) #marks 5'ETS
#5001+(4007-1) = 9007

kpRect(kp, chr = 'rDNA_locus', x0 = 9008, x1 = 10877, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram", borders= NA) #marks 18S
#9008+(1870-1) = 10877

kpRect(kp, chr = 'rDNA_locus', x0 = 10877, x1 = 11877, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram", borders= NA) #marks ITS1
#10877+(1000-1) = 11877

kpRect(kp, chr = 'rDNA_locus', x0 = 11878, x1 = 12034 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram", borders= NA) #marks 5.8S
#11878+(157-1) = 12034

kpRect(kp, chr = 'rDNA_locus', x0 = 12035, x1 = 13122, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram", borders= NA) #marks ITS2
#12035+(1088-1) = 13122

kpRect(kp, chr = 'rDNA_locus', x0 = 13123 , x1 = 17852, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram", borders= NA) #marks 28S
#13123+(4730-1)= 17852

kpRect(kp, chr = 'rDNA_locus', x0 = 17853, x1 = 18403, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram", borders= NA) #marks 3'ETS
#17853+(551-1) = 18403

kpRect(kp, chr = 'rDNA_locus', x0 = 18404, x1 = 20403, y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks 3'ETS


kpPlotRegions(kp, data=entire_rdna6_template, col="#1414E1", r0= -0.5, r1= -1.0) #-1.5 to make blue with more width

kpPlotCoverage(kp, data=entire_rdna6_template, col = "#1414E1")
kpPlotRegions(kp, data=entire_rdna6_template, data.panel=2, col = "#1414E1")
dev.off()




