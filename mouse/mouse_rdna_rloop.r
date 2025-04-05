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
  scale_fill_manual(values= rev(c("#F5FEFB", "#E21515", "#5AAA46", "#F36017", "#6B1519", 
                                  "#818689", "#ECE612", "#E07F80", "#DE9A22")))+
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
       plot = RLFS_norm_5000igs_nojuntn, width=18,height=10, dpi=150)



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


ggsave( "Normalized_template_RLFS_distribution_in_mouse_rDNA_subcomponents_AR.tiff", 
        plot = rlfs_template, width=18,height=10, dpi=150)





######All of these are incorrect due to wrong input entry

#here strand was not considered and promoter region was not considered
{

##I need to make graph with how many rlfs are found in each compartment

mouse_rdna_rlfs<- fread("rDNA_mouse_2013_qmrlfs.out.bed", header = FALSE, sep = "\t")
mouse_rdna_rlfs<- mouse_rdna_rlfs %>% select(V1, V2, V3, V4, V5, V6)

colnames(mouse_rdna_rlfs)<- c("rdna_region", "rlfs_start", "rlfs_end", "rlfs_name", "rlfs_score", "rlfs_strand")
mouse_rdna_rlfs_summary<- mouse_rdna_rlfs %>% 
  group_by(rdna_region) %>% 
  dplyr::count()

mouse_rdna_rlfs_summary<-separate(mouse_rdna_rlfs_summary, rdna_region, into = c("mouse_rdna_region", "org"), sep = "_")

mouse_rdna_rlfs_summary<- mouse_rdna_rlfs_summary %>% select(mouse_rdna_region, n)

colnames(mouse_rdna_rlfs_summary)<- c("mouse_rdna_region","rlfs_count")

no_rlfs<- data.frame(mouse_rdna_region = c("18S", "5.8S"), 
                     rlfs_count = c(0,0))
mouse_rdna_rlfs_summary<- rbind(mouse_rdna_rlfs_summary, no_rlfs)


##Find normalised frequency
mouse_rdna_rlfs_summary<- mouse_rdna_rlfs_summary %>% 
  mutate(norm_rlfs_count= rlfs_count/sum(rlfs_count)) %>% 
  mutate(norm_rlfs_count=round(norm_rlfs_count, 2))




mouse_rdna_rlfs_summary$mouse_rdna_region<- factor(mouse_rdna_rlfs_summary$mouse_rdna_region, 
                                             levels = c("5'ETS", "18S", "ITS1", "5.8S", 
                                                        "ITS2", "28S", "3'ETS", "IGS"))

fwrite (mouse_rdna_rlfs_summary, 
        file= "rDNA_BK000964_mouse_2013_RLFS_graphinput.csv")

mouse_rdna_rlfs_summary$rDNA_region <- factor(mouse_rdna_rlfs_summary$rDNA_region, 
                                   levels = c("5'ETS", "18S", "ITS1", "5.8S", "ITS2","28S", "3'ETS", "IGS" ))


stack<- ggplot(mouse_rdna_rlfs_summary, aes(rDNA_region, y=norm_RLFS_count,fill= rDNA_region))+
  geom_bar(stat = 'identity', color= "black")+
  labs(title= "Normalised RLFS distribution in the Mouse rDNA locus", 
       x= "mouse rDNA region", #(NCBI accession = BK000964), 
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

ggsave( "RLFS distribution in mouse rDNA subcomponents.tiff", 
        plot = stack, width=11,height=10, dpi=100)


norm_stack<- ggplot(mouse_rdna_rlfs_summary, aes(mouse_rdna_region, y=norm_rlfs_count,fill= mouse_rdna_region))+
  geom_bar(stat = 'identity', color= "black")+
  labs(title= "RLFS distribution in mouse rDNA subcomponents", 
       x= " mouse rDNA region (NCBI accession = BK000964)", 
       y= "RLFS count")+
  geom_text(aes(label= rlfs_count, vjust= -0.5))+
  scale_fill_brewer(palette = "Set3")+
  theme_minimal()+
  theme(axis.text.x = element_blank())+ 
  theme(panel.grid = element_blank(), 
        axis.line=element_line(color = "black"))

ggsave( "Normalised RLFS distribution in mouse rDNA subcomponents.tiff", 
        plot = norm_stack, width=11,height=10, dpi=600)
                                   

##for junction, 

##I REMOVED THE BK000964.3 TPA_exp: Mus musculus ribosomal DNA, complete repeating unit AND ONLY RETIANED BK000964 IN HEADER

#python QmRLFS-finder.py -i BK000964_mouse_rDNA_2013.fasta -o rDNA_BK000964_mouse_2013_qmrlfs
#python QmRLFS-finder.py -bed -i BK000964_mouse_rDNA_2013.fasta -o rDNA_BK000964_mouse_2013_qmrlfs

rdna_2013<- fread("rDNA_BK000964_mouse_2013_qmrlfs.out.bed", header = FALSE, sep = "\t")
rdna_2013_clean<- rdna_2013 %>% select(1:6)


rdna_2013_clean$v7 <- "junction"

rdna_2013_clean$v7[rdna_2013_clean$V3 < 4007] <- "5'ETS"
sum(rdna_2013_clean=="5'ETS")
#12 ##here, we see identified 12 RLFS in 5'ETS, compared to 41 RLFS in 5'ETS when boundary were defined
## suggesting that few RLFS are formed at junction 
rdna_2013_clean$v7[rdna_2013_clean$V3 > 4007 & rdna_2013_clean$V2 < 4007 ] <- "5'ETS and 18S junction"
sum(rdna_2013_clean=="5'ETS and 18S junction") #0

rdna_2013_clean$v7[rdna_2013_clean$V2 > 4008 & rdna_2013_clean$V3 < 5877 ] <- "18S"
sum(rdna_2013_clean=="18S") #0

rdna_2013_clean$v7[rdna_2013_clean$V3 > 5877 & rdna_2013_clean$V2 < 5877 ] <- "18S and ITS1 junction"
sum(rdna_2013_clean=="18S and ITS1 junction") #0

rdna_2013_clean$v7[rdna_2013_clean$V2 > 5877 & rdna_2013_clean$V3 < 6877 ] <- "ITS1"
sum(rdna_2013_clean=="ITS1") #4

rdna_2013_clean$v7[rdna_2013_clean$V3 > 6877 & rdna_2013_clean$V2 < 6877 ] <- "ITS1 and 5.8S junction"
sum(rdna_2013_clean=="ITS1 and 5.8S junction") #0

rdna_2013_clean$v7[rdna_2013_clean$V2 > 6877 & rdna_2013_clean$V3 < 7034 ] <- "5.8S"
sum(rdna_2013_clean=="5.8S") #0

rdna_2013_clean$v7[rdna_2013_clean$V3 > 7034 & rdna_2013_clean$V2 < 7034 ] <- "5.8S and ITS2 junction"
sum(rdna_2013_clean=="5.8S and ITS2 junction") #1


rdna_2013_clean$v7[rdna_2013_clean$V2 > 7034 & rdna_2013_clean$V3 < 8122 ] <- "ITS2"
sum(rdna_2013_clean=="ITS2") #7

rdna_2013_clean$v7[rdna_2013_clean$V3 > 8122 & rdna_2013_clean$V2 < 8122 ] <- "ITS2 and 28S junction"
sum(rdna_2013_clean=="ITS2 and 28S junction") #0

rdna_2013_clean$v7[rdna_2013_clean$V2 > 8122 & rdna_2013_clean$V3 < 12852 ] <- "28S"
sum(rdna_2013_clean=="28S") #36

rdna_2013_clean$v7[rdna_2013_clean$V3 > 12852 & rdna_2013_clean$V2 < 12852 ] <- "28S and 3'ETS junction"
sum(rdna_2013_clean=="28S and 3'ETS junction") #0

rdna_2013_clean$v7[rdna_2013_clean$V2 > 12852 & rdna_2013_clean$V3 < 13403 ] <- "3'ETS"
sum(rdna_2013_clean=="3'ETS") #2

rdna_2013_clean$v7[rdna_2013_clean$V3 > 13403 & rdna_2013_clean$V2 < 13403 ] <- "3'ETS and IGS junction"
sum(rdna_2013_clean=="3'ETS and IGS junction") #1


rdna_2013_clean$v7[rdna_2013_clean$V2 > 13403 & rdna_2013_clean$V3 < 45306 ] <- "IGS"
sum(rdna_2013_clean== "IGS") #5

rdna_2013_clean<- rdna_2013_clean %>% mutate(v8 = V3-V2) # calculate R-loop lengthcolnames(rdna_2018_clean_v1)<- c("GenBank_Accession", "RLFS_start", "RLFS_end", "details", "score", "strand", "rDNA_region", "RLFS_length")
colnames(rdna_2013_clean) <- c("GenBank_Accession", "RLFS_start", "RLFS_end", "details", "score", "strand", "rDNA_region", "RLFS_length")

fwrite(rdna_2013_clean, "rDNA_BK000964_mouse_2013_RLFS_at_junctn.csv", sep = ",")


rdna_summary<- rdna_2013_clean %>% group_by(rDNA_region) %>% dplyr::count()
no_rlfs<- data.frame(rDNA_region = c("18S","18S and ITS1 junction","ITS1 and 5.8S junction", 
                                     "5.8S", "ITS2 and 28S junction", "28S and 3'ETS junction", "5'ETS and 18S junction"), 
                     n = c(0,0))

rdna_summary<- rbind(rdna_summary, no_rlfs)
#rows order is lost 
rdna_summary[,1]

#1 28S, 
#2 3'ETS                 
#3 3'ETS and IGS junction
#4 5'ETS                 
# 5 5.8S and ITS2 junction
# 6 IGS                   
# 7 ITS1                  
# 8 ITS2                  
# 9 18S                   
# 10 18S and ITS1 junction 
# 11 ITS1 and 5.8S junction
# 12 5.8S                  
# 13 ITS2 and 28S junction 
# 14 28S and 3'ETS junction
# 15 5'ETS and 18S junction

rdna_summary <- rdna_summary [c(4,15,9,10, 7, 11, 12, 5, 8, 13, 1, 14, 2, 3, 6),]

names(rdna_summary)[2] <- "RLFS_count"


fwrite(rdna_summary, "rDNA_BK000964_mouse_2013_RLFS_at_junctn_graphinput.csv", sep = ",")
#Manually added zero 


rdna_summary<- fread("rDNA_BK000964_mouse_2013_RLFS_at_junctn_graphinput.csv", sep = ",", header = TRUE)

rdna_summary<- rdna_summary %>% mutate(norm_RLFS_count = RLFS_count/sum(RLFS_count)) %>% 
  mutate(norm_RLFS_count= round(norm_RLFS_count, 2))

fwrite(rdna_summary, "rDNA_BK000964_mouse_2013_RLFS_at_junctn_graphinput.csv")
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
  labs(title= "RLFS distribution in the Mouse rDNA locus", 
       x= "Mouse rDNA region", 
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
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5))+# Center Y-axis title
  theme(axis.ticks.y = element_line(color = "black"))
      

ggsave( "RLFS distribution in mouse rDNA subcomponents junctn.tiff", 
        plot = rDNA_rlfs_graph, width=15,height=11, dpi=600)

rDNA_nrlfs_graph<- ggplot(rdna_summary, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalised RLFS distribution in the Mouse rDNA locus", 
       x= "Mouse rDNA region", 
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
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +  # Center Y-axis title
  theme(axis.ticks.y = element_line(color = "black"))

ggsave( "Normalised RLFS distribution in mouse rDNA subcomponents junctn.tiff", 
        plot = rDNA_nrlfs_graph, width=20,height=11, dpi=600)



}


##plotting with IGS on start and promoter at end.
###Plot4 showing entire rDNA and begining of next but this is incorrect because the promoter RLFS was present in IGS


{
#open terminal
#(python3.11) jyotiadala@Jyotis-MacBook-Pro Downloads % conda activate python2.7
#(python2.7) jyotiadala@Jyotis-MacBook-Pro Downloads % python QmRLFS-finder.py -bed -i BK000964_added_3500nt_IGS_upstream_nontemplate.fasta -o BK000964_added_3500nt_IGS_upstream_qmrlfs
#QmRLFS-finder.py (version v1.5)
#run on Sun Feb 23 2025 05:01:50 
#command line: python QmRLFS-finder.py -bed -i BK000964_added_3500nt_IGS_upstream_nontemplate.fasta -o BK000964_added_3500nt_IGS_upstream_qmrlfs

#Time used: 0.36 mins

#read the rlfs that overlapped with rdna locus
entire_rdna<- fread("BK000964_added_3500nt_IGS_upstream_qmrlfs.out.bed", sep = "\t", header = FALSE) #71
entire_rdna$V1= "rDNA_locus"
entire_rdna6<- entire_rdna %>% select(1:6)
colnames(entire_rdna6)<- c("chr", "start", "end", "name", "score", "strand")

##separate as per strand
entire_rdna6_nontemplate<- entire_rdna6 %>% filter(strand=="+") #46
#because in NCBI keep nontemplate sequence.


entire_rdna6_template<- entire_rdna6 %>% filter(strand=="-")#25


##plotting begins
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=48806))
#end is 48806 because I added 3500 to 45306

##make separate data.table for each rdna components

#"yellowgreen", "red", "green", "orange", "brown", "pink", "lightyellow", "salmon", "peachpuff"

kp <- plotKaryotype(genome=custom_genome, plot.type = 2)

kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =1299 , y0 = 0, y1 = 1, col = "#DE9A22", data.panel = "ideogram", borders= NA) #marks last 1298bp from IGS representing previous rdna 
kpRect(kp, chr = 'rDNA_locus', x0 = 1300, x1 =3500 , y0 = 0, y1 = 1, col = "#F5FEFB", data.panel = "ideogram", borders= NA) #marks 2200 bp of  promoter

kpRect(kp, chr = 'rDNA_locus', x0 = 3501, x1 = 7507 , y0 = 0, y1 = 1, col = "#E21515", data.panel = "ideogram", borders= NA) #marks 5'ETS (3501+(3657-1))
#3501+(4007-1) = 7507

kpRect(kp, chr = 'rDNA_locus', x0 = 7508, x1 = 9377, y0 = 0, y1 = 1, col = "#5AAA46", data.panel = "ideogram", borders= NA) #marks 18S
#7508+(1870-1) = 9377

kpRect(kp, chr = 'rDNA_locus', x0 = 9378, x1 = 10377, y0 = 0, y1 = 1, col = "#F36017", data.panel = "ideogram", borders= NA) #marks ITS1
#9378+(1000-1) = 10377

kpRect(kp, chr = 'rDNA_locus', x0 = 10378, x1 = 10534, y0 = 0, y1 = 1, col = "#6B1519", data.panel = "ideogram", borders= NA) #marks 5.8S
#10378+(157-1) = 10534

kpRect(kp, chr = 'rDNA_locus', x0 = 10535, x1 = 11622, y0 = 0, y1 = 1, col = "#818689", data.panel = "ideogram", borders= NA) #marks ITS2
#10535+(1088-1) = 11622

kpRect(kp, chr = 'rDNA_locus', x0 = 11623, x1 = 16352, y0 = 0, y1 = 1, col = "#ECE612", data.panel = "ideogram", borders= NA) #marks 28S
#11623+(4730-1)= 16352

kpRect(kp, chr = 'rDNA_locus', x0 = 16353, x1 = 16903, y0 = 0, y1 = 1, col = "#E07F80", data.panel = "ideogram", borders= NA) #marks 3'ETS
#16353+(551-1) = 16903

kpRect(kp, chr = 'rDNA_locus', x0 = 16904, x1 = 48806, y0 = 0, y1 = 1, col = "#DE9A22", data.panel = "ideogram", borders= NA) #marks IGS
#16904+(31903-1)= 48338


kpRect(kp, chr = 'rDNA_locus', x0 = 46606, x1 = 48806, y0 = 0, y1 = 1, col = "#F5FEFB", data.panel = "ideogram", borders= NA) #marks 2202bp of promoter
#48806-2200

kpPlotRegions(kp, data=entire_rdna6_template, col="maroon3", r0= -0.5, r1= -1.3)
kpPlotRegions(kp, data=entire_rdna6_nontemplate, col="royalblue", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width


#use zoom option, took screenshot and edited in powerpoint




##i want plot bar graph after the rule for RLFSs that has 3500 bp added to 5ETS and also has promoter at the end. 

##basically input file will be output_RLFS_KY962518_added_3500nt_IGS_upstream_humanrDNA.txt

##teh rule: Counting the presence of RLFSs where it is first detected. For example, if RLFSs start towards the
# end of 5'ETS but stretches to 18S then it will counted under 5'ETS. 
# if we would have defined boundary in the first place then we would have lost this junction in our counting. 


entire_RLFSs_rdna<- fread("BK000964_added_3500nt_IGS_upstream_qmrlfs.out.bed", sep = "\t", header = FALSE) #71
entire_RLFSs_rdna<- entire_RLFSs_rdna %>% select(1:6)
entire_RLFSs_rdna<- entire_RLFSs_rdna %>% filter(V2>1299)#68
entire_RLFSs_rdna$v7 <- "junction"

entire_RLFSs_rdna$v8<- entire_RLFSs_rdna$V3-entire_RLFSs_rdna$V2


entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 > 1300 & entire_RLFSs_rdna$V2 < 3500] <- "Promoter"
sum(entire_RLFSs_rdna$v7=="Promoter")
#0

entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 > 3500] <- "Promoter and 5'ETS junction"


entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 > 3501] <- "5'ETS"

entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 > 7507 ] <- "5'ETS and 18S junction"

entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 >= 7508] <- "18S"
entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 > 9377 ] <- "18S and ITS1 junction"


entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 >= 9378   ] <- "ITS1"
entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 > 10377] <- "ITS1 and 5.8S junction"

entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 >= 10378 ] <- "5.8S"
entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 > 10534 ] <- "5.8S and ITS2 junction"

entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 >= 10535] <- "ITS2"
entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 > 11622] <- "ITS2 and 28S junction"

entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 >= 11623 ] <- "28S"
entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 > 16352 ] <- "28S and 3'ETS junction"

entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 >= 16353] <- "3'ETS"
entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 > 16903 ] <- "3'ETS and IGS junction"

entire_RLFSs_rdna$v7[entire_RLFSs_rdna$V2 > 16904 & entire_RLFSs_rdna$V2 < 48806] <- "IGS"


colnames(entire_RLFSs_rdna) <- c("GenBank_Accession", "RLFS_start", "RLFS_end", "RLFS_sequence", "RLFS_name", "strand", "rDNA_region", "RLFS_length")

fwrite(entire_RLFSs_rdna, "RLFS_BK000964_added_3500nt_IGS_upstream_at_junctn_details.csv", sep = ",")

entire_RLFSs_rdna_summary<- entire_RLFSs_rdna %>% group_by(rDNA_region) %>% count()
#rows order is lost 

sum(entire_RLFSs_rdna_summary$n)
#[1] 68

entire_RLFSs_rdna_summary[,1]

#1 28S         3'ETS       5'ETS    5.8S    
#5 IGS          ITS1       ITS2       



new_rows<- data.table(rDNA_region = c("Promoter", "Promoter and 5'ETS junction", "5'ETS and 18S junction", "18S", 
                                      "18S and ITS1 junction","ITS1 and 5.8S junction",
                                       "5.8S and ITS2 junction","ITS2 and 28S junction",
                                      "28S and 3'ETS junction", "3'ETS and IGS junction"))

new_rows$n<- 0


entire_RLFSs_rdna_summary<- rbind(entire_RLFSs_rdna_summary, new_rows)
entire_RLFSs_rdna_summary$rDNA_region
#[1] "28S"                         "3'ETS"                      
#[3] "5'ETS"                       "5.8S"                       
#[5] "IGS"                         "ITS1"                       
#[7] "ITS2"                        "Promoter"                   
#[9] "Promoter and 5'ETS junction" "5'ETS and 18S junction"     
#[11] "18S"                         "18S and ITS1 junction"      
#[13] "ITS1 and 5.8S junction"      "5.8S and ITS2 junction"     
#[15] "ITS2 and 28S junction"       "28S and 3'ETS junction"     
#[17] "3'ETS and IGS junction"      

entire_RLFSs_rdna_summary <- entire_RLFSs_rdna_summary[c(8,9,3,10,11,12,6,13,4,14,7, 15, 1,16,2,17,5),]

names(entire_RLFSs_rdna_summary)[2] <- "RLFS_count"


entire_RLFSs_rdna_summary<- entire_RLFSs_rdna_summary %>% mutate(norm_RLFS_count = RLFS_count/sum(entire_RLFSs_rdna_summary$RLFS_count)) %>% 
  mutate(norm_RLFS_count= round(norm_RLFS_count, 2))

fwrite(entire_RLFSs_rdna_summary, "RLFS_BK000964_added_3500nt_IGS_upstream_at_junctn_graphinput.csv", sep = ",")

entire_RLFSs_rdna_summary$rDNA_region <- factor(entire_RLFSs_rdna_summary$rDNA_region, 
                                                levels = c("Promoter","Promoter and 5'ETS junction", "5'ETS", "5'ETS and 18S junction", 
                                                           "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                           "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                           "28S and 3'ETS junction", "3'ETS", "3'ETS and IGS junction", "IGS" ))


RLFS_norm_3500igs<- ggplot(entire_RLFSs_rdna_summary, aes(x= rDNA_region, y = norm_RLFS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized RLFS distribution in the Mouse rDNA locus", 
       x= "Mouse rDNA region", 
       y= "Normalized RLFS count")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.60))+
  geom_text(aes(label= RLFS_count, vjust= -0.5, size= 50))+
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

ggsave( "Normalized_RLFS_distribution_in_mouse_rDNA_subcomponents_after_rule.tiff", 
        plot = RLFS_norm_3500igs, width=18,height=10, dpi=150)
}


