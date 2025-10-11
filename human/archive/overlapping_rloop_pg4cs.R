#comparing entire rdna region as a whole for count of pG4CS and RLFS 
# next checking how many G4s overlap with RIZ, REZ, LINKER, RLFS 


library(tidyverse)
library(data.table)

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/coorelation_and_length_distribution_human/files")
rlfs<- fread("RLFS_KY962518_added_3500nt_IGS_upstream_master_qmrlfs_table_after_rule.csv", sep = ",", header = TRUE) #195


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/g4s_and_rdna/human/pG4CS_at_rdna_output/files")
pg4cs<- fread("pG4CS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", sep = ",", header = TRUE) #210

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/entire_rdna_rloop_vs_g4s_human/")


rlfs_filt<- rlfs %>% select("#Name", actual_RLFS_start, actual_RLFS_end, location, model, strand, rDNA_region, length_RLFS)
rlfs_filt<- rlfs_filt %>% mutate(RLFS_start = ifelse(rlfs_filt$strand== "+", actual_RLFS_start, actual_RLFS_end))
rlfs_filt<- rlfs_filt %>% mutate(RLFS_end = ifelse(rlfs_filt$strand== "+", actual_RLFS_end, actual_RLFS_start))
rlfs_filt<- rlfs_filt %>% select("#Name", RLFS_start, RLFS_end, location, model, strand, rDNA_region, length_RLFS)



riz_filt<- rlfs %>% select("#Name", start_RIZ, end_RIZ, location, model, strand, rDNA_region, length_RIZ)
rez_filt<- rlfs %>% select("#Name", start_REZ, end_REZ, location, model, strand,rDNA_region, length_REZ)
linker<- rlfs %>% mutate(start_Linker = ifelse(rlfs$strand=="+", end_RIZ, end_REZ))
linker<- linker %>% mutate(end_Linker = ifelse(rlfs$strand=="+", start_REZ, start_RIZ))
linker_filt<- linker %>% select("#Name", start_Linker, end_Linker, location, model, strand,rDNA_region, Linker) %>% filter(Linker>0)


pg4cs_filt<- pg4cs %>% select(GenBank_Accession, pG4CS_start, pG4CS_end, pG4CS_name,pG4CS_sequence, strand,rDNA_region, pG4CS_length)


overlap_summary<- data.frame(matrix(nrow=0, ncol=20))

list_of_data<- list(RIZ = riz_filt, linker=linker_filt, REZ=rez_filt, RLFS= rlfs_filt)
for (i in 1:length(list_of_data)){
  colnames(list_of_data[[i]]) [1]<- "GenBank_Accession" # change the column name first
  
  
  
  pg4cs_filt<- as.data.table(pg4cs_filt)%>% 
    setkey(GenBank_Accession,pG4CS_start,pG4CS_end) #need last two numeric column 
  
  pg4cs_count<- nrow(pg4cs)
  pg4cs_nontemplate<- sum(pg4cs$strand == "+")
  pg4cs_template<-  sum(pg4cs$strand == "-") 
  
  
  expt_data<- as.data.table(list_of_data[[i]]) #experimental data
  expt_count<- nrow(expt_data)
  expt_nontemplate<- sum(expt_data$strand == "+")
  expt_template<-  sum(expt_data$strand == "-") 
  
  
  expt_data<- setkeyv(expt_data, colnames(expt_data)[c(1,2,3)]) #need last two numeric column 
  
  
  
  overlapped<- foverlaps(expt_data, pg4cs_filt, maxgap = 0L, type = "any") %>% 
    na.omit()%>% filter(strand == i.strand) 
  
  #not saved the other strand overlapped data as it is not revelant 
  #notoverlapped<- foverlaps(expt_data, pg4cs_filt, maxgap = 0L, type = "any") %>% 
    #na.omit()%>% filter(strand != i.strand) please note the few of teh identifier also matched with overlapped ones.
  #test<-  overlapped %>% select(GenBank_Accession, pG4CS_start, pG4CS_end,pG4CS_name, pG4CS_sequence, strand, rDNA_region)
  #A<-setdiff(pg4cs_filt,test) #which is same as subtraction of pG4CS overlapped  from pG4CS
  
  overlapped<- overlapped %>% mutate(overlap_length = (pmin(pG4CS_end,overlapped[[10]])-pmax(pG4CS_start,overlapped[[9]]))) #column 9 and 10 are start and end respectively after overlapping. 
  overlapped<- overlapped %>% mutate(pG4CS_overlaping_perc = round((overlap_length/ pG4CS_length)*100))
  overlapped<- overlapped %>% mutate(i.overlaping_perc = round((overlap_length/ overlapped[[15]])*100))
  
  pg4cs_ovlp_count<- length(unique(overlapped$pG4CS_name))
  pg4cs_ovlp_count_perc<-  round((pg4cs_ovlp_count/pg4cs_count)*100) #rounded to nearest integer
  
  pg4cs_ovlp_nontemplate<- length(unique(overlapped$pG4CS_name[overlapped$strand=="+"]))
  pg4cs_ovlp_nontemplate_perc<-  round((pg4cs_ovlp_nontemplate/ pg4cs_nontemplate)*100)
  
  
  pg4cs_ovlp_template<- length(unique(overlapped$pG4CS_name[overlapped$strand=="-"]))
  pg4cs_ovlp_template_perc<-  round((pg4cs_ovlp_template/pg4cs_template)*100)
  
  
  expt_ovlp_count<- length(unique(overlapped$location)) #11 is the unique identifier
  expt_ovlp_count_perc<-  round((expt_ovlp_count/expt_count)*100) 
  
  expt_ovlp_nontemplate<- length(unique(overlapped$location[overlapped$i.strand=="+"]))
  expt_ovlp_nontemplate_perc<-  round((expt_ovlp_nontemplate/ expt_nontemplate)*100)
  
  
  expt_ovlp_template<- length(unique(overlapped$location[overlapped$strand=="-"]))
  expt_ovlp_template_perc<-  round((expt_ovlp_template/expt_template)*100)
  
  a<- "pG4CS"
  b<- names(list_of_data)[i]
  
  new_row<- c(a, pg4cs_count, pg4cs_ovlp_count, pg4cs_ovlp_count_perc, pg4cs_nontemplate, pg4cs_ovlp_nontemplate, pg4cs_ovlp_nontemplate_perc, pg4cs_template,pg4cs_ovlp_template, pg4cs_ovlp_template_perc, 
              b, expt_count,  expt_ovlp_count, expt_ovlp_count_perc, expt_nontemplate, expt_ovlp_nontemplate, expt_ovlp_nontemplate_perc, expt_template,expt_ovlp_template, expt_ovlp_template_perc)
  
  overlap_summary<- rbind(overlap_summary, new_row)
  colnames(overlap_summary)<- c("Dataset", "pg4cs_count", "pg4cs_ovlp_count", "pg4cs_ovlp_count_perc", "pg4cs_nontemplate_count","pg4cs_ovlp_nontemplate_count","pg4cs_ovlp_nontemplate_perc", "pg4cs_template_count", "pg4cs_ovlp_template_count","pg4cs_ovlp_template_perc", 
                                "Experimental_Dataset","expt_count",  "expt_ovlp_count", "expt_ovlp_count_perc", "expt_nontemplate_count", "expt_ovlp_nontemplate_count", "expt_ovlp_nontemplate_perc", "expt_template_count", "expt_ovlp_template_count","expt_ovlp_template_perc")
  
  fwrite(overlapped, paste("pG4CS_overlapped_with_", names(list_of_data)[i], ".csv", sep = ""), sep = ",")
  
  
  
  #fwrite(notoverlapped, paste("pG4CS_notoverlapped_with_", names(list_of_data)[i], ".csv", sep = ""), sep = ",")
  #also note that overlapped %>% group_by(i.rDNA_region) %>% count() or overlapped %>% group_by(rDNA_region) %>% count() is same.
}

fwrite(overlap_summary, "pG4CS_overlapping_summary_with_RIZ_REZ_LINKER_RLFS.csv")
 

#visualization


# Plot with text labels
overlap_summary$Experimental_Dataset <- factor(overlap_summary$Experimental_Dataset, 
                               levels = c("RIZ", "linker", "REZ", "RLFS"))



p1<- ggplot(overlap_summary, aes(x = Experimental_Dataset, y = pg4cs_ovlp_count_perc, fill = Experimental_Dataset)) +
  geom_bar(stat = "identity", fill="cornflowerblue") + 
  geom_text(aes(label = pg4cs_ovlp_count), vjust = -0.5, size = 5) +  # Adding percentage labels above bars
  labs(title = "pG4CS Overlap with Different Datasets", 
       x = "Dataset", 
       y = "Number of Overlapping pG4CS") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        #panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))

ggsave( "entire_pG4CS_overlap_with_RLFS_and_its_components.tiff", 
        plot = p1, width=15,height=10, dpi=150)



p2<- ggplot(overlap_summary, aes(x = Experimental_Dataset, y = expt_ovlp_count_perc, fill = Experimental_Dataset)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = expt_ovlp_count), vjust = -0.5, size = 5) +  
  labs(title = "RLFS and its components Dataset Overlap", 
       x = "Dataset", 
       y = "Percentage of Overlapping Elements") +
  scale_fill_manual(values= c("RLFS"= "darkgrey","RIZ"= "pink", "linker"= "cyan","REZ"= "steelblue"))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank())
  

ggsave( "entire_RLFS_and_its_components_ovlp_pG4CS.tiff", 
        plot = p2, width=15,height=10, dpi=150)



#make for nontemplate
nt_p1<- ggplot(overlap_summary, aes(x = Experimental_Dataset, y = pg4cs_ovlp_nontemplate_perc, fill = Experimental_Dataset)) +
  geom_bar(stat = "identity", fill="cornflowerblue") + 
  geom_text(aes(label = pg4cs_ovlp_nontemplate_count), vjust = -0.5, size = 5) +  # Adding percentage labels above bars
  labs(title = "Non-template pG4CS Overlap with Different Datasets", 
       x = "Dataset", 
       y = "Number of Overlapping pG4CS") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        #panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))

ggsave( "Nontemplate_entire_pG4CS_overlap_with_RLFS_and_its_components.tiff", 
        plot = nt_p1, width=15,height=10, dpi=150)



nt_p2<- ggplot(overlap_summary, aes(x = Experimental_Dataset, y = expt_ovlp_nontemplate_perc, fill = Experimental_Dataset)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = expt_ovlp_nontemplate_count), vjust = -0.5, size = 5) +  
  labs(title = "Non-template RLFS and its components Dataset Overlap", 
       x = "Dataset", 
       y = "Percentage of Overlapping Elements") +
  scale_fill_manual(values= c("RLFS"= "darkgrey","RIZ"= "pink", "linker"= "cyan","REZ"= "steelblue"))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank())


ggsave( "Nontemplate_entire_RLFS_and_its_components_ovlp_pG4CS.tiff", 
        plot = nt_p2, width=15,height=10, dpi=150)


  
#make for template
t_p1<- ggplot(overlap_summary, aes(x = Experimental_Dataset, y = pg4cs_ovlp_template_perc, fill = Experimental_Dataset)) +
  geom_bar(stat = "identity", fill="cornflowerblue") + 
  geom_text(aes(label = pg4cs_ovlp_template_count), vjust = -0.5, size = 5) +  # Adding percentage labels above bars
  labs(title = "Template pG4CS Overlap with Different Datasets", 
       x = "Dataset", 
       y = "Number of Overlapping pG4CS") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        #panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))

ggsave( "Template_entire_pG4CS_overlap_with_RLFS_and_its_components.tiff", 
        plot = t_p1, width=15,height=10, dpi=150)



t_p2<- ggplot(overlap_summary, aes(x = Experimental_Dataset, y = expt_ovlp_template_perc, fill = Experimental_Dataset)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = expt_ovlp_template_count), vjust = -0.5, size = 5) +  
  labs(title = "Template RLFS and its components Dataset Overlap", 
       x = "Dataset", 
       y = "Percentage of Overlapping Elements") +
  scale_fill_manual(values= c("RLFS"= "darkgrey","RIZ"= "pink", "linker"= "cyan","REZ"= "steelblue"))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank())


ggsave( "Template_entire_RLFS_and_its_components_ovlp_pG4CS.tiff", 
        plot = t_p2, width=15,height=10, dpi=150)




