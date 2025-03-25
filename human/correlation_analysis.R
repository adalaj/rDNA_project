#Task is to find coorelation between alot of parameters in human 

#1) Correlation between average length of RLFS,RIZ or pG4CS and same length of GC skew strandwise?
#2) Correlation between RLFS and pG4CS in human rDNA locus (strandwise)
#3) Correlation between RIZ (m1 and m2) and pG4CS in human rDNA locus (strandwise)


# to start with first #1) Correlation between average length of RLFS,RIZ or pG4CS and same length of GC skew?
# i need to first get table format of RLFS to get information where RIZ is starting and ending. 
# find average length of RIZ, RLFS and pG4CS in each subsection of rDNA (as there is nothing found initiating at junction). 

#open terminal # conda activate python2.7>cd Downloads 
#(python2.7) python QmRLFS-finder.py -i KY962518_added_3500nt_IGS_upstream_nontemplate.fasta -o KY962518_added_3500nt_IGS_upstream_nontemplate_qmrlfs
# deafault table format is saved. 
#Time used: 0.72 mins


library(data.table)
library(tidyverse)

rdna_rlfs_table<- fread("KY962518_added_3500nt_IGS_upstream_nontemplate_qmrlfs.out.table.txt", sep="\t", header= TRUE) #208 
# this some part of IGS at start and some part of promoter at end.


# interesting observation of table format over bed format
# the column location in table format is three component: identifier:numberX-numberY
# the numberX-numberY is the  case when strand is plus (+) tell where RLFS starts+1 bp and RLFS end position respectively. 
# However, in case when strand is reversed then numberX is where RLFS end + 1 and numberY becomes RLFS start.

#so my approach to filter startRIZ with rdna after rule boundary is problematic. 
# because in some case rdna_rlfs$rdna_region[rdna_rlfs$start_RIZ > 7157 & rdna_rlfs$start_RIZ <= 9026] <- "18S"
# it will detect 3 RLFS however in reality there is none. 


#the approach would be to calculate RLFS start and RLFS end considering strand then doing further analysis

#for plus strand which is non-template in this case
# actual start of RLFS = start_RIZ
# actual end of RLFS = end_REZ

#for minus strand which is template in this case
# actual start of RLFS = end_RIZ 
# actual end of RLFS = start_REZ

rdna_rlfs_table<- rdna_rlfs_table %>% mutate(actual_RLFS_start = ifelse(rdna_rlfs_table$strand== "+", start_RIZ, end_RIZ))
rdna_rlfs_table<- rdna_rlfs_table %>% mutate(actual_RLFS_end = ifelse(rdna_rlfs_table$strand== "+", end_REZ, start_REZ))
rdna_rlfs_table$length_RLFS <- abs(rdna_rlfs_table$actual_RLFS_start-rdna_rlfs_table$actual_RLFS_end)


#also create another column with length of this RLFS.
# same this as master table

rdna_rlfs_table_filt<- rdna_rlfs_table[rdna_rlfs_table$actual_RLFS_start >1299 & rdna_rlfs_table$actual_RLFS_start < 46137] #195 
rdna_rlfs_table_filt$rDNA_region <- "junction"



#strand specificity will matter here now we want to allocate RLFS to rdna region sub components  based on their direction
# to do that, i will be creating a column that will have actual RLFS start based on strand specificity

rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start > 1299 & rdna_rlfs_table_filt$actual_RLFS_start < 3500] <- "Promoter"
sum(rdna_rlfs_table_filt$rDNA_region=="Promoter")#13

rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start > 3500] <- "Promoter and 5'ETS junction"

rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start >= 3501] <- "5'ETS"
rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start > 7157] <- "5'ETS and 18S junction"

rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start >= 7158] <- "18S"
rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start > 9026] <- "18S and ITS1 junction"



rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start >= 9027] <- "ITS1"
rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start > 10096] <- "ITS1 and 5.8S junction"

rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start >= 10097] <- "5.8S"
rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start > 10253] <- "5.8S and ITS2 junction"

rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start >= 10254] <- "ITS2"
rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start > 11420] <- "ITS2 and 28S junction"

rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start >= 11421] <- "28S"
rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start > 16471] <- "28S and 3'ETS junction"

rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start >= 16472] <- "3'ETS"
rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start > 16832] <- "3'ETS and IGS junction"

rdna_rlfs_table_filt$rDNA_region[rdna_rlfs_table_filt$actual_RLFS_start > 16833 & rdna_rlfs_table_filt$actual_RLFS_start < 46137] <- "IGS"

fwrite(rdna_rlfs_table_filt, "RLFS_KY962518_added_3500nt_IGS_upstream_master_table_after_rule.csv", sep = ",")



# I want to plot average length of RIZ, RLFS and G4s in rdna region 
# avearge length of RIZ,RLFS and G4s will be compared again similar sliding window of GC skew

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/g4s_and_rdna/human/pG4CS_at_rdna_output/files")


entire_g4s_rdna <- fread("pG4CS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", header = TRUE, sep = ",") #210

avg_g4cs<- round(mean(entire_g4s_rdna$pG4CS_length),2) #23.85
avg_rlfs<- round(mean(rdna_rlfs_table_filt$length_RLFS),2) #842.89
avg_riz <- round(mean(rdna_rlfs_table_filt$length_RIZ),2) #17.03
avg_rez<- round(mean(rdna_rlfs_table_filt$length_REZ),2) #817.38

length_distribution <- data.frame(matrix(nrow = 0, ncol = 6))
colnames(length_distribution)<- c("rDNA_region", "length_RIZ", "Linker", "length_REZ","length_RLFS", "strand")

for ( i in unique(rdna_rlfs_table_filt$rDNA_region)){
  rdna<- rdna_rlfs_table_filt %>% 
                filter (rdna_rlfs_table_filt$rDNA_region == i) %>% 
                select(rDNA_region, length_RIZ, Linker, length_REZ,length_RLFS, strand)
  length_distribution <- rbind(length_distribution, rdna)
}

template_length<- length_distribution %>% filter(strand == "-")
nontemplate_length <- length_distribution %>% filter(strand == "+")

list_of_data<- list(length_distribution, template_length, nontemplate_length)
names(list_of_data)<- c("entire_rDNA_RLFS", "template_rDNA_RLFS", "nontemplate_rDNA_RLFS")


for ( i in 1: length(list_of_data)){
  length_distribution_graph<- data.frame(matrix(nrow = 0, ncol=5))
  colnames(length_distribution_graph)<- c("rDNA_region", "Mean_RLFS_length", "Mean_RIZ_length", "Mean_linker_length", "Mean_REZ_length")
  
   data_of_interest <- as.data.frame(list_of_data[[i]])
  for (j in unique(data_of_interest$rDNA_region)) {
    tmp<- data_of_interest %>% filter(data_of_interest$rDNA_region == j)
    region_of_interest<- j
    avg_rlfs<- round(mean(tmp$length_RLFS),2)  #842.89
    avg_riz <- round(mean(tmp$length_RIZ),2) #17.03
    avg_linker <- round(mean(tmp$Linker),2) #17.03
    avg_rez<- round(mean(tmp$length_REZ),2) #817.38

    
    length_distribution_graph<- rbind(length_distribution_graph, 
                                      data.frame(rDNA_region = region_of_interest,
                                                 Mean_RLFS_length = avg_rlfs,
                                                 Mean_RIZ_length = avg_riz,
                                                 Mean_linker_length = avg_linker,
                                                 Mean_REZ_length = avg_rez,
                                                 stringsAsFactors = FALSE))
    
  }
  fwrite(length_distribution_graph, paste("length_distribution_graph_input_", names(list_of_data)[i], ".csv", sep = ""))
}
      
  
length_distribution_entire_rdna<- fread("length_distribution_graph_input_entire_rDNA_RLFS.csv", header = TRUE, sep = ",")


data_long <- length_distribution_entire_rdna %>%
  pivot_longer(cols = -rDNA_region, 
               names_to = "Property", 
               values_to = "Length")


ggplot(data_long, aes(x= rDNA_region, y = Length, fill= Property)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized RLFS strandwise distribution in the Human rDNA locus", 
       x= "Human rDNA region", 
       y= "Normalized RLFS count", 
       fill= "RLFS strand")+
  #scale_y_continuous(breaks= seq(0, 0.30, by = 0.1), limits =c(0,0.30))+
  geom_text(aes(label= Length), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  scale_fill_manual(values= c("#F5FEFB", "#E21515", "#5AAA46", "#F36017", "#6B1519", 
                                  "#818689", "#ECE612", "#E07F80", "#DE9A22"))+
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
        axis.ticks.y = element_line(color = "black"))
  
  
  
  
  
  
  
  












