#Task is to find coorelation between alot of parameters in human 

# There are 8 parameters which I am considering for correlation heatmap
# length of RIZ, LINKER, REZ, RLFS, G4S, count of RIZ, counts G4s, GC skew
# there will be three heatmaps total, template and non template



# to start with first #1) Calculation of average length of RLFS,RIZ,rex, linker or pG4CS?
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

fwrite(rdna_rlfs_table_filt, "RLFS_KY962518_added_3500nt_IGS_upstream_master_qmrlfs_table_after_rule.csv", sep = ",")



# I want to plot average length of RIZ, RLFS and G4s in rdna region 
# average length of RIZ,RLFS and G4s will be compared again similar sliding window of GC skew


#2) comparing pG4CS count, RIZ count
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/output/files")

riz_count_entire <- fread ("RLFS_KY962518_added_3500nt_IGS_upstream_at_junctn_after_rule_graphinput.csv", sep = ",", header = TRUE)
riz_count_entire<- riz_count_entire[!grepl("junction", riz_count_entire$rDNA_region),] #nrow=9
riz_count_entire$non_canonical_str <- "RIZ"
colnames(riz_count_entire)<- c("rDNA_region", "count", "norm_count", "non_canonical_str")


riz_count_strandwise<- fread("RLFS_KY962518_added_3500nt_IGS_upstream_no_junctn_strandwise_AR_graphinput.csv", sep = ",", header = TRUE)
riz_count_strandwise$non_canonical_str<- "RIZ"
colnames(riz_count_strandwise)<- c("rDNA_region", "strand", "count", "norm_count", "non_canonical_str")

riz_nontemplate<-riz_count_strandwise %>% filter(strand == "+") #nrow=9
riz_template<- riz_count_strandwise %>% filter(strand == "-") #nrow=9


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/g4s_and_rdna/human/pG4CS_at_rdna_output/files")

g4s_count<- fread("pG4CS_KY962518_added_3500nt_IGS_upstream_at_junctn_graphinput.csv", sep = ",", header = TRUE)
g4s_count<- g4s_count[!grepl("junction", g4s_count$rDNA_region),] #nrow=9
g4s_count$non_canonical_str<- "pG4CS"
colnames(g4s_count)<- c("rDNA_region", "count", "norm_count", "non_canonical_str")

g4s_count_strandwise <- fread("pG4CS_KY962518_added_3500nt_IGS_upstream_no_junctn_strandwise_graphinput.csv", sep = ",", header = TRUE)
g4s_count_strandwise$non_canonical_str<- "pG4CS"
colnames(g4s_count_strandwise)<- c("rDNA_region","strand","count", "norm_count", "non_canonical_str")


g4s_nontemplate<-g4s_count_strandwise %>% filter(strand == "+") #nrow=9
g4s_template<- g4s_count_strandwise %>% filter(strand == "-") #nrow=9



entire_riz_vs_g4s<- rbind(riz_count_entire, g4s_count) #nrow=18
nontemplate_riz_vs_g4s<- rbind(riz_nontemplate,g4s_nontemplate) #nrow=18
template_riz_vs_g4s<- rbind(riz_template,g4s_template) #nrow=18



setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/coorelation_and_length_distribution_human/files")

count_list<- list(entire_riz_vs_g4s, template_riz_vs_g4s, nontemplate_riz_vs_g4s)
names(count_list)<- c("entire_RIZ_vs_pG4CS", "template_RIZ_vs_pG4CS", "nontemplate_RIZ_vs_pG4CS")

for (i in 1:length(count_list)){
  print(paste("count_distribution_graphinput_", names(count_list)[i], ".csv", sep = ""))
fwrite(count_list[[i]], paste("count_distribution_graphinput_", names(count_list)[i], ".csv", sep = ""), sep = ",")
}



#to plot graphs simultaneously:

for(i in names(count_list)){
  filename<- paste("count_distribution_graphinput_", i, ".csv", sep = "")
  graph_data<- fread(filename, header = TRUE, sep = ",")
  print(paste(i, "_count_distribution_in_human_rdna.tiff", sep = ""))
  
  graph_data$rDNA_region <- factor(graph_data$rDNA_region, 
                                   levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                              "ITS2","28S", "3'ETS", "IGS" ))
  graph_data$non_canonical_str<- factor(graph_data$non_canonical_str, 
                                             levels = c("RIZ", "pG4CS"))
  
  count_graph<- ggplot(graph_data, aes(x= rDNA_region, y = norm_count, fill= non_canonical_str)) + 
    geom_bar(stat= "identity", position ="dodge", color = "black") +
    labs(title= paste(i, "_count_distribution", sep = ""), 
         x= "Human rDNA region", 
         y= "Normalized count", 
         fill= "Non Canonical Structures")+
    scale_y_continuous(breaks= seq(0, 0.30, by = 0.1), limits =c(0,0.40))+
    geom_text(aes(label= count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
    scale_fill_manual(values= c("RIZ" = "pink", "pG4CS" = "cornflowerblue"))+
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
  
  ggsave(paste(i, "_count_distribution_in_human_rdna.tiff", sep = ""), 
         plot = count_graph, width = 18, height = 10, dpi = 150)
  
           }




# 3) comparing length of RLFS with RIZ, Linker and REZ
length_distribution <- data.frame(matrix(nrow = 0, ncol = 6))
colnames(length_distribution)<- c("rDNA_region", "length_RIZ", "Linker", "length_REZ","length_RLFS", "strand")

for ( i in unique(rdna_rlfs_table_filt$rDNA_region)){
  rdna<- rdna_rlfs_table_filt %>% 
                filter (rdna_rlfs_table_filt$rDNA_region == i) %>% 
                select(rDNA_region, length_RIZ, Linker, length_REZ,length_RLFS, strand)
  length_distribution <- rbind(length_distribution, rdna)
  
}

#as we also want to see avlues in 18s and 5.8 we will create new rows

new_rows <- data.frame(rDNA_region= c("18S", "5.8S", "5.8S"),
                       length_RIZ=c(0, 0, 0),
                       Linker= c(0, 0, 0),
                       length_REZ = c(0, 0, 0),
                       length_RLFS = c(0, 0, 0),
                       strand = c("+", "+", "-"))
  
  
length_distribution <- rbind(length_distribution, new_rows)
fwrite(length_distribution, "RLFS_Length_distribution_entire_human_rDNA_master.csv", sep = ",")

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
  fwrite(length_distribution_graph, paste("length_distribution_", names(list_of_data)[i], ".csv", sep = ""))
  data_long <- length_distribution_graph %>%
    pivot_longer(cols = -rDNA_region, 
                 names_to = "Property", 
                 values_to = "Length")
  
  data_long<- data_long %>% mutate(round_length = round(data_long$Length))
  fwrite(data_long, paste("length_distribution_data_long_graphinput_", names(list_of_data)[i], ".csv", sep = ""))
}
      
  
#to plot graphs simultaneously:

for(i in names(list_of_data)){
  filename<- paste("length_distribution_data_long_graphinput_", i, ".csv", sep = "")
  graph_data<- fread(filename, header = TRUE, sep = ",")
  graph_data$Property <- factor(graph_data$Property, 
                                           levels = c("Mean_RLFS_length", 
                                                      "Mean_RIZ_length", 
                                                      "Mean_REZ_length",
                                                      "Mean_linker_length"))
  
  graph_data$rDNA_region <- factor(graph_data$rDNA_region, 
                                              levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                         "ITS2","28S", "3'ETS", "IGS" ))
  
  graph<- ggplot(graph_data, aes(x= rDNA_region, y = round_length, fill= Property)) + 
    geom_bar(stat= "identity", position ="dodge", color = "black") +
    labs(title= paste(i,"_Length_distribution", sep = ""), 
         x= "Human rDNA region", 
         y= "Length (Rounded to Nearest Integer)", 
         fill= "Parameter")+
    scale_y_continuous(breaks= seq(0, 2100, by = 500), limits =c(0,2100))+
    geom_text(aes(label= round_length), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
    scale_fill_manual(values= c("Mean_RLFS_length"= "darkgrey","Mean_RIZ_length"= "pink", "Mean_linker_length"= "cyan","Mean_REZ_length"= "steelblue"))+
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
  
  ggsave(paste(i, "_length_distribution_in_human_rdna.tiff", sep = ""), 
         plot = graph, width = 18, height = 10, dpi = 150)
  
}


#4) Next, want to make similar length distribution comparison for G4s


entire_g4s_rdna <- fread("pG4CS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", header = TRUE, sep = ",") #210

g4s_length_distribution <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(g4s_length_distribution)<- c("rDNA_region", "length_pG4CS", "strand")

for ( i in unique(entire_g4s_rdna$rDNA_region)){
  rdna<- entire_g4s_rdna %>% 
    filter (entire_g4s_rdna$rDNA_region == i) %>% 
    select(rDNA_region, length_pG4CS, strand)
  g4s_length_distribution <- rbind(g4s_length_distribution, rdna)
  
}

new_rows<- data.table(rDNA_region = c( "Promoter", "ITS1", "18S","18S", "5.8S","5.8S"),
                      length_pG4CS = c(0, 0, 0,0,0,0),
                      strand= c("+", "-", "+", "-", "+", "-"))

g4s_length_distribution <- rbind(g4s_length_distribution, new_rows)
fwrite(g4s_length_distribution, "pG4CS_length_distribution_entire_human_rDNA_master.csv", sep = ",")

g4s_template_length<- g4s_length_distribution %>% filter(strand == "-")
g4s_nontemplate_length <- g4s_length_distribution %>% filter(strand == "+")

list_of_data<- list(g4s_length_distribution, g4s_template_length, g4s_nontemplate_length)
names(list_of_data)<- c("entire_rDNA_pG4CS", "template_rDNA_pG4CS", "nontemplate_rDNA_pG4CS")


for ( i in 1: length(list_of_data)){
  g4s_length_distribution_graph<- data.frame(matrix(nrow = 0, ncol=5))
  colnames(length_distribution_graph)<- c("rDNA_region", "Mean_pG4CS_length")
  
  data_of_interest <- as.data.frame(list_of_data[[i]])
  for (j in unique(data_of_interest$rDNA_region)) {
    tmp<- data_of_interest %>% filter(data_of_interest$rDNA_region == j)
    region_of_interest<- j
    avg_pG4CS<- round(mean(tmp$length_pG4CS),2)
    
    
    g4s_length_distribution_graph<- rbind(g4s_length_distribution_graph, 
                                      data.frame(rDNA_region = region_of_interest,
                                                 Mean_pG4CS_length = avg_pG4CS,
                                                 stringsAsFactors = FALSE))
    
  }
  
  
  fwrite(g4s_length_distribution_graph, paste("length_distribution_", names(list_of_data)[i], ".csv", sep = ""))
  data_long <- g4s_length_distribution_graph %>%
    pivot_longer(cols = -rDNA_region, 
                 names_to = "Property", 
                 values_to = "Length")
  
  data_long<- data_long %>% mutate(round_length = round(data_long$Length))
  fwrite(data_long, paste("length_distribution_data_long_graphinput_", names(list_of_data)[i], ".csv", sep = ""))
}


#to plot graphs simultaneously:

for(i in names(list_of_data)){
  filename<- paste("length_distribution_data_long_graphinput_", i, ".csv", sep = "")
  graph_data<- fread(filename, header = TRUE, sep = ",")
  
  graph_data$rDNA_region <- factor(graph_data$rDNA_region, 
                                   levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                              "ITS2","28S", "3'ETS", "IGS" ))
  
  graph<- ggplot(graph_data, aes(x= rDNA_region, y = round_length, fill= Property)) + 
    geom_bar(stat= "identity", position ="dodge", color = "black") +
    labs(title= paste(i,"_length_distribution"), 
         x= "Human rDNA region", 
         y= "Length (Rounded to Nearest Integer)", 
         fill= "Parameter")+
    scale_y_continuous(breaks= seq(0, 30, by = 10), limits =c(0,30))+
    geom_text(aes(label= round_length), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
    scale_fill_manual(values= c("Mean_pG4CS_length"= "cornflowerblue"))+
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
  
  ggsave(paste(i, "_length_distribution_in_human_rdna.tiff", sep = ""), 
         plot = graph, width = 18, height = 10, dpi = 150)
  
}



#5) Next, plot length of RIZ and pG4CS next to each other 

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/coorelation_and_length_distribution_human/files")

entire_rlfs<- fread("length_distribution_entire_rDNA_RLFS.csv", sep = ",", header = TRUE)
entire_riz<- entire_rlfs %>% select(rDNA_region, Mean_RIZ_length)
colnames(entire_rlfs)<- c("rDNA_region", "length")
entire_rlfs$non_canonical_str <- "RIZ"


entire_g4s<- fread("length_distribution_entire_rDNA_pG4CS.csv", sep = ",", header = TRUE)
colnames(entire_g4s)<- c("rDNA_region", "length")
entire_g4s$non_canonical_str<- "pG4CS"

entire_rdna_riz_g4s<- rbind(entire_riz, entire_g4s)




nontemplate_riz<- fread("length_distribution_nontemplate_rDNA_RLFS.csv", sep = ",", header = TRUE)
nontemplate_riz<- nontemplate_riz %>% select(rDNA_region, Mean_RIZ_length)
colnames(nontemplate_riz)<- c("rDNA_region", "length")
nontemplate_riz$non_canonical_str <- "RIZ"

nontemplate_g4s<- fread("length_distribution_nontemplate_rDNA_pG4CS.csv", sep = ",", header = TRUE)
colnames(nontemplate_g4s)<- c("rDNA_region", "length")
nontemplate_g4s$non_canonical_str<- "pG4CS"

nontemplate_rdna_riz_g4s<- rbind(nontemplate_riz, nontemplate_g4s)


template_riz<- fread("length_distribution_template_rDNA_RLFS.csv", sep = ",", header = TRUE)
template_riz<- template_riz %>% select(rDNA_region, Mean_RIZ_length)
colnames(template_riz)<- c("rDNA_region", "length")
template_riz$non_canonical_str <- "RIZ"

template_g4s<- fread("length_distribution_template_rDNA_pG4CS.csv", sep = ",", header = TRUE)
colnames(template_g4s)<- c("rDNA_region", "length")
template_g4s$non_canonical_str<- "pG4CS"

template_rdna_riz_g4s<- rbind(template_riz, template_g4s)


data_list<- list(entire_rdna_riz_g4s,nontemplate_rdna_riz_g4s,template_rdna_riz_g4s)
names(data_list)<- c("entire_rdna_riz_g4s","nontemplate_rdna_riz_g4s","template_rdna_riz_g4s")

for (i in 1:length(data_list)){
  
  data <- as.data.frame(data_list[[i]])
  data$rDNA_region <- factor(data$rDNA_region, 
                                   levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                              "ITS2","28S", "3'ETS", "IGS" ))
  data$non_canonical_str <- factor(data$non_canonical_str, 
                             levels = c("RIZ", "pG4CS"))
  
  data<- data %>% mutate(round_length = round(data$length))
  fwrite(data, paste("length_distribution_",names(data_list)[i], "_graphinput.csv", sep = ""), sep = ",")
  
  
  graph_riz_g4<- ggplot(data, aes(x= rDNA_region, y = length, fill= non_canonical_str)) + 
    geom_bar(stat= "identity", position ="dodge", color = "black") +
    labs(title= paste(names(data_list)[i], "_length_distribution", sep = ""), 
         x= "Human rDNA region", 
         y= "Length (Rounded to Nearest Integer)", 
         fill= "Non Canonical structure")+
    scale_y_continuous(breaks= seq(0, 30, by = 10), limits =c(0,30))+
    geom_text(aes(label= round_length), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
    scale_fill_manual(values= c("RIZ" = "pink", "pG4CS" = "cornflowerblue"))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
          axis.ticks.y = element_line(color = "black"))
  
  ggsave(paste(names(data_list)[i], "_length_distribution_in_human_rdna.tiff", sep = ""), 
         plot = graph_riz_g4, width = 18, height = 10, dpi = 150)
}






#6) Make a correlation matrix between 


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/g4s_and_rdna/human/pG4CS_at_rdna_output/files")
g4s_count<- fread("pG4CS_KY962518_added_3500nt_IGS_upstream_at_junctn_graphinput.csv", sep = ",", header = TRUE)
g4s_count<- g4s_count[!grepl("junction", g4s_count$rDNA_region),] #nrow=9

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/output/files")
riz_count_entire <- fread ("RLFS_KY962518_added_3500nt_IGS_upstream_at_junctn_after_rule_graphinput.csv", sep = ",", header = TRUE)
riz_count_entire<- riz_count_entire[!grepl("junction", riz_count_entire$rDNA_region),] #nrow=9

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/coorelation_and_length_distribution_human/files")

rlfs_length<- fread("length_distribution_entire_rDNA_RLFS.csv")
pG4CS_length<- fread("length_distribution_entire_rDNA_pG4CS.csv")


correlation_entire<- rlfs_length %>% left_join(pG4CS_length, by=join_by(rDNA_region)) #only takes two dataset at a time
correlation_entire<- correlation_entire %>% left_join(g4s_count, by=join_by(rDNA_region))
correlation_entire<- correlation_entire %>% left_join(riz_count_entire, by=join_by(rDNA_region))

#factor/levels doesnt matter here

correlation_entire_input<- correlation_entire %>% select(-rDNA_region,-pG4CS_count, -RLFS_count)
fwrite(correlation_entire_input, "entire_rdna_inputfile_for_correlation_matrix.csv", sep = ",")

install.packages("corrplot")
library(corrplot)

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/coorelation_and_length_distribution_human/files")

cor_matrix= cor(correlation_entire_input, method = c("spearman"))

fwrite(cor_matrix, "spearman_correlation_matrix_graphinput_entire_human_rdna.csv", sep=",")

png("correlation_plot_entire_human_rdna.png", width = 1600, height = 1200, res=150)
corrplot(cor_matrix, type = "lower", method = "square", 
                    addCoef.col = "black",  # Color of correlation coefficient labels
                    col = colorRampPalette(c("blue", "white", "red"))(200),  # Custom color palette
                    tl.col = "black",       # Text label color
                    tl.cex = 1.2,           # Increase text label size
                    number.cex = 1.2,       # Increase correlation coefficient size
                    tl.srt = 45,            # Rotate axis labels
                    title = "Correlation Matrix of Parameters in entire human rDNA (not strandwise)", 
                    mar=c(0,0,1,0))
         # Add title
                         # Adjust margins to accommodate title

dev.off()


#prepare for non-template data and template
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/g4s_and_rdna/human/pG4CS_at_rdna_output/files")

g4s_count_strandwise<- fread("pG4CS_KY962518_added_3500nt_IGS_upstream_no_junctn_strandwise_graphinput.csv")
g4s_count_nontemplate<- g4s_count_strandwise %>% filter(strand == "+")
g4s_count_template<- g4s_count_strandwise %>% filter(strand == "-")

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/output/files")
riz_count_strandwise<- fread("RLFS_KY962518_added_3500nt_IGS_upstream_no_junctn_strandwise_AR_graphinput.csv")
riz_count_nontemplate<- riz_count_strandwise %>% filter(strand == "+")
riz_count_template<- riz_count_strandwise %>% filter(strand == "-")


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/coorelation_and_length_distribution_human/files")

rlfs_nontemplate<- fread("length_distribution_nontemplate_rDNA_RLFS.csv", header = TRUE, sep = ",")
rlfs_template<- fread("length_distribution_template_rDNA_RLFS.csv", header = TRUE, sep = ",")

pG4CS_nontemplate<- fread("length_distribution_nontemplate_rDNA_pG4CS.csv", header = TRUE, sep = ",")
pG4CS_template<- fread("length_distribution_template_rDNA_pG4CS.csv", header = TRUE, sep = ",")


correlation_filenames <- list(nontemplate = c("nontemplate"), 
                              template = c("template"))

for ( i in names(correlation_filenames)){
  filename<- paste("correlation_",i, sep = "")
  rlfs_length<- get(paste("rlfs_", i, sep = ""))
  pG4CS_length<- get(paste("pG4CS_", i, sep = ""))
  g4s_count<- get(paste("g4s_count_", i, sep = ""))
  riz_count<- get(paste("riz_count_", i, sep = ""))

  filename<- rlfs_length %>% left_join(pG4CS_length, by=join_by(rDNA_region))

#only takes two dataset at a time
  filename<- filename %>% left_join(g4s_count, by=join_by(rDNA_region))
  filename<- filename %>% left_join(riz_count, by=join_by(rDNA_region))

  filename2<- paste(filename, "_input", sep = "")
  filename2<- filename %>% select(-rDNA_region,-pG4CS_count, -RLFS_count,-strand.x, -strand.y) #bcoz it only recognise numbers, no character needed
fwrite(filename2, paste(i, "_inputfile_for_correlation_matrix.csv", sep= ""), sep = ",")



cor_matrix= cor(filename2, method = c("spearman")) #cor needs matrix meaning the row names as to be changed otherwise it will print as 1,2,3 on vertical side
png(paste("correlation_plot_", i, "human_rdna.png",sep = ""), width = 1600, height = 1200, res=150)
corrplot(cor_matrix, type = "lower", method = "square", 
         addCoef.col = "black",  # Color of correlation coefficient labels
         col = colorRampPalette(c("blue", "white", "red"))(200),  # Custom color palette
         tl.col = "black",       # Text label color
         tl.cex = 1.2,           # Increase text label size
         number.cex = 1.2,       # Increase correlation coefficient size
         tl.srt = 45,            # Rotate axis labels
         title = paste("Correlation Matrix of Parameters in", i, "human rDNA",sep=" "), 
         mar=c(0,0,1,0))
# Add title
# Adjust margins to accommodate title

dev.off()

fwrite(cor_matrix, paste("spearman_correlation_matrix_graphinput_",i, ".csv", sep=""), sep=",")
}



#make similar for strand specific 




# I wanted to try this but it didnt work! 

{
  #Mainly becuase on X-axis it is counting RIZ and on Y-axis it counting G4. 
  #and each dot represents the bins. 
  #Each dot in your scatter plot represents one genomic bin (e.g., 30 bp or 100 bp region, depending on your bin size).
  
  #For example:
    
  #If you see four dots at (0, y-value) on the x-axis (G4s count = 0), it means that four different bins have zero RIZs but different G4s counts on the y-axis.
  #This tells you that some bins have G4s without RIZ.
  
  #Similarly:
  #If a dot appears at (x, 0) on the y-axis, it means that there are bins with RIZ but no G4s.
  #If a dot is at (0,0), it means there are bins with neither G4s nor RIZ
  
  #If many dots are clustered around (0,0), most bins have low or no G4s and RIZ.
  #If dots are spread diagonally, there is a positive correlation (bins with high G4s also tend to have high RIZ).
  #If dots are all over with no pattern, there’s no strong correlation between G4s and RIZ.
  
  
#On a scatter diagram, the closer the points lie to a straight line, the stronger the linear relationship between two variables. 
#To quantify the strength of the relationship, we can calculate the correlation coefficient.
#A value of the correlation coefficient close to +1 indicates a strong positive linear relationship (i.e. one variable increases with the other). 
#A value close to -1 indicates a strong negative linear relationship (i.e. one variable decreases as the other increases
#A value close to 0 indicates no linear relationship; however, there could be a nonlinear relationship (a hyperbola type) between the variables
# refer: https://pmc.ncbi.nlm.nih.gov/articles/PMC374386/

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/g4s_and_rdna/human/pG4CS_at_rdna_output/files")
entire_g4s_rdna <- fread("pG4CS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", header = TRUE, sep = ",") #210

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/coorelation_and_length_distribution_human/files")
rlfs<- fread("RLFS_KY962518_added_3500nt_IGS_upstream_master_qmrlfs_table_after_rule.csv", header = TRUE, sep = ",") #195

#as total rdna length is 44838
range(entire_g4s_rdna$actual_pG4CS_start)
#[1]  2197 44857

range(rlfs$actual_RLFS_start)
#[1]  2267 43029

bin_size <- c(30)#, 50, 75, 100, 200, 500, 1000)
cor_results <- data.frame(bin_size = numeric(), rho = numeric(), p_value = numeric())



for (i in bin_size){
  bin_size_new<- i +1
  bin_break <- seq(1, 44857, length.out= bin_size_new)

g4s_hist<- hist(entire_g4s_rdna[["actual_pG4CS_start"]], breaks = bin_break, plot = FALSE)
riz_hist <- hist(rlfs[["actual_RLFS_start"]], breaks = bin_break, plot = FALSE)

combined_data<- data.frame(
  bin_midpoints = (head(bin_break,-1) + tail(bin_break,-1))/2,
  pG4CS_counts = g4s_hist$counts,
  RIZ_counts = riz_hist$counts
)


fwrite(combined_data, paste("correlation_riz_g4_input_file_bin_size", i,".csv", sep=""), sep=",")


print(i)
print(cor.test(combined_data$pG4CS_counts, combined_data$RIZ_counts, method = "spearman"))

cor_result <- cor.test(combined_data$pG4CS_counts, combined_data$RIZ_counts, method = "spearman")

cor_results <- rbind(cor_results, data.frame(
  bin_size = i,
  rho = cor_result$estimate,
  p_value = cor_result$p.value
))


ggplot(cor_results, aes(x = bin_size, y = rho)) +
  geom_point() + geom_line() +
  labs(title = "Spearman's Correlation vs. Bin Size",
       x = "Bin Size (bp)", y = "Spearman's rho") +
  theme_minimal()


test<- ggplot(combined_data, aes(x = `RIZ_counts`, y = `pG4CS_counts`)) +
  geom_point(position = position_jitter(width = 2.0, height=0.2)) +
  geom_smooth(method = "lm", level= 0.95) +
  labs(x = "riz", y = "g4s", title = paste("Scatter plot of RIZ vs pG4CS")) + 
  theme_minimal()

# what i see is that rho value increase with decrease in bin size, for example its 0.7 when bins are only 30. 
# in either case this graph is not helping me to visualise anything or interpret anything. 
# maybe rdna regionwise will help instead of bins


ggsave(paste("correlation_riz_g4_input_file_bin_size", i, ".tiff", sep = ""), 
       plot = test, width = 18, height = 10, dpi = 150)
}

}




