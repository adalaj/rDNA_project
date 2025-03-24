##Find rloop in mouse rdna region
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rloop_and_rdna/mouse")
library(Biostrings)
library(stringr)
library(data.table)
library(tidyverse)
library(seqinr)


#ref Sayers E.W., Bolton E.E., Brister J.R., Canese K., Chan J., Comeau D.C.et al.
#Database resources of the national center for biotechnology information.
#Nucleic Acids Res. 2022; 50: D20-D26

#Human KY962518.1 (44,838 nt) and mouse Genbank BK000964.3 (45,306 nt) reference rDNA sequences are available from NCBI Genbank (13) and are used in the literature as consensus reference rDNA sequences for their respective species (8, 14, 15, 16)

#I found mouse link https://www.ncbi.nlm.nih.gov/nuccore/BK000964?report=genbank
#rdna region is in 12, 15, 17,18 and 19 in mouse. I guess the sequency is in chr 17 in mouse ncbi genome assembly.

#I clicked on fasta, then select send to file and hit enter. I renamed from sequence.fasta to BK000964_mouse_rDNA_2013.fasta.
#BK000964_mouse_rDNA_2013.fastacontains one mouse rDNA repeat from chr17 i believe.

mouse_rDNA<- read.fasta(file = "BK000964_mouse_rDNA_2013.fasta", forceDNAtolower = FALSE, as.string = TRUE) #belongs to package seqinr

#forceDNAtolower false is used to keep dna seq in uppercase, as. string to entire seq as one string
mouse_rDNA_seq<- mouse_rDNA[[1]]
nchar(mouse_rDNA_seq)
#45306


#Now i want to separate the rdna1 into subsections described by the authors
#1..4007- 5’ external transcribed spacer
#4008..5877 - 18S ribosomal RNA
#5878..6877 - internal transcribed spacer 1
#6878..7034 - 5.8S ribosomal RNA
#7035..8122 - internal transcribed spacer 2
#8123..12852- 28S ribosomal RNA
#12853..13403 - 3' external transcribed spacer
#13404 ..45306 - intergenic spacer; IGS


ets5_mouse<- str_sub(mouse_rDNA_seq, start = 1, end = 4007)
s18_mouse<- str_sub(mouse_rDNA_seq, start= 4008, end = 5877)
its1_mouse<- str_sub(mouse_rDNA_seq, start= 5878, end = 6877)
s5.8_mouse <- str_sub(mouse_rDNA_seq, start= 6878, end = 7034)
its2_mouse<- str_sub(mouse_rDNA_seq, start= 7035, end = 8122)
s28_mouse<- str_sub(mouse_rDNA_seq, start= 8123, end = 12852)
ets3_mouse<- str_sub(mouse_rDNA_seq, start= 12853, end = 13403)
igs_mouse<- str_sub(mouse_rDNA_seq, start= 13404, end = 45306)


Name<- c("5'ETS_mouse", "18S_mouse", "ITS1_mouse", "5.8S_mouse", 
         "ITS2_mouse", "28S_mouse", "3'ETS_mouse", "IGS_mouse")


Sequences <- c(ets5_mouse, s18_mouse, its1_mouse, s5.8_mouse, 
               its2_mouse, s28_mouse, ets3_mouse, igs_mouse)

Details <- c("1_4007", "4008_5877", "5878_6877", "6878_7034",
             "7035_8122", "8123_12852", "12853_13403","13404_45306")



rdna_mouse_dataset<- data.frame(Name, Sequences, Details)

attempt1<- data.frame((matrix(nrow = 0, ncol=5)) )
for (j in 1: nrow(rdna_mouse_dataset)){
  a<- nchar (rdna_mouse_dataset[j,2])
  b<- str_count(rdna_mouse_dataset[j,2], "A")
  c<- str_count(rdna_mouse_dataset[j,2], "T")
  d<- str_count(rdna_mouse_dataset[j,2], "G")
  e<- str_count(rdna_mouse_dataset[j,2], "C")
  f <- c(a, b, c, d, e)
  attempt1[nrow(attempt1)+1,] <- f
}

colnames(attempt1) <- c("Total_nucleotides", "A", "T", "G", "C")
rdna_mouse_dataset_v1<- cbind(rdna_mouse_dataset, attempt1)


rdna_mouse_dataset_details_v1 <- rdna_mouse_dataset_v1 %>% 
  mutate(across(all_of(c(5,6,7,8)), function(x) x/Total_nucleotides*100, .names = "{col}%")) %>% 
  mutate(across(all_of(c(9,10,11,12)), function (x) round(x, 2))) %>% 
  select(1,2,3,4,5,9,6,10,7,11,8,12)


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rloop_and_rdna/mouse/output")

fwrite(rdna_mouse_dataset_details_v1, 
       file = "rdna_mouse_2013_dataset_details_v1.csv")

##filter rows that contain a certain string

rdna_mouse_dataset_sequences<- rdna_mouse_dataset_details_v1 %>% select(1,2)


##create FASTA format
for (j in 1: nrow(rdna_mouse_dataset_sequences)){
  
  write(paste(">",rdna_mouse_dataset_sequences[j,1], sep=''),                                            
        file = "rDNA_mouse_2013_sequence_fasta_format.txt",
        append = TRUE)
  
  write(rdna_mouse_dataset_sequences[j,2],                                            
        file = "rDNA_mouse_2013_sequence_fasta_format.txt",
        append = TRUE)
  
  write("\n",                                            
        file = "rDNA_mouse_2013_sequence_fasta_format.txt",
        append = TRUE)
  
}

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rloop_and_rdna/mouse/output")
rdna_2013<- fread("rdna_mouse_2013_dataset_details_v1.csv", header = TRUE, sep = ",")
# select only nucleotide percent and identifier 
nucleotide<- rdna_2013 %>% select(Name,"A%","T%","G%","C%")
nucleotide_new<- separate(nucleotide, Name, "Name", sep = "_")


nucleotide_reshape <- pivot_longer(nucleotide_new, -Name, names_to = "Nucleotide", values_to = "Percent")
fwrite(nucleotide_reshape,  
       file= "nucleotide_rdna_mouse_2013_dataset_sequences_graph_input.csv")

#Please note below code do not work, because you can only drop one column 
#nar_reshape <- pivot_longer(nar, -Name, -Sequences,-Details, names_to = "Variable", values_to = "Value")

nucleotide_reshape$Name <- factor(nucleotide_reshape$Name, 
                                  levels = c("5'ETS", "18S", "ITS1", "5.8S", "ITS2","28S", "3'ETS", "IGS" ))




#The fill aesthetic is used to color the bars based on the the X variable.

rdna_2013_sequences_nucleotide_distribution<- ggplot(nucleotide_reshape, aes(x = Name, y = Percent, fill = Nucleotide)) + 
  geom_bar(stat= 'identity', color = "black") +
  theme(axis.text.x = element_text(angle = 45, size = 10)) +
  labs(title= " Nucleotide distribution percent in mouse rDNA locus", x= "rDNA region", y= "Nucleotide distribution percent")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        theme(panel.grid = element_blank()),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25)) +  # Center Y-axis title
  theme(axis.ticks.y = element_line(color = "black"))

ggsave( "rdna_mouse_2013_sequences_nucleotide_distribution.tiff", 
        plot = rdna_2013_sequences_nucleotide_distribution, width=15,height=10, dpi=150)



##put these sequences in qmrlfs_finder and save results in QmRLFS results folder 2018
#Open terminal, go to downloads where qmrlfs.py file exits, copy paste rDNA_mouse_2013_sequence_fasta_format.TXT
#in downloads and run the code 
#default output is always txt file which is "\t" separated file 
#python QmRLFS-finder.py -bed -i rDNA_mouse_2013_sequence_fasta_format.txt -o rDNA_mouse_2013_qmrlfs
#out will be saved as rDNA_mouse_2013_qmrlfs.out.bed

#similarly run this code to download table format
#python QmRLFS-finder.py -i rDNA_mouse_2013_sequence_fasta_format.txt -o rDNA_mouse_2013_qmrlfs

#i also saved the rloop_db pdf file 



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






##plotting with IGS on start and promoter at end.
###Plot4 showing entire rDNA and begining of next



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



