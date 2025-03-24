##Aim: Find the total nucleotide distribution percent in KY962518_2018 NCBI entry that is for rDNA locus. 
## this calculation as been previously done in my earliar code. 
##however, after alot of back and forth we redefined the promoter and IGS boundary. 
## we simply took 2202 from IGS end side, because these were matching with UCSC hg38 nucleotides and added them to front of 5'ETS. We call them as promoter. 

library(Biostrings)
library(seqinr)
library(stringr)
library(data.table)
library(tidyverse)


##I wanted to add new promoter and redefined IGS to rDNA
hg38_rDNA<- read.fasta(file = "KY962518_added_3500nt_IGS_upstream_nontemplate.fasta", forceDNAtolower = FALSE, as.string = TRUE) #belongs to package seqinr

getwd#forceDNAtolower false is used to keep dna seq in uppercase, as. string to entire seq as one string
hg38_rDNA_seq<- hg38_rDNA[[1]]
nchar(hg38_rDNA_seq)
#48338


#Now i want to separate the rdna1 into subsections described by the authors
#1..3657- 5’ external transcribed spacer
#3658..5526 - 18S ribosomal RNA
#5527..6596 - internal transcribed spacer 1
#6597..6753 - 5.8S ribosomal RNA
#6754..7920 - internal transcribed spacer 2
#7921..12971- 28S ribosomal RNA
#12972..13332 - external transcribed spacer
#I assume 13333 ..44838 - intergenic spacer; IGS


promoter<- str_sub(hg38_rDNA_seq, start= 1299, end = 3500)
ets5<- str_sub(hg38_rDNA_seq, start = 3501, end = 7157)
s18<- str_sub(hg38_rDNA_seq, start= 7158, end = 9026)
its1<- str_sub(hg38_rDNA_seq, start= 9027, end = 10096)
s5.8 <- str_sub(hg38_rDNA_seq, start= 10097, end = 10253)
its2<- str_sub(hg38_rDNA_seq, start= 10254, end = 11420)
s28<- str_sub(hg38_rDNA_seq, start= 11421, end = 16471)
ets3<- str_sub(hg38_rDNA_seq, start= 16472, end = 16832)
igs<- str_sub(hg38_rDNA_seq, start= 16833, end = 48338)

Name<- c("Promoter_KY962518", "5'ETS_KY962518", "18S_KY962518", "ITS1_KY962518", "5.8S_KY962518", "ITS2_KY962518", 
         "28S_KY962518", "3'ETS_KY962518", "IGS_KY962518")


Sequences <- c(promoter,ets5, s18, its1, s5.8, its2, s28, ets3,igs)

Details <- c("1299_3500","3501_7157", "7158_9026", "9027_10096", "10097_10253",
             "10254_11420", "11421_16471", "16472_16832","16833_48338")



rdna_hg38_dataset<- data.frame(Name, Sequences, Details)

attempt1<- data.frame((matrix(nrow = 0, ncol=5)))
for (j in 1: nrow(rdna_hg38_dataset)){
  a<- nchar (rdna_hg38_dataset[j,2])
  b<- str_count(rdna_hg38_dataset[j,2], "A")
  c<- str_count(rdna_hg38_dataset[j,2], "T")
  d<- str_count(rdna_hg38_dataset[j,2], "G")
  e<- str_count(rdna_hg38_dataset[j,2], "C")
  f <- c(a, b, c, d, e)
  attempt1[nrow(attempt1)+1,] <- f
}

colnames(attempt1) <- c("Total_nucleotides", "A", "T", "G", "C")
rdna_hg38_dataset_v2<- cbind(rdna_hg38_dataset, attempt1)


rdna_hg38_dataset_details_v2 <- rdna_hg38_dataset_v2 %>% 
  mutate(across(all_of(c(5,6,7,8)), function(x) x/Total_nucleotides*100, .names = "{col}%")) %>% 
  mutate(across(all_of(c(9,10,11,12)), function (x) round(x, 2))) %>% 
  select(1,2,3,4,5,9,6,10,7,11,8,12)


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/output")

fwrite(rdna_hg38_dataset_details_v2, 
       file = "rdna_hg38_chr21_2018_dataset_details_v2.csv")

##filter rows that contain a certain string

rdna_hg38_dataset_sequences<- rdna_hg38_dataset_details_v2 %>% select(1,2)


##create FASTA format
for (j in 1: nrow(rdna_hg38_dataset_sequences)){
  
  write(paste(">",rdna_hg38_dataset_sequences[j,1], sep=''),                                            
        file = "rDNA_KY962518_2018_sequence_fasta_format.txt",
        append = TRUE)
  
  write(rdna_hg38_dataset_sequences[j,2],                                            
        file = "rDNA_KY962518_2018_sequence_fasta_format.txt",
        append = TRUE)
  
  write("\n",                                            
        file = "rDNA_KY962518_2018_sequence_fasta_format.txt",
        append = TRUE)
  
}


#
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rloop_and_rdna/human/one_rDNA_seq/output")
rdna_2018<- fread("rdna_hg38_chr21_2018_dataset_details_v2.csv", header = TRUE, sep = ",")
# select only nucleotide percent and identifier 
nucleotide<- rdna_2018 %>% select(Name,"A%","G%","C%", "T%")
nucleotide_new<- separate(nucleotide, Name, "Name", sep = "_")


#to convert columns to rows pivot_longer can work 
#pivot_longer() "lengthens" data, increasing the number of rows and decreasing the number of columns. The inverse transformation is pivot_wider()
##all columns except Name column will be reshaped

nucleotide_reshape <- pivot_longer(nucleotide_new, -Name, names_to = "Nucleotide", values_to = "Percent")
fwrite(nucleotide_reshape,  
       file= "nucleotide_rdna_hg38_chr21_dataset_sequences_graph_input2.csv")

#Please note below code do not work, because you can only drop one column 
#nar_reshape <- pivot_longer(nar, -Name, -Sequences,-Details, names_to = "Variable", values_to = "Value")

nucleotide_reshape$Name <- factor(nucleotide_reshape$Name, 
                                  levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2","28S", "3'ETS", "IGS" ))




#The fill aesthetic is used to color the bars based on the the X variable.

rdna_2018_sequences_nuceotide_distribution<- ggplot(nucleotide_reshape, aes(x = Name, y = Percent, fill = Nucleotide)) + 
  geom_bar(stat= 'identity', color = "black") +
  theme(axis.text.x = element_text(angle = 45, size = 20, hjust=1)) +
  labs(title= " Nucleotide distribution percent in human rDNA locus", 
       subtitle = "KY962518",
       x= "rDNA region", 
       y= "Nucleotide distribution percent")+
  geom_text(aes(label = Percent), position = position_stack(vjust = 0.5), size = 6)+ # Adjusted label positioning
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=20))



ggsave( "rdna_2018_sequences_nuceotide_distribution.tiff", 
        plot = rdna_2018_sequences_nuceotide_distribution, width=15,height=10, dpi=150)




