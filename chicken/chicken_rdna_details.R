# ------------------------------------------------------------------------------
# This code is part of paper: In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus.
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj

# Purpose: This R code is designed to characterize nucleotide distribution of chicken rDNA sequence 
# Inputs:fasta file from https://www.ncbi.nlm.nih.gov/nuccore/KT445934 
# the code takes the input chicken rDNA fasta file, divide sequence into specific rDNA regions and count their nucleotide distribution
# Outputs: Nucleotide distribution csv file for chicken rDNA
# ------------------------------------------------------------------------------



#to gather information on chicken rDNA

library(Biostrings)
library(data.table)
library(tidyverse)

chicken_rdna<- readDNAStringSet("KT445934_chicken_rDNA_2017.fasta", format = "fasta", use.names = FALSE)
#default is fastq and use.names= FALSE will drop the header in the fasta file

chicken_rdna_seq<- as.character(chicken_rdna[[1]])
nchar(chicken_rdna_seq)
#11863

# from the NCBI entry KT445934.2

#Now i want to separate the rdna1 into subsections described by the authors
#1..1836- 5’ external transcribed spacer
#1837..3659 - 18S ribosomal RNA
#3660..6189 - internal transcribed spacer 1
#6190..6346 - 5.8S ribosomal RNA
#6347..7079 - internal transcribed spacer 2
#7080..11520- 28S ribosomal RNA
#11521..11863 - external transcribed spacer

#There is no IGS reported for chicken 

ets5<- subseq(chicken_rdna_seq, start= 1, end = 1836) 
s18<- subseq(chicken_rdna_seq, start= 1837, end = 3659) 
its1<- subseq(chicken_rdna_seq, start= 3660, end = 6189)
s5.8 <- subseq(chicken_rdna_seq, start= 6190, end = 6346)
its2<- subseq(chicken_rdna_seq, start= 6347, end = 7079)
s28<- subseq(chicken_rdna_seq, start= 7080, end = 11520)
ets3<- subseq(chicken_rdna_seq, start= 11521, end = 11863)

Name<- c("5'ETS_KT445934", "18S_KT445934", "ITS1_KT445934", "5.8S_KT445934", 
         "ITS2_KT445934", "28S_KT445934", "3'ETS_KT445934")


Sequences <- c(ets5, s18, its1, s5.8, its2, s28, ets3)

Details <- c("1_1836","1837_3659", "3660_6189", "6190_6346", 
             "6347_7079","7080_11520", "11521_11863")



rdna_chicken_dataset<- data.frame(Name, Sequences, Details)

attempt1<- data.frame((matrix(nrow = 0, ncol=5)))
for (j in 1: nrow(rdna_chicken_dataset)){
  a<- nchar (rdna_chicken_dataset[j,2])
  b<- str_count(rdna_chicken_dataset[j,2], "A")
  c<- str_count(rdna_chicken_dataset[j,2], "T")
  d<- str_count(rdna_chicken_dataset[j,2], "G")
  e<- str_count(rdna_chicken_dataset[j,2], "C")
  f <- c(a, b, c, d, e)
  attempt1[nrow(attempt1)+1,] <- f
}

colnames(attempt1) <- c("Total_nucleotides", "A", "T", "G", "C")
rdna_chicken_dataset_v2<- cbind(rdna_chicken_dataset, attempt1)


rdna_chicken_dataset_details_v1 <- rdna_chicken_dataset_v2 %>% 
  mutate(across(all_of(c("A","T","G","C")), function(x) x/Total_nucleotides*100, .names = "{col}%")) %>% 
  mutate(across(all_of(c(9,10,11,12)), function (x) round(x, 2))) %>% 
  select(Name,Sequences, Details, Total_nucleotides, "A","A%", "T","T%", "G", "G%","C","C%")



fwrite(rdna_chicken_dataset_details_v1, 
       file = "rdna_chicken_dataset_details_v1.csv")

#filter rows that contain a certain string

rdna_chicken_dataset_sequences<- rdna_chicken_dataset_details_v1 %>% select(1,2)


#create FASTA format
for (j in 1: nrow(rdna_chicken_dataset_sequences)){
  
  write(paste(">",rdna_chicken_dataset_sequences[j,1], sep=''),                                            
        file = "rDNA_KT445934_chicken_2017_sequence_fasta_format.txt",
        append = TRUE)
  
  write(rdna_chicken_dataset_sequences[j,2],                                            
        file = "rDNA_KT445934_chicken_2017_sequence_fasta_format.txt",
        append = TRUE)
  
  write("\n",                                            
        file = "rDNA_KT445934_chicken_2017_sequence_fasta_format.txt",
        append = TRUE)
  
}



rdna_chicken_dataset_details_v1<- fread("rdna_chicken_dataset_details_v1.csv", header = TRUE, sep = ",")
# select only nucleotide percent and identifier 
nucleotide<- rdna_chicken_dataset_details_v1 %>% select(Name,"A%","G%","C%", "T%")
nucleotide_new<- separate(nucleotide, Name, "Name", sep = "_")


#to convert columns to rows pivot_longer can work 
#pivot_longer() "lengthens" data, increasing the number of rows and decreasing the number of columns. The inverse transformation is pivot_wider()
#all columns except Name column will be reshaped

nucleotide_reshape <- pivot_longer(nucleotide_new, -Name, names_to = "Nucleotide", values_to = "Percent")
fwrite(nucleotide_reshape,  
       file= "nucleotide_rdna_chicken_dataset_sequences_graph_input1.csv")

#Please note below code do not work, because you can only drop one column 
#nar_reshape <- pivot_longer(nar, -Name, -Sequences,-Details, names_to = "Variable", values_to = "Value")

nucleotide_reshape$Name <- factor(nucleotide_reshape$Name, 
                                  levels = c("5'ETS", "18S", "ITS1", "5.8S", "ITS2","28S", "3'ETS"))




#The fill aesthetic is used to color the bars based on the the X variable.

rdna_chicken_nuceotide_distribution<- ggplot(nucleotide_reshape, aes(x = Name, y = Percent, fill = Nucleotide)) + 
  geom_bar(stat= 'identity', color = "black") +
  theme(axis.text.x = element_text(angle = 45, size = 20, hjust=1)) +
  labs(title= " Nucleotide distribution percent in chicken rDNA locus", 
       subtitle = "KT445934",
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



ggsave( "rdna_KT445934_chicken_nuceotide_distribution.tiff", 
        plot = rdna_chicken_nuceotide_distribution, width=15,height=10, dpi=150)

