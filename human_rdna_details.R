##Aim: Find the total length of human rDNA repeat region and find the annotation of each of the rDNA element.
# go to link https://www.ncbi.nlm.nih.gov/nuccore/U13369.1?report=fasta
# copy the sequence in sublime text and remove all space 
# saved this sublime text as .fasta file
# here form which chr they have taken are not specified but when i blast few sequences it seems chr21.


library(Biostrings)
library(seqinr)
library(stringr)
library(data.table)
library(tidyverse)


test<- read.fasta(file = "test.fasta", forceDNAtolower = FALSE, as.string = TRUE) #belongs to package seqinr
test_seq<- test[[1]]
test_int<- str_remove_all(test_seq, "N")


rDNA<- read.fasta(file = "U13369_1_humanrDNA_2010.fasta", forceDNAtolower = FALSE, as.string = TRUE) #belongs to package seqinr

#forceDNAtolower false is used to keep dna seq in uppercase, as. string to entire seq as one string
rDNA_seq<- rDNA[[1]]
nchar(rDNA_seq)
#42999

#Now i want to separate the rdna1 into subsections described by the authors
#1 to 3656- 5’ ETS (5' external transcribed spacer)
#1 to 694bp - may be promoter
#if we think promoter is 1000, then imagine if we add 306 upstream to 1 to 694. this 306 i can remove from intergenic spacer element
#695 to 3656 – Site of Transcription
#3657 to 5527 - 18S ribosomal RNA
#5528 to 6622 – Internal Transcriber spacer 1 (ITS1)
#6623 to 6779 – 5.8S ribosomal RNA
#6780 to 7934 - Internal transcribed spacer 2 (ITS2)
#7935 to 12969- 28S ribosomal RNA
#12970 to 13314- 3' external transcribed spacer (3’ETS)
#13315 to 42999- intergenic spacer; IGS
#13334 to 13350- terminator
#42694 to 42999- 306bp which is added to upsteam of promoter

ets5<- str_sub(rDNA_seq, start = 1, end = 3656)
tss<- str_sub(rDNA_seq, start= 695, end = 3656)
s18<- str_sub(rDNA_seq, start= 3657, end = 5527)
its1<- str_sub(rDNA_seq, start= 5528, end = 6622)
s5.8 <- str_sub(rDNA_seq, start= 6623, end = 6779)
its2<- str_sub(rDNA_seq, start= 6780, end = 7934)
s28<- str_sub(rDNA_seq, start= 7935, end = 12969)
ets3<- str_sub(rDNA_seq, start= 12970, end = 13314)
igs<- str_sub(rDNA_seq, start= 13315, end = 42999) # contains NNN, Y, R like characters
p1<- str_sub(rDNA_seq, start= 42694, end = 42999) #total 306 bp
p2<- str_sub(rDNA_seq, start= 1, end= 694) # total 694 bp
promoter<- paste(p1, p2, sep= "") 

Name<- c("5'ETS", "TSS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS", "IGS", "P1", "P2", "Promoter")


Sequences <- c(ets5, tss, s18, its1, s5.8, its2, s28, ets3, igs, p1, p2, promoter)

Details <- c("1_3656", "695_3656", "3657_5527", "5528_6622", "6623_6779",
             "6980_7934", "7935_12969", "12970_13314","13315_42999",
             "42694_42999", "1_694", "p1_to_p2")



rdna_dataset<- data.frame(Name, Sequences, Details)

attempt1<- data.frame((matrix(nrow = 0, ncol=5)))
for (j in 1: nrow(rdna_dataset)){
  a<- nchar (rdna_dataset[j,2])
  b<- str_count(rdna_dataset[j,2], "A")
  c<- str_count(rdna_dataset[j,2], "T")
  d<- str_count(rdna_dataset[j,2], "G")
  e<- str_count(rdna_dataset[j,2], "C")
  f <- c(a, b, c, d, e)
  attempt1[nrow(attempt1)+1,] <- f
}

colnames(attempt1) <- c("Total_nucleotides", "A", "T", "G", "C")
rdna_dataset_v1<- cbind(rdna_dataset, attempt1)


rdna_dataset_details_v1 <- rdna_dataset_v1 %>% 
  mutate(across(all_of(c(5,6,7,8)), function(x) x/Total_nucleotides*100, .names = "{col}%")) %>% 
  mutate(across(all_of(c(9,10,11,12)), function (x) round(x, 2))) %>% 
  select(1,2,3,4,5,9,6,10,7,11,8,12)


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rloop_and_rdna/one_rDNA_seq/output")

fwrite(rdna_dataset_details_v1, 
       file = "rdna_dataset_details_2010_v1.csv")

##filter rows that contain a certain string

rdna_dataset_sequences<- rdna_dataset_details_v1 %>% select(1,2)


##create FASTA format
for (j in 1: nrow(rdna_dataset_sequences)){
  
  write(paste(">",rdna_dataset_sequences[j,1], sep=''),                                            
        file = "rDNA_sequence_2010_fasta_format.txt",
        append = TRUE)
  
  write(rdna_dataset_sequences[j,2],                                            
        file = "rDNA_sequence_2010_fasta_format.txt",
        append = TRUE)
  
  write("\n",                                            
        file = "rDNA_sequence_2010_fasta_format.txt",
        append = TRUE)
  
}


#I found another link https://www.ncbi.nlm.nih.gov/nuccore/KY962518.1?report=genbank
#where the defined rDNA in chr21 in human using reference hg38

#I copy pasted sequence in sublime text, removed spaces and saved as fasta format KY962518_1_humanrDNA_2018.fasta.
#KY962518_1_humanrDNA_2018.fasta contains one rDNA repeat from chr21.

hg38_rDNA<- read.fasta(file = "KY962518_1_humanrDNA_2018.fasta", forceDNAtolower = FALSE, as.string = TRUE) #belongs to package seqinr

#forceDNAtolower false is used to keep dna seq in uppercase, as. string to entire seq as one string
hg38_rDNA_seq<- hg38_rDNA[[1]]
nchar(hg38_rDNA_seq)
#44838


#Now i want to separate the rdna1 into subsections described by the authors
#1..3657- 5’ external transcribed spacer
#3658..5526 - 18S ribosomal RNA
#5527..6596 - internal transcribed spacer 1
#6597..6753 - 5.8S ribosomal RNA
#6754..7920 - internal transcribed spacer 2
#7921..12971- 28S ribosomal RNA
#12972..13332 - external transcribed spacer
#I assume 13333 ..44838 - intergenic spacer; IGS


ets5_hg38<- str_sub(hg38_rDNA_seq, start = 1, end = 3657)
s18_hg38<- str_sub(hg38_rDNA_seq, start= 3658, end = 5526)
its1_hg38<- str_sub(hg38_rDNA_seq, start= 5527, end = 6596)
s5.8_hg38 <- str_sub(hg38_rDNA_seq, start= 6597, end = 6753)
its2_hg38<- str_sub(hg38_rDNA_seq, start= 6754, end = 7920)
s28_hg38<- str_sub(hg38_rDNA_seq, start= 7921, end = 12971)
ets3_hg38<- str_sub(hg38_rDNA_seq, start= 12972, end = 13332)
igs_hg38<- str_sub(hg38_rDNA_seq, start= 13333, end = 44838)
p1_hg38<- str_sub(hg38_rDNA_seq, start= 44339, end = 44838)
p2_hg38<- str_sub(hg38_rDNA_seq, start= 1, end= 500)
promoter_hg38<- paste(p1_hg38, p2_hg38, sep= "")

Name<- c("5'ETS_hg38", "18S_hg38", "ITS1_hg38", "5.8S_hg38", "ITS2_hg38", 
         "28S_hg38", "3'ETS_hg38", "IGS_hg38", "P1_hg38", "P2_hg38", "Promoter_hg38")


Sequences <- c(ets5_hg38, s18_hg38, its1_hg38, s5.8_hg38, its2_hg38, s28_hg38, ets3_hg38, 
               igs_hg38, p1_hg38, p2_hg38, promoter_hg38)

Details <- c("1_3657", "3658_5526", "5527_6596", "6597_6753",
             "6754_7920", "7921_12971", "12972_13332","13333_44838",
             "44339_44838", "1_500", "add p1_hg38 to p2_hg38")



rdna_hg38_dataset<- data.frame(Name, Sequences, Details)

attempt1<- data.frame((matrix(nrow = 0, ncol=5)) )
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
rdna_hg38_dataset_v1<- cbind(rdna_hg38_dataset, attempt1)


rdna_hg38_dataset_details_v1 <- rdna_hg38_dataset_v1 %>% 
  mutate(across(all_of(c(5,6,7,8)), function(x) x/Total_nucleotides*100, .names = "{col}%")) %>% 
  mutate(across(all_of(c(9,10,11,12)), function (x) round(x, 2))) %>% 
  select(1,2,3,4,5,9,6,10,7,11,8,12)


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rloop_and_rdna/output")

fwrite(rdna_hg38_dataset_details_v1, 
       file = "rdna_hg38_chr21_2018_dataset_details_v1.csv")

##filter rows that contain a certain string

rdna_hg38_dataset_sequences<- rdna_hg38_dataset_details_v1 %>% select(1,2)


##create FASTA format
for (j in 1: nrow(rdna_hg38_dataset_sequences)){
  
  write(paste(">",rdna_hg38_dataset_sequences[j,1], sep=''),                                            
        file = "rDNA_hg38_chr21_2018_sequence_fasta_format.txt",
        append = TRUE)
  
  write(rdna_hg38_dataset_sequences[j,2],                                            
        file = "rDNA_hg38_chr21_2018_sequence_fasta_format.txt",
        append = TRUE)
  
  write("\n",                                            
        file = "rDNA_hg38_chr21_2018_sequence_fasta_format.txt",
        append = TRUE)
  
}




##put these sequences in qmrlfs_finder and save results in QmRLFS results folder 2018
#Open terminal, go to downloads where qmrlfs.py file exits, copy paste rDNA_hg38_chr21_2018_sequence_fasta_format.txt
#in downloads and run the code 
#default output is always txt file which is "\t" separated file 
#python QmRLFS-finder.py -bed -i rDNA_hg38_chr21_2018_sequence_fasta_format.txt -o rDNA_hg38_chr21_2018_qmrlfs
#out will be saved as rDNA_hg38_chr21_2018_qmrlfs.out.bed

#similarly run this code to download table format
#python QmRLFS-finder.py -i rDNA_hg38_chr21_2018_sequence_fasta_format.txt -o rDNA_hg38_chr21_2018_qmrlfs

#i also saved the rloop_db pdf file in qmrlfs results 2018 folder


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rloop_and_rdna/human/one_rDNA_seq/output")
rdna_2018<- fread("rdna_hg38_chr21_2018_dataset_details_v1.csv", header = TRUE, sep = ",")
# select only nucleotide percent and identifier 
nucleotide<- rdna_2018 %>% select(Name,"A%","T%","G%","C%")
nucleotide<- nucleotide[1:8,]
nucleotide_new<- separate(nucleotide, Name, "Name", sep = "_")


#to convert columns to rows pivot_longer can work 
#pivot_longer() "lengthens" data, increasing the number of rows and decreasing the number of columns. The inverse transformation is pivot_wider()
##all columns except Name column will be reshaped

nucleotide_reshape <- pivot_longer(nucleotide_new, -Name, names_to = "Nucleotide", values_to = "Percent")
fwrite(nucleotide_reshape,  
       file= "nucleotide_rdna_hg38_chr21_dataset_sequences_graph_input.csv")

#Please note below code do not work, because you can only drop one column 
#nar_reshape <- pivot_longer(nar, -Name, -Sequences,-Details, names_to = "Variable", values_to = "Value")

nucleotide_reshape$Name <- factor(nucleotide_reshape$Name, 
                                            levels = c("5'ETS", "18S", "ITS1", "5.8S", "ITS2","28S", "3'ETS", "IGS" ))




#The fill aesthetic is used to color the bars based on the the X variable.

rdna_2018_sequences_nuceotide_distribution<- ggplot(nucleotide_reshape, aes(x = Name, y = Percent, fill = Nucleotide)) + 
  geom_bar(stat= 'identity', color = "black") +
  theme(axis.text.x = element_text(angle = 45, size = 10)) +
  labs(title= " Nucleotide distribution percent in human rDNA locus", x= "rDNA region", y= "Nucleotide distribution percent")+
   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      text = element_text(size = 30),
      axis.line = element_line(color = "black"),
      theme(panel.grid = element_blank()),
      axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25)) +  # Center Y-axis title
  theme(axis.ticks.y = element_line(color = "black"))

ggsave( "rdna_2018_sequences_nuceotide_distribution.tiff", 
        plot = rdna_2018_sequences_nuceotide_distribution, width=15,height=10, dpi=150)



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

attempt1<- data.frame((matrix(nrow = 0, ncol=5)) )
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




##miscellaneous activity 
##i just want to see if rDNA exits in ensemble genelist or not

{
gene<-fread("Ensemble_gene.txt", sep="\t", header = TRUE) # taken from ojashpi_rnaseq_folder, copy pasted and then delete from this folder
#280306
gene1<- gene %>% group_by(`Gene type`) %>% count()
fwrite(gene1, "ensemble_groupby_table.csv", sep = ",")

gene_rRNA<- gene %>% filter(`Gene type`== "rRNA")
colnames(gene_rRNA)[1] <- "chr"

fwrite(gene1, "ensemble_groupby_table.csv", sep = ",")

table(gene_rRNA$chr)

#1      11   14   17    2   20   21   22    3    4    5    6    8    9    X 
#1734    2    1    3    1    5  843    1    1    2    1    1    1    1    2 

filter_values<- c("13", "14", "15","21", "22")
filtered_rRNA <- gene_rRNA  %>% filter(chr %in% filter_values)
filtered_rRNA_coi<- filtered_rRNA %>% select(-c(7))
filtered_rRNA_coi2<- unique(filtered_rRNA_coi)

fwrite(gene_rRNA, "ensemble_gene_rDNA.csv", sep = ",")
fwrite(filtered_rRNA_coi2, "ensemble_filtered_rDNA.csv", sep = ",")
}
