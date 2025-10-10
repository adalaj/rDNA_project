
#to gather information on monkey rDNA
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/monkey")


library(Biostrings)
library(data.table)
library(tidyverse)

monkey_rdna<- readDNAStringSet("KX061890_monkey_entire_rdna_2018_nontemplate.fasta", format = "fasta", use.names = FALSE)
#https://www.ncbi.nlm.nih.gov/nuccore/KX061890
#the header i removed KX061890.1 Macaca mulatta precursor 48S ribosomal RNA gene, external transcribed spacer, 18S ribosomal RNA gene, internal transcribed spacer 1, 5.8S ribosomal RNA gene, internal transcribed spacer 2, 28S ribosomal RNA gene, and external transcribed spacer, complete sequence
#default is fastq and use.names= FALSE will drop the header in the fasta file

monkey_rdna_seq<- as.character(monkey_rdna[[1]])
nchar(monkey_rdna_seq)
#41735

# from the NCBI entry KX061890.1

#Now i want to separate the rdna1 into subsections described by the authors
#1..3640- 5’ external transcribed spacer
#3641..5508 - 18S ribosomal RNA
#5509..6535 - internal transcribed spacer 1
#6536..6692 - 5.8S ribosomal RNA
#6693..7863 - internal transcribed spacer 2
#7864..12648- 28S ribosomal RNA
#12649..12979 - 3' external transcribed spacer
#12980..41735- "intergenic spacer; IGS"

#For the moment i dont want to include IGS for monkey 

ets5<- subseq(monkey_rdna_seq, start= 1, end = 3640) 
s18<- subseq(monkey_rdna_seq, start= 3641, end = 5508) 
its1<- subseq(monkey_rdna_seq, start= 5509, end = 6535)
s5.8 <- subseq(monkey_rdna_seq, start= 6536, end = 6692)
its2<- subseq(monkey_rdna_seq, start= 6693, end = 7863)
s28<- subseq(monkey_rdna_seq, start= 7864, end = 12648)
ets3<- subseq(monkey_rdna_seq, start= 12649, end = 12979)

# I found https://www.ncbi.nlm.nih.gov/nuccore/NR_146166.1?report=genbank
# it has 18s sequence of monkey that was added to NCBI on 2020. 

test<- readDNAStringSet("NR_146166_1_moneky_18S_2020.fasta", format = "fasta", use.names = FALSE)
test_seq<- as.character(test[[1]])
nchar(test_seq)
#1869 bp same as s18

subject<- DNAString(s18)
pattern<- DNAString(test)
alignment <- pairwiseAlignment(subject, pattern, substitutionMatrix = NULL, gapOpening = -5, gapExtension = -2)
pid <- pid(alignment)  # returns percentage identity
print(pid)
#99.9465

s18_NR_146166.1<- test





Name<- c("5'ETS_KX061890", "18S_KX061890", "ITS1_KX061890", "5.8S_KX061890", 
         "ITS2_KX061890", "28S_KX061890", "3'ETS_KX061890", "s18_NR_146166.1")


Sequences <- c(ets5, s18, its1, s5.8, its2, s28, ets3, s18_NR_146166.1)

Details <- c("1_3640","3641_5508", "5509_6535", "6536_6692", 
             "6693_7863","7864_12648", "12649_12979", "s18_NR_146166.1")



rdna_monkey_dataset<- data.frame(Name, Sequences, Details)

attempt1<- data.frame((matrix(nrow = 0, ncol=5)))
for (j in 1: nrow(rdna_monkey_dataset)){
  a<- nchar (rdna_monkey_dataset[j,2])
  b<- str_count(rdna_monkey_dataset[j,2], "A")
  c<- str_count(rdna_monkey_dataset[j,2], "T")
  d<- str_count(rdna_monkey_dataset[j,2], "G")
  e<- str_count(rdna_monkey_dataset[j,2], "C")
  f <- c(a, b, c, d, e)
  attempt1[nrow(attempt1)+1,] <- f
}

colnames(attempt1) <- c("Total_nucleotides", "A", "T", "G", "C")
rdna_monkey_dataset_v2<- cbind(rdna_monkey_dataset, attempt1)


rdna_monkey_dataset_details_v1 <- rdna_monkey_dataset_v2 %>% 
  mutate(across(all_of(c("A","T","G","C")), function(x) x/Total_nucleotides*100, .names = "{col}%")) %>% 
  mutate(across(all_of(c(9,10,11,12)), function (x) round(x, 2))) %>% 
  select(Name,Sequences, Details, Total_nucleotides, "A","A%", "T","T%", "G", "G%","C","C%")



fwrite(rdna_monkey_dataset_details_v1, 
       file = "rdna_monkey_dataset_details_v1.csv")

##filter rows that contain a certain string

rdna_monkey_dataset_sequences<- rdna_monkey_dataset_details_v1 %>% select(1,2)


##create FASTA format
for (j in 1: nrow(rdna_monkey_dataset_sequences)){
  
  write(paste(">",rdna_monkey_dataset_sequences[j,1], sep=''),                                            
        file = "rDNA_KX061890_monkey_2018_sequence_fasta_format.txt",
        append = TRUE)
  
  write(rdna_monkey_dataset_sequences[j,2],                                            
        file = "rDNA_KX061890_monkey_2018_sequence_fasta_format.txt",
        append = TRUE)
  
  write("\n",                                            
        file = "rDNA_KX061890_monkey_2018_sequence_fasta_format.txt",
        append = TRUE)
  
}


#created two fasta format from 5ets to 3'ets one from KX061890 
ets5_to_ets3<- subseq(monkey_rdna_seq, start = 1, end=12979)
nchar(ets5_to_ets3)
#12979


monkey_nontemplate <- DNAStringSet(ets5_to_ets3)

names(monkey_nontemplate) <- "nontemplate_KX061890"
writeXStringSet(monkey_nontemplate, "nontemplate_monkey_5ets_KX061890_3ets.fasta")


#another form replacing the existing 18s with new 18s from NR_146166.1
monkey_modified<- paste0(ets5, s18_NR_146166.1, its1, s5.8, its2, s28, ets3)
nchar(monkey_modified)
#[1] 12980

monkey_modified_nontemplate <- DNAStringSet(monkey_modified)

names(monkey_modified_nontemplate) <- "nontemplate_KX061890_and_NR_146166"
writeXStringSet(monkey_modified_nontemplate, "nontemplate_monkey_5ets_KX061890_and_NR_146166_3ets.fasta")


#going forward i will use both 18s to see if a see a presence of RLFS, G4 or imotif




###############################pending######################
#I have not run the below lines 

rdna_monkey_dataset_details_v1<- fread("rdna_monkey_dataset_details_v1.csv", header = TRUE, sep = ",")
# select only nucleotide percent and identifier 
nucleotide<- rdna_monkey_dataset_details_v1 %>% select(Name,"A%","G%","C%", "T%")
nucleotide_new<- separate(nucleotide, Name, "Name", sep = "_")


#to convert columns to rows pivot_longer can work 
#pivot_longer() "lengthens" data, increasing the number of rows and decreasing the number of columns. The inverse transformation is pivot_wider()
##all columns except Name column will be reshaped

nucleotide_reshape <- pivot_longer(nucleotide_new, -Name, names_to = "Nucleotide", values_to = "Percent")
fwrite(nucleotide_reshape,  
       file= "nucleotide_rdna_monkey_dataset_sequences_graph_input1.csv")

#Please note below code do not work, because you can only drop one column 
#nar_reshape <- pivot_longer(nar, -Name, -Sequences,-Details, names_to = "Variable", values_to = "Value")

nucleotide_reshape$Name <- factor(nucleotide_reshape$Name, 
                                  levels = c("5'ETS", "18S", "ITS1", "5.8S", "ITS2","28S", "3'ETS"))




#The fill aesthetic is used to color the bars based on the the X variable.

rdna_monkey_nuceotide_distribution<- ggplot(nucleotide_reshape, aes(x = Name, y = Percent, fill = Nucleotide)) + 
  geom_bar(stat= 'identity', color = "black") +
  theme(axis.text.x = element_text(angle = 45, size = 20, hjust=1)) +
  labs(title= " Nucleotide distribution percent in monkey rDNA locus", 
       subtitle = "KX061890",
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



ggsave( "rdna_KX061890_monkey_nuceotide_distribution.tiff", 
        plot = rdna_monkey_nuceotide_distribution, width=15,height=10, dpi=150)

