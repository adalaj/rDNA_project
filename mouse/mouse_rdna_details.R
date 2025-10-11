# ------------------------------------------------------------------------------
# This code accompanies the paper:
# "In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus"
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
# This R script processes and visualizes R-loop forming sequence (RLFS) predictions 
# across the mouse rDNA locus obtained from the QmRLFS-finder algorithm.
# It assigns each predicted RLFS to defined rDNA subregions and quantifies their 
# strand-specific and normalized distributions for both analytical and visual representation.
#
# Specifically, it:
#   1. Runs QmRLFS-finder (v1.5) on the mouse rDNA FASTA sequence 
#      (GenBank: BK000964) to identify RLFSs in BED format.
#   2. Parses and processes the RLFS output to calculate start–end coordinates, 
#      strand specificity, and RLFS lengths.
#   3. Assigns each RLFS to rDNA subcomponents (5′ETS, 18S, ITS1, 5.8S, ITS2, 28S, 3′ETS) 
#      using a rule-based system that counts RLFSs in the region where they initiate, 
#      ensuring junction-spanning RLFSs are not lost to boundary truncation.
#   4. Computes raw and normalized RLFS counts per region, including junctional regions, 
#      and generates tabular outputs suitable for statistical or graphical analyses.
#   5. Visualizes RLFS distributions using bar plots (overall, junction-excluded, 
#      and strandwise) and custom ideograms with karyoploteR for template and 
#      non-template strands.
#
# Inputs:
#   - mouse rDNA FASTA sequence (GenBank: BK000964)
#   - BED output file from QmRLFS-finder (v1.5)
#
# Outputs:
#   - Annotated CSVs summarizing RLFS counts per rDNA region (with and without junctions)
#   - Strandwise normalized RLFS count tables
#   - Publication-quality TIFF plots showing RLFS distributions
#   - KaryoploteR ideograms illustrating RLFS positions on both strands
#
# ------------------------------------------------------------------------------


library(Biostrings)
library(seqinr)
library(stringr)
library(data.table)
library(tidyverse)


#ref Sayers E.W., Bolton E.E., Brister J.R., Canese K., Chan J., Comeau D.C.et al.
#Database resources of the national center for biotechnology information.
#Nucleic Acids Res. 2022; 50: D20-D26


#I found mouse link https://www.ncbi.nlm.nih.gov/nuccore/BK000964?report=genbank
#rdna region is in 12, 15, 17,18 and 19 in mouse. I guess the sequency is in chr 17 in mouse ncbi genome assembly.

#I clicked on fasta, then select send to file and hit enter. I renamed from sequence.fasta to BK000964_mouse_rDNA_2013.fasta.
#BK000964_mouse_rDNA_2013.fastacontains one mouse rDNA repeat from chr17 i believe.




mouse_rDNA<- read.fasta(file = "BK000964_mouse_rDNA_2013.fasta", forceDNAtolower = FALSE, as.string = TRUE) #belongs to package seqinr
# previously I manually added 3500 bp from IGS to 5'ETS. out of which 2202 bp is dedicated to promoter. 
# however I can see RLFS in IGS boreder next to promoter in visualisation using karyoplote
# so decided to add 5000nucleotides from IGS to front, which will smooth my calculation and out of which I will assign 3000 to promoter this time.


#forceDNAtolower false is used to keep dna seq in uppercase, as. string to entire seq as one string
mouse_rDNA_seq<- mouse_rDNA[[1]]
nchar(mouse_rDNA_seq)
#45306

#i performed similar calculation of mouse 45306bp rDNA which you will find in git commit history for mouse_rdna_rloop.r

igs_and_promoter<- str_sub(mouse_rDNA_seq,40307,45306)
nchar(igs_and_promoter)
#[1] 5000

new_mouse_rDNA_seq <- str_c(igs_and_promoter, mouse_rDNA_seq)
nchar(new_mouse_rDNA_seq)

write.fasta(
  sequences = new_mouse_rDNA_seq,
  names = "BK000964",
  file.out = "BK000964_added_5000nt_IGS_upstream_nontemplate.fasta")





#Now i want to separate the rdna1 into subsections described by the authors
#1..4007- 5’ external transcribed spacer
#4008..5877 - 18S ribosomal RNA
#5878..6877 - internal transcribed spacer 1
#6878..7034 - 5.8S ribosomal RNA
#7035..8122 - internal transcribed spacer 2
#8123..12852- 28S ribosomal RNA
#12853..13403 - 3' external transcribed spacer
#13404 ..45306 - intergenic spacer; IGS




promoter_mouse<- str_sub(new_mouse_rDNA_seq, start= 2001, end = 5000) #3000
ets5_mouse<- str_sub(new_mouse_rDNA_seq, start = 5001, end = 9007) #4007
s18_mouse<- str_sub(new_mouse_rDNA_seq, start= 9008, end = 10877)
its1_mouse<- str_sub(new_mouse_rDNA_seq, start= 10878, end = 11877)
s5.8_mouse <- str_sub(new_mouse_rDNA_seq, start= 11878, end = 12034)
its2_mouse<- str_sub(new_mouse_rDNA_seq, start= 12035, end = 13122)
s28_mouse<- str_sub(new_mouse_rDNA_seq, start= 13123, end = 17852)
ets3_mouse<- str_sub(new_mouse_rDNA_seq, start= 17853, end = 18403)
igs_mouse<- str_sub(new_mouse_rDNA_seq, start= 18404, end = 50306)


Name<- c("Promoter_BK000964","5'ETS_BK000964", "18S_BK000964", "ITS1_BK000964", "5.8S_BK000964", 
         "ITS2_BK000964", "28S_BK000964", "3'ETS_BK000964", "IGS_BK000964")


Sequences <- c(promoter_mouse,ets5_mouse, s18_mouse, its1_mouse, s5.8_mouse, 
               its2_mouse, s28_mouse, ets3_mouse, igs_mouse)

Details <- c("2001_5000", "5001_9007", "9008_10877", "10878_11877", "11878_12034",
             "12035_13122", "13123_17852", "17853_18403","18404_50306")



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


rdna_mouse_dataset_details_v2 <- rdna_mouse_dataset_v1 %>% 
  mutate(across(all_of(c(5,6,7,8)), function(x) x/Total_nucleotides*100, .names = "{col}%")) %>% 
  mutate(across(all_of(c(9,10,11,12)), function (x) round(x, 2))) %>% 
  select(1,2,3,4,5,9,6,10,7,11,8,12)



fwrite(rdna_mouse_dataset_details_v2, 
       file = "rdna_mouse_2013_dataset_details_v2.csv")

##filter rows that contain a certain string

rdna_mouse_dataset_sequences<- rdna_mouse_dataset_details_v2 %>% select(1,2)


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


rdna_2013<- fread("rdna_mouse_2013_dataset_details_v2.csv", header = TRUE, sep = ",")
# select only nucleotide percent and identifier 
nucleotide<- rdna_2013 %>% select(Name,"A%","T%","G%","C%")
nucleotide_new<- separate(nucleotide, Name, "Name", sep = "_")


nucleotide_reshape <- pivot_longer(nucleotide_new, -Name, names_to = "Nucleotide", values_to = "Percent")
fwrite(nucleotide_reshape,  
       file= "nucleotide_rdna_mouse_2013_added_sequences_graph_input.csv")

#Please note below code do not work, because you can only drop one column 
#nar_reshape <- pivot_longer(nar, -Name, -Sequences,-Details, names_to = "Variable", values_to = "Value")

nucleotide_reshape$Name <- factor(nucleotide_reshape$Name, 
                                  levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2","28S", "3'ETS", "IGS" ))




#The fill aesthetic is used to color the bars based on the the X variable.

rdna_2013_sequences_nucleotide_distribution<- ggplot(nucleotide_reshape, aes(x = Name, y = Percent, fill = Nucleotide)) + 
  geom_bar(stat= 'identity', color = "black") +
  theme(axis.text.x = element_text(angle = 45, size = 10)) +
  labs(title= " Nucleotide distribution percent in mouse rDNA locus", 
       subtitle = "BK000964",
       x= "rDNA region", 
       y= "Nucleotide distribution percent")+
  geom_text(aes(label = Percent), position = position_stack(vjust = 0.5), size = 6)+ # Adjusted label positioning
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        theme(panel.grid = element_blank()),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25)) +  # Center Y-axis title
  theme(axis.ticks.y = element_line(color = "black"))

ggsave( "rdna_mouse_2013_sequences_nucleotide_distribution.tiff", 
        plot = rdna_2013_sequences_nucleotide_distribution, width=15,height=10, dpi=150)



#(python2.7) Downloads % python QmRLFS-finder.py -bed -i BK000964_added_5000nt_IGS_upstream_nontemplate.fasta -o BK000964_added_5000nt_IGS_upstream_nontemplate_qmrlfs
#QmRLFS-finder.py (version v1.5)
#run on Fri Apr 04 2025 13:14:27 
#command line: python QmRLFS-finder.py -bed -i BK000964_added_5000nt_IGS_upstream_nontemplate.fasta -o BK000964_added_5000nt_IGS_upstream_nontemplate_qmrlfs

#Time used: 0.36 mins



#(python2.7) Downloads % python QmRLFS-finder.py -i BK000964_added_5000nt_IGS_upstream_nontemplate.fasta -o BK000964_added_5000nt_IGS_upstream_nontemplate_qmrlfs 
#QmRLFS-finder.py (version v1.5)
#run on Fri Apr 04 2025 13:15:12 
#command line: python QmRLFS-finder.py -i BK000964_added_5000nt_IGS_upstream_nontemplate.fasta -o BK000964_added_5000nt_IGS_upstream_nontemplate_qmrlfs


# % conda activate python3.11
#(python3.11) % python g4_canonical_finder_3.11python.py BK000964_added_5000nt_IGS_upstream_nontemplate.fasta >output_BK000964_added_5000nt_IGS_upstream_nontemplate.txt
