# ------------------------------------------------------------------------------
#This code is part of paper: In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus.
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
#   This script calculates sequence identity between human and monkey rDNA
#   using global (Needleman–Wunsch) pairwise alignment. The analysis focuses
#   on the transcribed region (5′ETS–3′ETS) and its subcomponents to generate
#   an identity matrix summarizing cross-species conservation.
#
# Major Steps:
#   1. Load human (*KY962518*) and monkey (*KX061890*) rDNA FASTA sequences.
#   2. Extract the 5′ETS–3′ETS (coding) regions from both genomes.
#   3. Perform pairwise global alignment for each corresponding rDNA subregion.
#   4. Compute alignment statistics:
#        - Percent identity (pwalign::pid and NCBI-style)
#        - Match, mismatch counts, alignment score, and alignment length
#   5. Save all pairwise results and an identity matrix (human vs. monkey) to CSV files.
#
# Input:Fasta files for human and monkey entire and subregions. 
#
# Output:
#   - human_vs_monkey_rDNA_sequences_global_blast.csv
#       Pairwise alignment statistics between all human and monkey rDNA regions (Supplementary table 7 and Fig7 asterisk number)
#
#
# Notes:
#   - Alignments are computed using Biostrings::pairwiseAlignment (global mode)
#     with match = 2, mismatch = –3, gap opening = 5, gap extension = 2.
#   - Identity (%) is reported both via Biostrings (pid) and NCBI-style formula.
#   - Output matrices can be visualized as heatmaps for conservation profiling.
#
# Dependencies:
#   tidyverse
#   data.table
#   Biostrings
#   pwalign
#
# ------------------------------------------------------------------------------

#load libraries:

library(tidyverse)
library(data.table)
library(Biostrings)


#read required files for human

human_entire_rdna<- readDNAStringSet("KY962518_added_3500nt_IGS_upstream_nontemplate.fasta", format = "fasta", use.names = FALSE)
human_entire_rdna_seq<- as.character(human_entire_rdna[[1]])
nchar(human_entire_rdna_seq)
#48338


human_coding_rdna<- subseq(human_entire_rdna_seq, start = 3501, end= 16832)
nchar(human_coding_rdna)
#13332

human<- fread("rdna_hg38_chr21_2018_dataset_details_v3.csv", sep = ",", header = TRUE)
#this file has all rDNA region wise sequence and their nucleotide distribution
#9 rows


#read required files for monkey

monkey_rdna_KX061890<- readDNAStringSet("nontemplate_monkey_5ets_KX061890_3ets.fasta", format = "fasta", use.names = FALSE)
#default is fastq and use.names= FALSE will drop the header in the fasta file

monkey_coding_rdna_KX061890<- as.character(monkey_rdna_KX061890[[1]])
nchar(monkey_coding_rdna_KX061890)
#12979



monkey_rdna_NR146166<- readDNAStringSet("nontemplate_monkey_5ets_KX061890_and_NR_146166_3ets.fasta", format = "fasta", use.names = FALSE)
monkey_coding_rdna_NR146166<- as.character(monkey_rdna_NR146166[[1]])
nchar(monkey_coding_rdna_NR146166)
#12980


monkey<- fread("rdna_monkey_dataset_details_v1.csv", sep = ",", header = TRUE)
#this file has all rDNA region wise sequence of KX061890 and NR_146166 and their nucleotide distribution
#8 rows


#calculate pair wise alignment btw monkey and human
human_tmp<- human[c(2:8), c(1,2)]
new_row <- data.frame(Name = "human_coding_rdna", Sequences = human_coding_rdna)
human_tmp<- rbind(human_tmp, new_row)


monkey_tmp<- monkey[c(1:8), c(1,2)]
new_row1 <- data.frame(Name = "monkey_coding_rdna_KX061890", Sequences = monkey_coding_rdna_KX061890)
new_row2 <- data.frame(Name = "monkey_coding_rdna_NR146166", Sequences = monkey_coding_rdna_NR146166)
monkey_tmp<- rbind(monkey_tmp, new_row1, new_row2)


sigma<- pwalign::nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = TRUE) #for glopbal alignment

attempt2<- data.frame((matrix(nrow = 0, ncol=10)))
colnames(attempt2)<- c("query", "subject", "identity_perc", "NW_score_NCBI", "aligned_length", 
                       "match","identity_perc_NCBI","mismatch","aligned_query","aligned_subject")


for (i in 1:nrow(human_tmp)){
  # Define query sequence and identifier
  query_id<- human_tmp$Name[i]
  query_seq <- human_tmp$Sequences[i]
  
  #Loop through the rows for subject sequences
  
  for (j in 1:nrow(monkey_tmp)){
    #if (i !=j){ #would make sense if i dont want compare identical query and subject
    # Define subject sequence and identifier
    subject_id <- monkey_tmp$Name[j]
    subject_seq <- monkey_tmp$Sequences[j]
    
    alignment<- pwalign::pairwiseAlignment(query_seq, subject_seq,substitutionMatrix = sigma,type="global", gapOpening = 5, gapExtension = 2)
    
    query<-  query_id
    subject<-  subject_id
    identity_perc<- round(pwalign::pid(alignment),2) #identity percent
    NW_score_NCBI<- pwalign::score(alignment) #NW score
    aligned_length<- length(pwalign::alignedPattern(alignment)[[1]]) #The 2 objects returned by alignedPattern(x) and alignedSubject(x) are guaranteed to have the same shape (i.e. same length() and width())
    match<- pwalign::nmatch(alignment)
    identity_perc_NCBI<- round((match/aligned_length)*100, 2)#as per ncbi
    mismatch<- pwalign::nmismatch(alignment)
    aligned_query<-as.character(pwalign::alignedPattern(alignment)) #aligned query
    aligned_subject<- as.character(pwalign::alignedSubject(alignment)) #aligned subject
    
    
    
    final<- data.frame(query, subject, identity_perc, NW_score_NCBI, aligned_length, match,
                       identity_perc_NCBI,mismatch,
                       aligned_query,aligned_subject)
    
    attempt2 <- rbind(attempt2, final)
  }
  
}

fwrite(attempt2, "human_vs_monkey_rDNA_sequences_global_blast.csv")
#8*10 =80


globalblast_4c<- attempt2 %>% select(query, subject, identity_perc_NCBI) #4c - 4 columns
globalblast_4c$identifier <- paste(globalblast_4c$query, globalblast_4c$subject, sep = "&")

identity_matrix <- data.frame(matrix(NA, nrow=0, ncol=10)) #bcoz 10 rows in monkey_tmp

colnames(identity_matrix)<- t(monkey_tmp[,1]) #1 is Name

test<- data.frame(matrix(NA, nrow = 0, ncol = 1)) #additional column is for identifier column. This column will  contain query name
colnames(test)<- c( "identifier")

identity_matrix <- cbind(test, identity_matrix) #colnames now increased to 91.


for (i in 1:nrow(human_tmp)){
  listofvalues<- list()
  listofvalues <- append(listofvalues, human_tmp$Name[i])
  for (j in 1:length(monkey_tmp$Name)){
    k <- paste(human_tmp$Name[i], monkey_tmp$Name[j], sep="&") 
    match<- globalblast_4c[globalblast_4c$identifier== k,] #match is a dataframe that only contain one row which provided condition
    listofvalues <- append(listofvalues, match[1,3])
    #if (nrow(match)<1){ #code debug
    #print(k)
    #print(match)}#append function works by adding things row wise 
    #listofvalues [[i]][i]<- match[1,3] # this didnt work but good to know how to add data points to first entry of list
  }
  
  my_row <- do.call(rbind, listofvalues) # do.call is used to execute a function with a list of arguments.
  identity_matrix[nrow(identity_matrix)+1,] <- my_row #rbind will not work because my myrow gives as column 
  
} 

fwrite(identity_matrix,
       "human_vs_monkey_rDNA_sequences_global_blast_identity_matrix.csv")
