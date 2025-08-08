#!/usr/bin/env Rscript


#load libraries:

library(tidyverse)
library(data.table)
library(Biostrings)

#read required files for human

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/input")
human_entire_rdna<- readDNAStringSet("KY962518_added_3500nt_IGS_upstream_nontemplate.fasta", format = "fasta", use.names = FALSE)
human_entire_rdna_seq<- as.character(human_entire_rdna[[1]])
nchar(human_entire_rdna_seq)
#48338


human_coding_rdna<- subseq(human_entire_rdna_seq, start = 3501, end= 16832)
nchar(human_coding_rdna)
#13332

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/output")
human<- fread("rdna_hg38_chr21_2018_dataset_details_v3.csv", sep = ",", header = TRUE)
#this file has all rDNA region wise sequence and their nucleotide distribution
#9 rows


#read required files for mouse
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/mouse/input")
mouse_entire_rdna<- readDNAStringSet("BK000964_added_5000nt_IGS_upstream_nontemplate.fasta", format = "fasta", use.names = FALSE)
mouse_entire_rdna_seq<- as.character(mouse_entire_rdna[[1]])
nchar(mouse_entire_rdna_seq)
#50306


mouse_coding_rdna<- subseq(mouse_entire_rdna_seq, start = 5001, end= 18403)
nchar(mouse_coding_rdna)
#13403

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/mouse/output")
mouse<- fread("rdna_mouse_2013_dataset_details_v2.csv", sep = ",", header = TRUE)
#this file has all rDNA region wise sequence and their nucleotide distribution
#9 rows


#calculate pair wise alignment btw mouse and human
human_tmp<- human[c(2:8), c(1,2)]
new_row <- data.frame(Name = "human_coding_rdna", Sequences = human_coding_rdna)
human_tmp<- rbind(human_tmp, new_row)


mouse_tmp<- mouse[c(2:8), c(1,2)]
new_row1 <- data.frame(Name = "mouse_coding_rdna", Sequences = mouse_coding_rdna)
mouse_tmp<- rbind(mouse_tmp, new_row1)


sigma<- pwalign::nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = TRUE)
#sigma_local<- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -2, baseOnly = TRUE)

attempt2<- data.frame((matrix(nrow = 0, ncol=10)))
colnames(attempt2)<- c("query", "subject", "identity_perc", "NW_score_NCBI", "aligned_length", 
                       "match","identity_perc_NCBI","mismatch","aligned_query","aligned_subject")


for (i in 1:nrow(human_tmp)){
  # Define query sequence and identifier
  query_id<- human_tmp$Name[i]
  query_seq <- human_tmp$Sequences[i]
  
  #Loop through the rows for subject sequences
  
  for (j in 1:nrow(mouse_tmp)){
    #if (i !=j){ #would make sense if i dont want compare identical query and subject
    # Define subject sequence and identifier
    subject_id <- mouse_tmp$Name[j]
    subject_seq <- mouse_tmp$Sequences[j]
    
    alignment<- pwalign::pairwiseAlignment(query_seq, subject_seq,substitutionMatrix = sigma,type="global", gapOpening = 5, gapExtension = 2)
    #alignment<- pwalign::pairwiseAlignment(query_seq, subject_seq,substitutionMatrix = sigma_local,type="local", gapOpening = 5, gapExtension = 2)
    
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

fwrite(attempt2, "human_vs_mouse_rDNA_sequences_global_blast.csv")
#8*8 =64


globalblast_4c<- attempt2 %>% select(query, subject, identity_perc_NCBI) #4c - 4 columns
globalblast_4c$identifier <- paste(globalblast_4c$query, globalblast_4c$subject, sep = "&")

identity_matrix <- data.frame(matrix(NA, nrow=0, ncol=8)) #bcoz 8 rows in human_tmp

colnames(identity_matrix)<- t(mouse_tmp[,1]) #1 is Name

test<- data.frame(matrix(NA, nrow = 0, ncol = 1)) #additional column is for identifier column. This column will contain query name
colnames(test)<- c( "identifier")

identity_matrix <- cbind(test, identity_matrix) #colnames now increased to 91.


for (i in 1:nrow(human_tmp)){
  listofvalues<- list()
  listofvalues <- append(listofvalues, human_tmp$Name[i])
  for (j in 1:length(mouse_tmp$Name)){
    k <- paste(human_tmp$Name[i], mouse_tmp$Name[j], sep="&") 
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
       "human_vs_mouse_rDNA_sequences_global_blast_identity_matrix.csv")
