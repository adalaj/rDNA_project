# ------------------------------------------------------------------------------
# This code is part of paper: In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus.
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
# This R script processes i-motif forming sequence (iMFS) predictions in the 
# mouse rDNA locus to assign each predicted iMFS to defined 
# rDNA subregions (promoter, 5′ETS, 18S, ITS1, 5.8S, ITS2, 28S, 3′ETS, and IGS).
#
# Specifically, it:
#   1. Generates a reverse-complement (template) strand FASTA sequence, since 
#      iM-seeker analyzes only the input (non-template) strand.
#   2. Merges iM-seeker output files from both strands.
#   3. Computes iMFS lengths and normalized start–end coordinates.
#   4. Annotates each iMFS with its corresponding rDNA subregion.
#
# Inputs:
#   - iM-seeker output CSV files for non-template and template strands
#     (https://im-seeker.org/)
#   - mouse rDNA FASTA sequence (GenBank: BK000964)
#
# Outputs:
#   - Annotated CSV file containing all predicted iMFSs with strand, 
#     coordinates, length, and assigned rDNA region.
#
# ------------------------------------------------------------------------------


#install libraries
library(data.table)
library(tidyverse)
library(Biostrings)
library(karyoploteR)

mouse_rDNA<- readDNAStringSet(file = "BK000964_added_5000nt_IGS_upstream_nontemplate.fasta") #belongs to package BIOSTRINGS
mouse_rDNA_seq<- mouse_rDNA[[1]]
nchar(mouse_rDNA_seq)
# 50306

template_mouse_rdna<- reverseComplement(mouse_rDNA_seq)
nchar(template_mouse_rdna)
# 50306

template_mouse_rdna_set <- DNAStringSet(template_mouse_rdna)

names(template_mouse_rdna_set) <- "template_BK000964"
writeXStringSet(template_mouse_rdna_set, "BK000964_added_3500nt_IGS_upstream_template.fasta")


#at this moment i only need 3'ETS to 3'ets

#So i am gonna subseq 
mouse_5ets_to_3ets <- subseq(mouse_rDNA_seq, start =5001, end = 18403)
nchar(mouse_5ets_to_3ets)
#13403
mouse_5ets_to_3ets_set <- DNAStringSet(mouse_5ets_to_3ets)


names(mouse_5ets_to_3ets_set) <- "nontemplate_BK000964"
writeXStringSet(mouse_5ets_to_3ets_set, "BK000964_5ets_to_3ets_nontemplate.fasta")


template_mouse_5ets_to_3ets<- reverseComplement(mouse_5ets_to_3ets)
nchar(template_mouse_5ets_to_3ets)
# 13403

template_mouse_5ets_to_3ets_set <- DNAStringSet(template_mouse_5ets_to_3ets)

names(template_mouse_5ets_to_3ets_set) <- "template_BK000964"
writeXStringSet(template_mouse_5ets_to_3ets_set, "BK000964_5ets_to_3ets_template.fasta")



# I added BK000964_5ets_to_3ets_nontemplate.fasta to predict>upload DNA FILE (FASTA FORMAT) (https://im-seeker.org/predict)>continue>end to end with default model>predict
# I got a csv file which i saved as mouse_imotif_prediction_end_to_end_default_setting_nontemplate.csv



setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/imotif/mouse")
template<- fread("mouse_imotif_prediction_end_to_end_prediction_default_setting_template.csv", sep = ",", header = TRUE)

template<- template %>%
  unite(imotif_sequence, c(ctract1, loop1, ctract2, loop2, ctract3, loop3, ctract4), sep = "", remove = FALSE)
template$strand<- "-"
template$chr = "BK000964"

nontemplate<- fread("mouse_imotif_prediction_end_to_end_prediction_default_setting_nontemplate.csv", sep = ",", header = TRUE)
nontemplate<- nontemplate %>%
  unite(imotif_sequence, c(ctract1, loop1, ctract2, loop2, ctract3, loop3, ctract4), sep = "", remove = FALSE)
nontemplate$strand <- "+"
nontemplate$chr = "BK000964"

fwrite(template, "mouse_imotif_prediction_end_to_end_prediction_default_setting_template_modified.csv")
fwrite(nontemplate, "mouse_imotif_prediction_end_to_end_prediction_default_setting_nontemplate_modified.csv")

master<- rbind(nontemplate, template)
master<- master %>% mutate(length= end-beg)
master<- master %>% mutate(actual_imotif_start = ifelse(master$strand == "+", beg, end))
master<- master %>% mutate(actual_imotif_end = ifelse(master$strand=="+", end, beg))


fwrite(master, "mouse_imotif_prediction_end_to_end_prediction_default_setting_master.csv")

master$rDNA_region <- "junction"

master$rDNA_region[master$actual_imotif_start >= 1] <- "5'ETS"

master$rDNA_region[master$actual_imotif_start > 4007 ] <- "5'ETS and 18S junction"

master$rDNA_region[master$actual_imotif_start >= 4008] <- "18S"
master$rDNA_region[master$actual_imotif_start > 5877 ] <- "18S and ITS1 junction"


master$rDNA_region[master$actual_imotif_start >= 5878   ] <- "ITS1"
master$rDNA_region[master$actual_imotif_start > 6877] <- "ITS1 and 5.8S junction"

master$rDNA_region[master$actual_imotif_start >= 6878 ] <- "5.8S"
master$rDNA_region[master$actual_imotif_start > 7034 ] <- "5.8S and ITS2 junction"

master$rDNA_region[master$actual_imotif_start >= 7035] <- "ITS2"
master$rDNA_region[master$actual_imotif_start > 8122] <- "ITS2 and 28S junction"

master$rDNA_region[master$actual_imotif_start >= 8123 ] <- "28S"
master$rDNA_region[master$actual_imotif_start > 12852 ] <- "28S and 3'ETS junction"

master$rDNA_region[master$actual_imotif_start >= 12853] <- "3'ETS"
master$rDNA_region[master$actual_imotif_start > 13403 ] <- "3'ETS and IGS junction"

fwrite(master, "mouse_imotif_prediction_end_to_end_prediction_default_setting_master.csv")



