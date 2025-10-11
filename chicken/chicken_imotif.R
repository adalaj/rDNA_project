# ------------------------------------------------------------------------------
# This code is part of paper: In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus.
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
# This R script processes i-motif forming sequence (iMFS) predictions in the 
# chicken (Gallus gallus) rDNA locus to assign each predicted iMFS to defined 
# rDNA subregions (5′ETS, 18S, ITS1, 5.8S, ITS2, 28S, 3′ETS, and IGS).
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
#   - Chicken rDNA FASTA sequence (GenBank: KT445934)
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
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/chicken/input")

chicken_rDNA<- readDNAStringSet(file = "KT445934_chicken_rDNA_2017.fasta") #JA:belongs to package BIOSTRINGS
chicken_rDNA_seq<- chicken_rDNA[[1]]
nchar(chicken_rDNA_seq)
#JA:11863


#JA:made template because iM seeker do not calculate for iMFS in opposite (template) strand
template_chicken_rdna<- reverseComplement(chicken_rDNA_seq)
nchar(template_chicken_rdna)
#JA:11863

template_chicken_rdna_set <- DNAStringSet(template_chicken_rdna)

names(template_chicken_rdna_set) <- "template_KT445934"
writeXStringSet(template_chicken_rdna_set, "KT445934_chicken_rDNA_2017_template.fasta")


#JA:steps i did to find IMFS
#JA: go to https://im-seeker.org/
#JA: I added KT445934_5ets_to_3ets_nontemplate.fasta to predict>upload DNA FILE (FASTA FORMAT) (https://im-seeker.org/predict)>continue>end to end with default model>predict
#JA: I got a csv file which i saved as chicken_imotif_prediction_end_to_end_default_setting_nontemplate.csv
#JA: did the same for temaplate and output was saved as chicken_imotif_prediction_end_to_end_prediction_default_setting_template.csv



#JA:combined both template and nontemplate iMFS in one file for downstream analysis
template<- fread("chicken_imotif_prediction_end_to_end_prediction_default_setting_template.csv", sep = ",", header = TRUE)
template<- template %>%
  unite(imotif_sequence, c(ctract1, loop1, ctract2, loop2, ctract3, loop3, ctract4), sep = "", remove = FALSE)
template$strand<- "-"
template$chr = "KT445934"

nontemplate<- fread("chicken_imotif_prediction_end_to_end_prediction_default_setting_nontemplate.csv", sep = ",", header = TRUE)
nontemplate<- nontemplate %>%
  unite(imotif_sequence, c(ctract1, loop1, ctract2, loop2, ctract3, loop3, ctract4), sep = "", remove = FALSE)
nontemplate$strand <- "+"
nontemplate$chr = "KT445934"

fwrite(template, "chicken_imotif_prediction_end_to_end_prediction_default_setting_template_modified.csv")
fwrite(nontemplate, "chicken_imotif_prediction_end_to_end_prediction_default_setting_nontemplate_modified.csv")

master<- rbind(nontemplate, template)
master<- master %>% mutate(length= end-beg)
master<- master %>% mutate(actual_imotif_start = ifelse(master$strand == "+", beg, end))
master<- master %>% mutate(actual_imotif_end = ifelse(master$strand=="+", end, beg))



master$rDNA_region <- "junction"

master$rDNA_region[master$actual_imotif_start >= 1] <- "5'ETS"

master$rDNA_region[master$actual_imotif_start > 1836 ] <- "5'ETS and 18S junction"

master$rDNA_region[master$actual_imotif_start >= 1837] <- "18S"
master$rDNA_region[master$actual_imotif_start > 3659 ] <- "18S and ITS1 junction"


master$rDNA_region[master$actual_imotif_start >= 3660   ] <- "ITS1"
master$rDNA_region[master$actual_imotif_start > 6189] <- "ITS1 and 5.8S junction"

master$rDNA_region[master$actual_imotif_start >= 6190 ] <- "5.8S"
master$rDNA_region[master$actual_imotif_start > 6346 ] <- "5.8S and ITS2 junction"

master$rDNA_region[master$actual_imotif_start >= 6347] <- "ITS2"
master$rDNA_region[master$actual_imotif_start > 7079] <- "ITS2 and 28S junction"

master$rDNA_region[master$actual_imotif_start >= 7080 ] <- "28S"
master$rDNA_region[master$actual_imotif_start > 11520 ] <- "28S and 3'ETS junction"

master$rDNA_region[master$actual_imotif_start >= 11521] <- "3'ETS"
master$rDNA_region[master$actual_imotif_start > 11863] <- "3'ETS and IGS junction"

fwrite(master, "chicken_imotif_prediction_end_to_end_prediction_default_setting_master.csv")



