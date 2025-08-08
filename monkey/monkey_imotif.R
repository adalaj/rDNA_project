#monkey imotif

# go to https://im-seeker.org/

#install libraries
library(data.table)
library(tidyverse)
library(Biostrings)
library(karyoploteR)


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/monkey")

monkey_rDNA<- readDNAStringSet(file = "nontemplate_monkey_5ets_KX065350_3ets.fasta") #belongs to package BIOSTRINGS
monkey_rDNA_seq<- monkey_rDNA[[1]]
nchar(monkey_rDNA_seq)
#12979

template_monkey_rdna<- reverseComplement(monkey_rDNA_seq)
nchar(template_monkey_rdna)
#12979

template_monkey_rdna_set <- DNAStringSet(template_monkey_rdna)

names(template_monkey_rdna_set) <- "template_KX065350"
writeXStringSet(template_monkey_rdna_set, "template_monkey_5ets_KX065350_3ets.fasta")


monkey_modified_rDNA<- readDNAStringSet(file = "nontemplate_monkey_5ets_KX065350_and_NR_146166_3ets.fasta") #belongs to package BIOSTRINGS
monkey_modified_rDNA_seq<- monkey_modified_rDNA[[1]]
nchar(monkey_modified_rDNA_seq)
#12980

template_monkey_modified_rdna<- reverseComplement(monkey_modified_rDNA_seq)
nchar(template_monkey_modified_rdna)
#12980

template_monkey_modified_rdna_set <- DNAStringSet(template_monkey_modified_rdna)

names(template_monkey_modified_rdna_set) <- "template_KX065350_and_NR_146166"
writeXStringSet(template_monkey_modified_rdna_set, "template_monkey_5ets_KX065350_and_NR_146166_3ets.fasta")



# I added KX065350_5ets_to_3ets_nontemplate.fasta to predict>upload DNA FILE (FASTA FORMAT) (https://im-seeker.org/predict)>continue>end to end with default model>predict
# I got a csv file which i saved as monkey_imotif_prediction_end_to_end_default_setting_nontemplate.csv



setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/imotif/monkey")

datasets<- list(
  KX065350 = "KX065350",
  KX065350_and_NR_146166 = "KX065350_and_NR_146166"
)





for (i in names(datasets)){
filename <- paste0("monkey_imotif_5ets_",i,"_3ets_prediction_end_to_end_default_setting_template.csv")
template<- fread(filename, sep = ",", header = TRUE)

template<- template %>%
  unite(imotif_sequence, c(ctract1, loop1, ctract2, loop2, ctract3, loop3, ctract4), sep = "", remove = FALSE)
template$strand<- "-"
template$chr = i



filename2 <- paste0("monkey_imotif_5ets_",i,"_3ets_prediction_end_to_end_default_setting_nontemplate.csv")
nontemplate<- fread(filename2, sep = ",", header = TRUE)
nontemplate<- nontemplate %>%
  unite(imotif_sequence, c(ctract1, loop1, ctract2, loop2, ctract3, loop3, ctract4), sep = "", remove = FALSE)
nontemplate$strand <- "+"
nontemplate$chr = i

fwrite(template, paste0("monkey_imotif_5ets_", i, "_3ets_prediction_end_to_end_prediction_default_setting_template_modified.csv"))
fwrite(nontemplate, paste0("monkey_imotif_5ets_", i, "_3ets_prediction_end_to_end_prediction_default_setting_nontemplate_modified.csv"))
master<- rbind(nontemplate, template)
master<- master %>% mutate(length= end-beg)
master<- master %>% mutate(actual_imotif_start = ifelse(master$strand == "+", beg, end))
master<- master %>% mutate(actual_imotif_end = ifelse(master$strand=="+", end, beg))



master$rDNA_region <- "junction"

master$rDNA_region[master$actual_imotif_start >= 1] <- "5'ETS"

master$rDNA_region[master$actual_imotif_start > 3640 ] <- "5'ETS and 18S junction"

master$rDNA_region[master$actual_imotif_start >= 3641] <- "18S"
master$rDNA_region[master$actual_imotif_start > 5508 ] <- "18S and ITS1 junction"


master$rDNA_region[master$actual_imotif_start >= 5509   ] <- "ITS1"
master$rDNA_region[master$actual_imotif_start > 6535] <- "ITS1 and 5.8S junction"

master$rDNA_region[master$actual_imotif_start >= 6536 ] <- "5.8S"
master$rDNA_region[master$actual_imotif_start > 6692 ] <- "5.8S and ITS2 junction"

master$rDNA_region[master$actual_imotif_start >= 6693] <- "ITS2"
master$rDNA_region[master$actual_imotif_start > 7863] <- "ITS2 and 28S junction"

master$rDNA_region[master$actual_imotif_start >= 7864 ] <- "28S"
master$rDNA_region[master$actual_imotif_start > 12648 ] <- "28S and 3'ETS junction"

master$rDNA_region[master$actual_imotif_start >= 12649] <- "3'ETS"

fwrite(master, paste0("monkey_imotif_5ets_",i,"_3ets_prediction_end_to_end_prediction_default_setting_master.csv"))

}

