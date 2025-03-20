

#Task is to create GC skew for human RDNA locus 
# and make coorelation with RLFS.

#first make GC for each region and then calculate correlation with RLFS in that region.


#next, I wanted to take average RLFS length and then calculate GC skew in them.

#load GC skew function
source("/Users/jyotiadala/Downloads/gc_skew_function.R")

setwd("/Users/jyotiadala/Downloads/")

library(stringr)
library(Biostrings)
library(tidyverse)
library(data.table)


rdna_human<- readDNAStringSet("KY962518_added_3500nt_IGS_upstream_nontemplate.fasta")
seq<- as.character(rdna_human[[1]])

gc_skew_data <- gc_skew(seq, window_size = 10)
fwrite(gc_skew_data, "gc_skew_data.csv", sep=",")


ggplot(gc_skew_data, aes(x = start, y = GC_skew_value)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "GC Skew Across DNA Sequence",
       x = "Position (bp)",
       y = "GC Skew") +
  theme_minimal()
