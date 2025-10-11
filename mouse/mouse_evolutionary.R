# ------------------------------------------------------------------------------
# This code is part of paper: In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus.
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj

# Purpose: This R code is designed to visualize R-loops, G4 and imotif forming sequences in mouse rDNA 
# the code takes the input data and makes plots all predicted non-canonical DNA sequences in mouse rDNA
# Inputs: G4FS_BK000964_added_5000nt_IGS_upstream_at_junctn_details.csv, RLFS_BK000964_added_5000nt_IGS_upstream_master_qmrlfs_table_after_rule.csv, 
#mouse_imotif_prediction_end_to_end_prediction_default_setting_master.csv (deposited in supplementary table)
# Outputs: PNG files of strand-specific coverage plots 
# ------------------------------------------------------------------------------


# load library
library(tidyverse)
library(data.table)
library(karyoploteR)

# load G-quadruplex forming sequences
entire_g4s_rdna <- fread("G4FS_BK000964_added_5000nt_IGS_upstream_at_junctn_details.csv", header = TRUE, sep = ",") #87
entire_g4s_rdna<- entire_g4s_rdna %>% select(chr, actual_G4FS_start, actual_G4FS_end, rDNA_region, strand)


# load R-loop forming sequences
rlfs<- fread("RLFS_BK000964_added_5000nt_IGS_upstream_master_qmrlfs_table_after_rule.csv", header = TRUE, sep = ",") #195
rlfs<- rlfs %>% select("#Name", actual_RLFS_start, actual_RLFS_end, rDNA_region, strand)

# load i-motif forming sequences
imotif<- fread("mouse_imotif_prediction_end_to_end_prediction_default_setting_master.csv", header = TRUE, sep = ",") #85
imotif<- imotif %>% select(chr, actual_imotif_start, actual_imotif_end, rDNA_region, strand, beg, end)


#prepare datasets:

datasets<- list(G4FS= entire_g4s_rdna,
                RLFS=rlfs)

nontemplate_filt_datasets <- list()  # initialize new list
template_filt_datasets <- list()

for (i in names(datasets)) {
  filt <- datasets[[i]] %>% 
    filter(!rDNA_region %in% c("Promoter", "IGS"))
  colnames(filt)<- c("chr", "actual_start", "actual_end", "rDNA_region", "strand")
  filt$name<- i
  filt$score<- 0
  filt<- filt %>% select(chr, actual_start, actual_end, name, score, strand)
  
  filt <- filt %>% mutate(actual_start= actual_start-5000)
  filt<- filt %>% mutate(actual_end = actual_end-5000)
  filt<- filt %>% mutate(start = ifelse(filt$strand == "+", actual_start, actual_end))
  filt<- filt %>% mutate(end = ifelse(filt$strand == "+", actual_end, actual_start))
  filt<- filt %>% select(chr, start, end, name, score, strand)
  
  nontemplate<- filt %>% filter(strand=="+")
  #because in NCBI keep nontemplate sequence.
  template<- filt %>% filter(strand=="-")
  
  
  nontemplate_filt_datasets[[paste0(i, "_filt")]] <- nontemplate
  template_filt_datasets[[paste0(i, "_filt")]] <- template
  
}


for (i in names(template_filt_datasets)){
  png(paste0("mouse_transcribed_rdna_template_", i , "_coverage.png"), width = 10, height = 10, units= "in", res = 600)
  
  
  custom_genome <- toGRanges(data.frame(chr="BK000964", start=1, end=13403))
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
  
  
  
  kpRect(kp, chr = 'BK000964', x0 = 1, x1 = 4007 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
  #5001+(4007-1) = 9007
  
  kpRect(kp, chr = 'BK000964', x0 = 4008, x1 = 5877, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
  #9008+(1870-1) = 10877
  
  kpRect(kp, chr = 'BK000964', x0 = 5878, x1 = 6877, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
  #10877+(1000-1) = 11877
  
  kpRect(kp, chr = 'BK000964', x0 = 6878, x1 = 7034 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
  #11878+(157-1) = 12034
  
  kpRect(kp, chr = 'BK000964', x0 = 7035, x1 = 8122, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
  #12035+(1088-1) = 13122
  
  kpRect(kp, chr = 'BK000964', x0 = 8123 , x1 = 12852, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
  #13123+(4730-1)= 17852
  
  kpRect(kp, chr = 'BK000964', x0 = 12853, x1 = 13403, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
  #17853+(551-1) = 18403
  
  kpPlotCoverage(kp, data=template_filt_datasets[[i]], col = "#1414E1", lwd = 6)
  dev.off()
  
}


for (i in names(nontemplate_filt_datasets)){
  png(paste0("mouse_transcribed_rdna_nontemplate_", i , "_coverage.png"), width = 10, height = 10, units= "in", res = 600)
  
  
  custom_genome <- toGRanges(data.frame(chr="BK000964", start=1, end=13403))
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
  
  
  kpRect(kp, chr = 'BK000964', x0 = 1, x1 = 4007 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
  #5001+(4007-1) = 9007
  
  kpRect(kp, chr = 'BK000964', x0 = 4008, x1 = 5877, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
  #9008+(1870-1) = 10877
  
  kpRect(kp, chr = 'BK000964', x0 = 5878, x1 = 6877, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
  #10877+(1000-1) = 11877
  
  kpRect(kp, chr = 'BK000964', x0 = 6878, x1 = 7034 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
  #11878+(157-1) = 12034
  
  kpRect(kp, chr = 'BK000964', x0 = 7035, x1 = 8122, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
  #12035+(1088-1) = 13122
  
  kpRect(kp, chr = 'BK000964', x0 = 8123 , x1 = 12852, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
  #13123+(4730-1)= 17852
  
  kpRect(kp, chr = 'BK000964', x0 = 12853, x1 = 13403, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
  #17853+(551-1) = 18403
  
  kpPlotCoverage(kp, data=nontemplate_filt_datasets[[i]], col = "#E21515", lwd = 6)
  dev.off()
  
}
  






#as for imotif i have not calculated promoter 
imotif$score<- 0
imotif_filt<- imotif %>% select(chr, beg, end,  rDNA_region, score,strand)
imotif_template<- imotif_filt %>% filter(strand=="-")
imotif_nontemplate<- imotif_filt %>% filter(strand=="+")



  png("mouse_transcribed_rdna_template_imotif_coverage.png", width = 10, height = 10, units= "in", res = 1000)
  
  
  custom_genome <- toGRanges(data.frame(chr="BK000964", start=1, end=13403))
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
  
  
  
  kpRect(kp, chr = 'BK000964', x0 = 1, x1 = 4007 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
  #5001+(4007-1) = 9007
  
  kpRect(kp, chr = 'BK000964', x0 = 4008, x1 = 5877, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
  #9008+(1870-1) = 10877
  
  kpRect(kp, chr = 'BK000964', x0 = 5878, x1 = 6877, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
  #10877+(1000-1) = 11877
  
  kpRect(kp, chr = 'BK000964', x0 = 6878, x1 = 7034 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
  #11878+(157-1) = 12034
  
  kpRect(kp, chr = 'BK000964', x0 = 7035, x1 = 8122, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
  #12035+(1088-1) = 13122
  
  kpRect(kp, chr = 'BK000964', x0 = 8123 , x1 = 12852, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
  #13123+(4730-1)= 17852
  
  kpRect(kp, chr = 'BK000964', x0 = 12853, x1 = 13403, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
  #17853+(551-1) = 18403
  
  kpPlotCoverage(kp, data=imotif_template, data.panel = 2, col ="#1414E1",  r0= -0.5, r1= -0.7,lwd = 6)
  dev.off()
  



  png("mouse_transcribed_rdna_nontemplate_imotif_coverage.png", width = 10, height = 10, units= "in", res = 1000)
  
  
  custom_genome <- toGRanges(data.frame(chr="BK000964", start=1, end=13403))
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
  
  
  kpRect(kp, chr = 'BK000964', x0 = 1, x1 = 4007 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
  #5001+(4007-1) = 9007
  
  kpRect(kp, chr = 'BK000964', x0 = 4008, x1 = 5877, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
  #9008+(1870-1) = 10877
  
  kpRect(kp, chr = 'BK000964', x0 = 5878, x1 = 6877, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
  #10877+(1000-1) = 11877
  
  kpRect(kp, chr = 'BK000964', x0 = 6878, x1 = 7034 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
  #11878+(157-1) = 12034
  
  kpRect(kp, chr = 'BK000964', x0 = 7035, x1 = 8122, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
  #12035+(1088-1) = 13122
  
  kpRect(kp, chr = 'BK000964', x0 = 8123 , x1 = 12852, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
  #13123+(4730-1)= 17852
  
  kpRect(kp, chr = 'BK000964', x0 = 12853, x1 = 13403, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
  #17853+(551-1) = 18403
  
  kpPlotCoverage(kp, data=imotif_nontemplate, data.panel = 2, col = "#E21515",  r0= -0.5, r1= -0.7, lwd = 6)
  dev.off()
  