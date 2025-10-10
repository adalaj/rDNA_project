#i wanted to show evolutionary conservation
#everything will be 5'ETS to 3'ETS
# thinking to write a code that will code and simultaneously save the graphs
#visualisation is identical for KX061890 and KX061890_and_NR_146166

#for KX061890

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/g4s_and_rdna/monkey")
entire_g4s_rdna <- fread("pG4CS_KX061890_and_NR_146166_monkey_junctn_details.csv", header = TRUE, sep = ",")#114
entire_g4s_rdna<- entire_g4s_rdna %>% select(chr, actual_pG4CS_start, actual_pG4CS_end, rDNA_region, strand)

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/monkey")
rlfs<- fread("RLFS_KX061890_monkey_junctn_details.csv", header = TRUE, sep = ",")#131
rlfs<- rlfs %>% select(GenBank_Accession, actual_rlfs_start, actual_rlfs_end, rDNA_region, strand)


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/imotif/monkey")
imotif<- fread("monkey_imotif_5ets_KX061890_3ets_prediction_end_to_end_prediction_default_setting_master.csv", header = TRUE, sep = ",")#53
imotif<- imotif %>% select(chr, actual_imotif_start, actual_imotif_end, rDNA_region, strand)


#prepare datasets:

datasets<- list(pG4CS=entire_g4s_rdna,
                RLFS=rlfs,
                imotif=imotif)

nontemplate_filt_datasets <- list()  # initialize new list
template_filt_datasets <- list()

for (i in names(datasets)) {
  filt <- datasets[[i]] 
  colnames(filt)<- c("chr", "actual_start", "actual_end", "rDNA_region", "strand")
  filt$chr<- "KX061890"
  filt$name<- i
  filt$score<- 0
  filt<- filt %>% select(chr, actual_start, actual_end, name, score, strand)
  filt<- filt %>% mutate(start = ifelse(filt$strand == "+", actual_start, actual_end))
  filt<- filt %>% mutate(end = ifelse(filt$strand == "+", actual_end, actual_start))
  filt<- filt %>% select(chr, start, end, name, score, strand)
  nontemplate<- filt %>% filter(strand=="+")
  #because in NCBI keep nontemplate sequence.
  template<- filt %>% filter(strand=="-")
  
  
  nontemplate_filt_datasets[[paste0(i, "_filt")]] <- nontemplate
  template_filt_datasets[[paste0(i, "_filt")]] <- template
  
}


for (i in names(template_filt_datasets)[1:2]){
  png(paste0("monkey_transcribed_rdna_template_KX061890_", i , "_coverage.png"), width = 10, height = 10, units= "in", res = 600)
  
  
  custom_genome <- toGRanges(data.frame(chr="KX061890", start=1, end=12979))
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
  
  
  
  kpRect(kp, chr = 'KX061890', x0 = 1, x1 = 3640 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
  #5001+(3640-1) = 9007
  
  kpRect(kp, chr = 'KX061890', x0 = 3641, x1 = 5508, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
  #9008+(1870-1) = 10877
  
  kpRect(kp, chr = 'KX061890', x0 = 5509, x1 = 6535, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
  #10877+(1000-1) = 11877
  
  kpRect(kp, chr = 'KX061890', x0 = 6536, x1 = 6692 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
  #11878+(157-1) = 12034
  
  kpRect(kp, chr = 'KX061890', x0 = 6693, x1 = 7863, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
  #12035+(1088-1) = 13122
  
  kpRect(kp, chr = 'KX061890', x0 = 7864 , x1 = 12648, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
  #13123+(4730-1)= 17852
  
  kpRect(kp, chr = 'KX061890', x0 = 12649, x1 = 12979, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
  #17853+(551-1) = 18403
  
  kpPlotCoverage(kp, data=template_filt_datasets[[i]], col = "#1414E1",lwd = 6)
  dev.off()
  
}


for (i in names(nontemplate_filt_datasets)[1:2]){
  png(paste0("monkey_transcribed_rdna_nontemplate_KX061890_", i , "_coverage.png"), width = 10, height = 10, units= "in", res = 600)
  
  
  custom_genome <- toGRanges(data.frame(chr="KX061890", start=1, end=12979))
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
  
  
  kpRect(kp, chr = 'KX061890', x0 = 1, x1 = 3640 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
  #5001+(3640-1) = 9007
  
  kpRect(kp, chr = 'KX061890', x0 = 3641, x1 = 5508, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
  #9008+(1870-1) = 10877
  
  kpRect(kp, chr = 'KX061890', x0 = 5509, x1 = 6535, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
  #10877+(1000-1) = 11877
  
  kpRect(kp, chr = 'KX061890', x0 = 6536, x1 = 6692 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
  #11878+(157-1) = 12034
  
  kpRect(kp, chr = 'KX061890', x0 = 6693, x1 = 7863, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
  #12035+(1088-1) = 13122
  
  kpRect(kp, chr = 'KX061890', x0 = 7864 , x1 = 12648, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
  #13123+(4730-1)= 17852
  
  kpRect(kp, chr = 'KX061890', x0 = 12649, x1 = 12979, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
  #17853+(551-1) = 18403
  
  kpPlotCoverage(kp, data=nontemplate_filt_datasets[[i]], col = "#E21515", lwd = 6)
  dev.off()
  
}


for (i in names(template_filt_datasets)[3]){
  png(paste0("monkey_transcribed_rdna_template_KX061890_", i , "_coverage.png"), width = 10, height = 10, units= "in", res = 1000)
  
  
  custom_genome <- toGRanges(data.frame(chr="KX061890", start=1, end=12979))
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
  
  
  
  kpRect(kp, chr = 'KX061890', x0 = 1, x1 = 3640 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
  #5001+(3640-1) = 9007
  
  kpRect(kp, chr = 'KX061890', x0 = 3641, x1 = 5508, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
  #9008+(1870-1) = 10877
  
  kpRect(kp, chr = 'KX061890', x0 = 5509, x1 = 6535, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
  #10877+(1000-1) = 11877
  
  kpRect(kp, chr = 'KX061890', x0 = 6536, x1 = 6692 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
  #11878+(157-1) = 12034
  
  kpRect(kp, chr = 'KX061890', x0 = 6693, x1 = 7863, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
  #12035+(1088-1) = 13122
  
  kpRect(kp, chr = 'KX061890', x0 = 7864 , x1 = 12648, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
  #13123+(4730-1)= 17852
  
  kpRect(kp, chr = 'KX061890', x0 = 12649, x1 = 12979, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
  #17853+(551-1) = 18403
  
  kpPlotCoverage(kp, data=template_filt_datasets[[i]],  data.panel = 2, col = "#1414E1", r0= -0.5, r1= -0.7,  lwd = 6)
  dev.off()
  
}


for (i in names(nontemplate_filt_datasets)[3]){
  png(paste0("monkey_transcribed_rdna_nontemplate_KX061890_", i , "_coverage.png"), width = 10, height = 10, units= "in", res = 1000)
  
  
  custom_genome <- toGRanges(data.frame(chr="KX061890", start=1, end=12979))
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
  
  
  kpRect(kp, chr = 'KX061890', x0 = 1, x1 = 3640 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
  #5001+(3640-1) = 9007
  
  kpRect(kp, chr = 'KX061890', x0 = 3641, x1 = 5508, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
  #9008+(1870-1) = 10877
  
  kpRect(kp, chr = 'KX061890', x0 = 5509, x1 = 6535, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
  #10877+(1000-1) = 11877
  
  kpRect(kp, chr = 'KX061890', x0 = 6536, x1 = 6692 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
  #11878+(157-1) = 12034
  
  kpRect(kp, chr = 'KX061890', x0 = 6693, x1 = 7863, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
  #12035+(1088-1) = 13122
  
  kpRect(kp, chr = 'KX061890', x0 = 7864 , x1 = 12648, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
  #13123+(4730-1)= 17852
  
  kpRect(kp, chr = 'KX061890', x0 = 12649, x1 = 12979, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
  #17853+(551-1) = 18403
  
  kpPlotCoverage(kp, data=nontemplate_filt_datasets[[i]], data.panel = 2, col = "#E21515",  r0= -0.5, r1= -0.7,  lwd = 5)
  dev.off()
  
}




# make new for KX061890_and_NR_146166, custom genome to 12980

  setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/g4s_and_rdna/monkey")
  entire_g4s_rdna <- fread("pG4CS_KX061890_and_NR_146166_monkey_junctn_details.csv", header = TRUE, sep = ",")
  entire_g4s_rdna<- entire_g4s_rdna %>% select(chr, actual_pG4CS_start, actual_pG4CS_end, rDNA_region, strand)
  
  setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/monkey")
  rlfs<- fread("RLFS_KX061890_and_NR_146166_monkey_junctn_details.csv", header = TRUE, sep = ",")
  rlfs<- rlfs %>% select(GenBank_Accession, actual_rlfs_start, actual_rlfs_end, rDNA_region, strand)
  
  
  setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/imotif/monkey")
  imotif<- fread("monkey_imotif_5ets_KX061890_and_NR_146166_3ets_prediction_end_to_end_prediction_default_setting_master.csv", header = TRUE, sep = ",")
  imotif<- imotif %>% select(chr, actual_imotif_start, actual_imotif_end, rDNA_region, strand)
  
  
  #prepare datasets:
  
  datasets<- list(pG4CS=entire_g4s_rdna,
                  RLFS=rlfs,
                  imotif=imotif)
  
  nontemplate_filt_datasets <- list()  # initialize new list
  template_filt_datasets <- list()
  
  for (i in names(datasets)) {
    filt <- datasets[[i]] 
      colnames(filt)<- c("chr", "actual_start", "actual_end", "rDNA_region", "strand")
      filt$chr<- "KX061890_and_NR_146166"
      filt$name<- i
      filt$score<- 0
      filt<- filt %>% select(chr, actual_start, actual_end, name, score, strand)
      filt<- filt %>% mutate(start = ifelse(filt$strand == "+", actual_start, actual_end))
      filt<- filt %>% mutate(end = ifelse(filt$strand == "+", actual_end, actual_start))
      filt<- filt %>% select(chr, start, end, name, score, strand)
      nontemplate<- filt %>% filter(strand=="+")
      #because in NCBI keep nontemplate sequence.
      template<- filt %>% filter(strand=="-")
      
      
      nontemplate_filt_datasets[[paste0(i, "_filt")]] <- nontemplate
      template_filt_datasets[[paste0(i, "_filt")]] <- template
      
  }
  
  
  for (i in names(template_filt_datasets)[1:2]){
    png(paste0("monkey_transcribed_rdna_template_KX061890_and_NR_146166_", i , "_coverage.png"), width = 10, height = 10, units= "in", res = 200)
    
    
    custom_genome <- toGRanges(data.frame(chr="KX061890_and_NR_146166", start=1, end=12980))
    kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
    
    
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 1, x1 = 3640 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
    #5001+(3640-1) = 9007
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 3641, x1 = 5508, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
    #9008+(1870-1) = 10877
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 5509, x1 = 6535, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
    #10877+(1000-1) = 11877
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 6536, x1 = 6692 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
    #11878+(157-1) = 12034
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 6693, x1 = 7863, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
    #12035+(1088-1) = 13122
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 7864 , x1 = 12648, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
    #13123+(4730-1)= 17852
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 12649, x1 = 12980, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
    #17853+(551-1) = 18403
    
    kpPlotCoverage(kp, data=template_filt_datasets[[i]], col = "#1414E1", lwd = 6)
    dev.off()
    
  }
  
  
  for (i in names(nontemplate_filt_datasets)[1:2]){
    png(paste0("monkey_transcribed_rdna_nontemplate_KX061890_and_NR_146166_", i , "_coverage.png"), width = 10, height = 10, units= "in", res = 200)
    
    
    custom_genome <- toGRanges(data.frame(chr="KX061890_and_NR_146166", start=1, end=12980))
    kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
    
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 1, x1 = 3640 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
    #5001+(3640-1) = 9007
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 3641, x1 = 5508, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
    #9008+(1870-1) = 10877
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 5509, x1 = 6535, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
    #10877+(1000-1) = 11877
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 6536, x1 = 6692 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
    #11878+(157-1) = 12034
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 6693, x1 = 7863, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
    #12035+(1088-1) = 13122
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 7864 , x1 = 12648, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
    #13123+(4730-1)= 17852
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 12649, x1 = 12980, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
    #17853+(551-1) = 18403
    
    kpPlotCoverage(kp, data=nontemplate_filt_datasets[[i]], col = "#E21515", lwd = 6)
    dev.off()
    
  }
  
  for (i in names(template_filt_datasets)[3]){
    png(paste0("monkey_transcribed_rdna_template_KX061890_and_NR_146166_", i , "_coverage.png"), width = 10, height = 10, units= "in", res = 200)
    
    
    custom_genome <- toGRanges(data.frame(chr="KX061890_and_NR_146166", start=1, end=12980))
    kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
    
    
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 1, x1 = 3640 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
    #5001+(3640-1) = 9007
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 3641, x1 = 5508, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
    #9008+(1870-1) = 10877
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 5509, x1 = 6535, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
    #10877+(1000-1) = 11877
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 6536, x1 = 6692 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
    #11878+(157-1) = 12034
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 6693, x1 = 7863, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
    #12035+(1088-1) = 13122
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 7864 , x1 = 12648, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
    #13123+(4730-1)= 17852
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 12649, x1 = 12980, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
    #17853+(551-1) = 18403
    
    kpPlotCoverage(kp, data=template_filt_datasets[[i]],  data.panel = 2, col = "#1414E1",  r0= -0.5, r1= -0.7, lwd = 6)
    dev.off()
    
  }
  
  
  for (i in names(nontemplate_filt_datasets)[3]){
    png(paste0("monkey_transcribed_rdna_nontemplate_KX061890_and_NR_146166_", i , "_coverage.png"), width = 10, height = 10, units= "in", res = 200)
    
    
    custom_genome <- toGRanges(data.frame(chr="KX061890_and_NR_146166", start=1, end=12980))
    kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
    
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 1, x1 = 3640 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
    #5001+(3640-1) = 9007
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 3641, x1 = 5508, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
    #9008+(1870-1) = 10877
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 5509, x1 = 6535, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
    #10877+(1000-1) = 11877
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 6536, x1 = 6692 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
    #11878+(157-1) = 12034
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 6693, x1 = 7863, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
    #12035+(1088-1) = 13122
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 7864 , x1 = 12648, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
    #13123+(4730-1)= 17852
    
    kpRect(kp, chr = 'KX061890_and_NR_146166', x0 = 12649, x1 = 12980, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
    #17853+(551-1) = 18403
    
    kpPlotCoverage(kp, data=nontemplate_filt_datasets[[i]], col = "#E21515", data.panel = 2, r0= -0.5, r1= -0.7,  lwd = 6)
    dev.off()
    
  }




