#i wanted to show evolutionary conservation
#everything will be 5'ETS to 3'ETS
# thinking to write a code that will code and simultaneously save the graphs


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/g4s_and_rdna/chicken/files")
entire_g4s_rdna <- fread("pG4CS_KT445934_chicken_junctn_details.csv", header = TRUE, sep = ",") #65
entire_g4s_rdna<- entire_g4s_rdna %>% select(chr, actual_pG4CS_start, actual_pG4CS_end, rDNA_region, strand)

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/chicken/output/files")
rlfs<- fread("RLFS_KT445934_chicken_junctn_details.csv", header = TRUE, sep = ",") #91
rlfs<- rlfs %>% select(GenBank_Accession, actual_rlfs_start, actual_rlfs_end, rDNA_region, strand)


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/imotif/chicken")
imotif<- fread("chicken_imotif_prediction_end_to_end_prediction_default_setting_master.csv", header = TRUE, sep = ",") #85
imotif<- imotif %>% select(chr, actual_imotif_start, actual_imotif_end, rDNA_region, strand, beg, end)


#prepare datasets:

datasets<- list(pG4CS=entire_g4s_rdna,
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
  png(paste0("chicken_transcribed_rdna_template_", i , "_coverage.png"), width = 10, height = 10, units= "in", res = 600)
  
  
  custom_genome <- toGRanges(data.frame(chr="KT445934", start=1, end=11863))
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
  
  
  
  kpRect(kp, chr = 'KT445934', x0 = 1, x1 = 1836 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
  #5001+(1836-1) = 9007
  
  kpRect(kp, chr = 'KT445934', x0 = 1837, x1 = 3659, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
  #9008+(1870-1) = 10877
  
  kpRect(kp, chr = 'KT445934', x0 = 3660, x1 = 6189, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
  #10877+(1000-1) = 11877
  
  kpRect(kp, chr = 'KT445934', x0 = 6190, x1 = 6346 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
  #11878+(157-1) = 12034
  
  kpRect(kp, chr = 'KT445934', x0 = 6347, x1 = 7079, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
  #12035+(1088-1) = 13122
  
  kpRect(kp, chr = 'KT445934', x0 = 7080 , x1 = 11520, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
  #13123+(4730-1)= 17852
  
  kpRect(kp, chr = 'KT445934', x0 = 11521, x1 = 11863, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
  #17853+(551-1) = 18403
  
  kpPlotCoverage(kp, data=template_filt_datasets[[i]], col = "#1414E1", lwd = 6 )
  dev.off()
  
}


for (i in names(nontemplate_filt_datasets)){
  png(paste0("chicken_transcribed_nontemplate_", i , "_coverage.png"), width = 10, height = 10, units= "in", res = 600)
  
  
  custom_genome <- toGRanges(data.frame(chr="KT445934", start=1, end=11863))
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
  
  
  kpRect(kp, chr = 'KT445934', x0 = 1, x1 = 1836 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
  #5001+(1836-1) = 9007
  
  kpRect(kp, chr = 'KT445934', x0 = 1837, x1 = 3659, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
  #9008+(1870-1) = 10877
  
  kpRect(kp, chr = 'KT445934', x0 = 3660, x1 = 6189, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
  #10877+(1000-1) = 11877
  
  kpRect(kp, chr = 'KT445934', x0 = 6190, x1 = 6346 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
  #11878+(157-1) = 12034
  
  kpRect(kp, chr = 'KT445934', x0 = 6347, x1 = 7079, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
  #12035+(1088-1) = 13122
  
  kpRect(kp, chr = 'KT445934', x0 = 7080 , x1 = 11520, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
  #13123+(4730-1)= 17852
  
  kpRect(kp, chr = 'KT445934', x0 = 11521, x1 = 11863, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
  #17853+(551-1) = 18403
  
  kpPlotCoverage(kp, data=nontemplate_filt_datasets[[i]], col = "#E21515", lwd = 6 )
  dev.off()
  
}







#as for imotif i have not calculated promoter 
imotif$score<- 0
imotif_filt<- imotif %>% select(chr, beg, end,  rDNA_region, score,strand)
imotif_template<- imotif_filt %>% filter(strand=="-")
imotif_nontemplate<- imotif_filt %>% filter(strand=="+")



png("chicken_transcribed_template_imotif_coverage.png", width = 10, height = 10, units= "in", res = 1000)


custom_genome <- toGRanges(data.frame(chr="KT445934", start=1, end=11863))
kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")



kpRect(kp, chr = 'KT445934', x0 = 1, x1 = 1836 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
#5001+(1836-1) = 9007

kpRect(kp, chr = 'KT445934', x0 = 1837, x1 = 3659, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
#9008+(1870-1) = 10877

kpRect(kp, chr = 'KT445934', x0 = 3660, x1 = 6189, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
#10877+(1000-1) = 11877

kpRect(kp, chr = 'KT445934', x0 = 6190, x1 = 6346 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
#11878+(157-1) = 12034

kpRect(kp, chr = 'KT445934', x0 = 6347, x1 = 7079, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
#12035+(1088-1) = 13122

kpRect(kp, chr = 'KT445934', x0 = 7080 , x1 = 11520, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
#13123+(4730-1)= 17852

kpRect(kp, chr = 'KT445934', x0 = 11521, x1 = 11863, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
#17853+(551-1) = 18403

kpPlotCoverage(kp, data=imotif_template, data.panel = 2, col ="#1414E1",  r0= -0.5, r1= -0.7, lwd = 6 )
dev.off()




png("chicken_transcribed_nontemplate_imotif_coverage.png", width = 10, height = 10, units= "in", res = 1000)


custom_genome <- toGRanges(data.frame(chr="KT445934", start=1, end=11863))
kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")


kpRect(kp, chr = 'KT445934', x0 = 1, x1 = 1836 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS
#5001+(1836-1) = 9007

kpRect(kp, chr = 'KT445934', x0 = 1837, x1 = 3659, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
#9008+(1870-1) = 10877

kpRect(kp, chr = 'KT445934', x0 = 3660, x1 = 6189, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1
#10877+(1000-1) = 11877

kpRect(kp, chr = 'KT445934', x0 = 6190, x1 = 6346 , y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
#11878+(157-1) = 12034

kpRect(kp, chr = 'KT445934', x0 = 6347, x1 = 7079, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
#12035+(1088-1) = 13122

kpRect(kp, chr = 'KT445934', x0 = 7080 , x1 = 11520, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
#13123+(4730-1)= 17852

kpRect(kp, chr = 'KT445934', x0 = 11521, x1 = 11863, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
#17853+(551-1) = 18403

kpPlotCoverage(kp, data=imotif_nontemplate, data.panel = 2, col = "#E21515",  r0= -0.5, r1= -0.7, lwd = 6 )
dev.off()

