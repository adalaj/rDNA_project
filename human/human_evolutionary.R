#i wanted to show evolutionary conservation
#everything will be 5'ETS to 3'ETS
# thinking to write a code that will code and simultaneously save the graphs


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/g4s_and_rdna/human/pG4CS_at_rdna_output/files")
entire_g4s_rdna <- fread("pG4CS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", header = TRUE, sep = ",") #210
entire_g4s_rdna<- entire_g4s_rdna %>% select(GenBank_Accession, actual_pG4CS_start, actual_pG4CS_end, rDNA_region, strand)

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/output")
rlfs<- fread("RLFS_KY962518_added_3500nt_IGS_upstream_master_qmrlfs_table_after_rule.csv", header = TRUE, sep = ",") #195
rlfs<- rlfs %>% select("#Name", actual_RLFS_start, actual_RLFS_end, rDNA_region, strand)


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/imotif/output/files")
imotif<- fread("imotif_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", header = TRUE, sep = ",") #85
imotif<- imotif %>% select(chr, actual_imotif_start, actual_imotif_end, rDNA_region, strand)


#prepare datasets:

datasets<- list(pG4CS=entire_g4s_rdna,
                RLFS=rlfs, 
                imotif=imotif)


nontemplate_filt_datasets <- list()  # initialize new list
template_filt_datasets <- list()

for (i in names(datasets)) {
  filt <- datasets[[i]] %>% 
    filter(!rDNA_region %in% c("Promoter", "IGS"))
  colnames(filt)<- c("chr", "actual_start", "actual_end", "rDNA_region", "strand")
  filt$name<- i
  filt$score<- 0
  filt<- filt %>% select(chr, actual_start, actual_end, name, score, strand)
  
  filt <- filt %>% mutate(actual_start= actual_start-3500)
  filt<- filt %>% mutate(actual_end = actual_end-3500)
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
png(paste0("rdna_template_", i , "_coverage.png"), width = 30, height= 30, units= "in", res = 200)


custom_genome <- toGRanges(data.frame(chr="KY962518", start=1, end=13332))
kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")



kpRect(kp, chr = 'KY962518', x0 = 1, x1 = 3657 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS (3501+(3657-1))
#3501+(3657-1) = 7157

kpRect(kp, chr = 'KY962518', x0 = 3658, x1 = 5526, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
#7158+ (1869-1) = 9026

kpRect(kp, chr = 'KY962518', x0 = 5527, x1 = 6596, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1S
#9027+ (1070-1) = 10096 

kpRect(kp, chr = 'KY962518', x0 = 6597, x1 = 6753, y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
#10097+ (157-1) = 10253

kpRect(kp, chr = 'KY962518', x0 = 6754, x1 = 7920, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
#10254+(1167-1) = 11420

kpRect(kp, chr = 'KY962518', x0 = 7921, x1 = 12971, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
#11421+(5051-1) = 16471

kpRect(kp, chr = 'KY962518', x0 = 12972, x1 = 13332, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
#16472+(361-1) = 16832

kpPlotCoverage(kp, data=template_filt_datasets[[i]], col = "#1414E1")
dev.off()

}


for (i in names(nontemplate_filt_datasets)[1:2]){
  png(paste0("rdna_nontemplate_", i , "_coverage.png"), width = 30, height= 30, units= "in", res = 200)
  
  
  custom_genome <- toGRanges(data.frame(chr="KY962518", start=1, end=13332))
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")

  
  kpRect(kp, chr = 'KY962518', x0 = 1, x1 = 3657 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS (3501+(3657-1))
  #3501+(3657-1) = 7157
  
  kpRect(kp, chr = 'KY962518', x0 = 3658, x1 = 5526, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
  #7158+ (1869-1) = 9026
  
  kpRect(kp, chr = 'KY962518', x0 = 5527, x1 = 6596, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1S
  #9027+ (1070-1) = 10096 
  
  kpRect(kp, chr = 'KY962518', x0 = 6597, x1 = 6753, y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
  #10097+ (157-1) = 10253
  
  kpRect(kp, chr = 'KY962518', x0 = 6754, x1 = 7920, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
  #10254+(1167-1) = 11420
  
  kpRect(kp, chr = 'KY962518', x0 = 7921, x1 = 12971, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
  #11421+(5051-1) = 16471
  
  kpRect(kp, chr = 'KY962518', x0 = 12972, x1 = 13332, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
  #16472+(361-1) = 16832
  
  kpPlotCoverage(kp, data=nontemplate_filt_datasets[[i]], col = "#E21515")
  dev.off()
  
}



# for i motif i want vertical lines to be shorter
for (i in names(nontemplate_filt_datasets)[3]){
  png(paste0("rdna_nontemplate_", i , "_coverage.png"), width = 30, height= 30, units= "in", res = 200)
  
  
  custom_genome <- toGRanges(data.frame(chr="KY962518", start=1, end=13332))
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
  
  
  kpRect(kp, chr = 'KY962518', x0 = 1, x1 = 3657 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS (3501+(3657-1))
  #3501+(3657-1) = 7157
  
  kpRect(kp, chr = 'KY962518', x0 = 3658, x1 = 5526, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
  #7158+ (1869-1) = 9026
  
  kpRect(kp, chr = 'KY962518', x0 = 5527, x1 = 6596, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1S
  #9027+ (1070-1) = 10096 
  
  kpRect(kp, chr = 'KY962518', x0 = 6597, x1 = 6753, y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
  #10097+ (157-1) = 10253
  
  kpRect(kp, chr = 'KY962518', x0 = 6754, x1 = 7920, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
  #10254+(1167-1) = 11420
  
  kpRect(kp, chr = 'KY962518', x0 = 7921, x1 = 12971, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
  #11421+(5051-1) = 16471
  
  kpRect(kp, chr = 'KY962518', x0 = 12972, x1 = 13332, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
  #16472+(361-1) = 16832
  
  kpPlotCoverage(kp, data=nontemplate_filt_datasets[[i]], data.panel = 2, col = "#E21515",  r0= -0.5, r1= -0.7)
  dev.off()
  
}



# for i motif i want vertical lines to be shorter
for (i in names(template_filt_datasets)[3]){
  png(paste0("rdna_template_", i , "_coverage.png"), width = 30, height= 30, units= "in", res = 200)
  
  
  custom_genome <- toGRanges(data.frame(chr="KY962518", start=1, end=13332))
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2, chromosomes = "all")
  
  
  kpRect(kp, chr = 'KY962518', x0 = 1, x1 = 3657 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram") #marks 5'ETS (3501+(3657-1))
  #3501+(3657-1) = 7157
  
  kpRect(kp, chr = 'KY962518', x0 = 3658, x1 = 5526, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram") #marks 18S
  #7158+ (1869-1) = 9026
  
  kpRect(kp, chr = 'KY962518', x0 = 5527, x1 = 6596, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram") #marks ITS1S
  #9027+ (1070-1) = 10096 
  
  kpRect(kp, chr = 'KY962518', x0 = 6597, x1 = 6753, y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram") #marks 5.8S
  #10097+ (157-1) = 10253
  
  kpRect(kp, chr = 'KY962518', x0 = 6754, x1 = 7920, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram") #marks ITS2
  #10254+(1167-1) = 11420
  
  kpRect(kp, chr = 'KY962518', x0 = 7921, x1 = 12971, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram") #marks 28S
  #11421+(5051-1) = 16471
  
  kpRect(kp, chr = 'KY962518', x0 = 12972, x1 = 13332, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram") #marks 3'ETS
  #16472+(361-1) = 16832
  
  kpPlotCoverage(kp, data=template_filt_datasets[[i]], data.panel = 2, col = "#1414E1",  r0= -0.5, r1= -0.7)
  dev.off()
  
}


