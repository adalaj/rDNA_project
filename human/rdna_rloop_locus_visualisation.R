##To make figure for R-loop forming sequences in the human rDNA locus and see how it looks

library(data.table)
library(tidyverse)
library(BiocManager)
BiocManager::install("karyoploteR")
library(karyoploteR)

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rloop_and_rdna/human/one_rDNA_seq/output/Qmrlfs_results_2018")

##Two ways to do it:
#1) using karyoploteR

#read the rlfs that overlapped with rdna locus
rdna<- fread("KY962518_1_humanrDNA_2018_qmrlfs.out.bed", sep = "\t", header = FALSE)
rdna$V1= "rDNA_locus"
rdna6<- rdna %>% select(1:6)
colnames(rdna6)<- c("chr", "start", "end", "name", "score", "strand")

##separate as per strand
rdna6_nontemplate<- rdna6 %>% filter(strand=="+")
#because in NCBI keep nontemplate sequence.
rdna6_template<- rdna6 %>% filter(strand=="-")



#Plus Strand R-loops (Non-template Strand):
#R-loops annotated as being on the plus strand (non-template strand) indicate RNA is hybridized to the minus strand (template strand).
#These R-loops are formed during sense transcription, where the RNA transcript is complementary to the DNA template strand.

#Minus Strand R-loops (Template Strand):
#R-loops annotated as being on the minus strand (template strand) indicate RNA is hybridized to the plus strand (non-template strand).
#These R-loops are formed during antisense transcription, where RNA transcripts are complementary to the non-template strand.



##Initiate the steps to save the plot
##search code

#png("rdna_rlfs_fig3.png", width = 20, height= 30, units= "in", res = 150)

#Step 1: create custom genome for entire rdna locus
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=44838))

##make separate data.table for each rdna components

#red", "green", "orange", "brown", "pink", "lightyellow", "salmon", "peachpuff"

kp <- plotKaryotype(genome=custom_genome, plot.type = 2)
kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 = 3657, y0 = 0, y1 = 1, col = "red", data.panel = "ideogram", borders= NA) #marks 5'ETS
kpRect(kp, chr = 'rDNA_locus', x0 = 3658, x1 = 5526, y0 = 0, y1 = 1, col = "green", data.panel = "ideogram", borders= NA) #marks 18S
kpRect(kp, chr = 'rDNA_locus', x0 = 5527, x1 = 6596, y0 = 0, y1 = 1, col = "orange", data.panel = "ideogram", borders= NA) #marks ITS1S
kpRect(kp, chr = 'rDNA_locus', x0 = 6597, x1 = 6753, y0 = 0, y1 = 1, col = "brown", data.panel = "ideogram", borders= NA) #marks 5.8S
kpRect(kp, chr = 'rDNA_locus', x0 = 6754, x1 = 7920, y0 = 0, y1 = 1, col = "pink", data.panel = "ideogram", borders= NA) #marks ITS2
kpRect(kp, chr = 'rDNA_locus', x0 = 7921, x1 = 12971, y0 = 0, y1 = 1, col = "lightyellow", data.panel = "ideogram", borders= NA) #marks 28S
kpRect(kp, chr = 'rDNA_locus', x0 = 12972, x1 = 13332, y0 = 0, y1 = 1, col = "salmon", data.panel = "ideogram", borders= NA) #marks 3'ETS
kpRect(kp, chr = 'rDNA_locus', x0 = 13333, x1 = 44838, y0 = 0, y1 = 1, col = "peachpuff", data.panel = "ideogram", borders= NA) #marks IGS

kpPlotRegions(kp, data=rdna6_template, col="#1414E1", r0= -0.5, r1= -1.9)
kpPlotRegions(kp, data=rdna6_nontemplate, col="#E21515", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width

# I open this plots in full screen in a monitor and then took the snip using snipping tool and paste in a powerpoint.
#dev.off()




## if you want to add labels 
# kpAddLabels(kp, labels = "rDNA components", side = "left", r0 = 0.1, r1 = 0.3, cex = 1.2)
# kpAddLabels(kp, labels = "RLFS Regions", side = "left", r0 = -0.5, r1 = -0.8, cex = 1.2)


###Plot2 including promoters
 #read the rlfs that overlapped with rdna locus
 promoter_rdna<- fread("KY962518_inclu_2kb_promoter_qmrlfs.out.bed", sep = "\t", header = FALSE)
 promoter_rdna$V1= "rDNA_locus"
 promoter_rdna6<- promoter_rdna %>% select(1:6)
 colnames(promoter_rdna6)<- c("chr", "start", "end", "name", "score", "strand")
 
 ##separate as per strand
 promoter_rdna6_nontemplate<- promoter_rdna6 %>% filter(strand=="+")
 #because in NCBI keep nontemplate sequence.
 promoter_rdna6_template<- promoter_rdna6 %>% filter(strand=="-")
 
 
 ##plotting begins
 custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=46838))
 
 ##make separate data.table for each rdna components
 
 #"yellowgreen", "red", "green", "orange", "brown", "pink", "lightyellow", "salmon", "peachpuff"
 
 kp <- plotKaryotype(genome=custom_genome, plot.type = 2)
 kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =2000 , y0 = 0, y1 = 1, col = "mintcream", data.panel = "ideogram", borders= NA) #marks 2kb promoter
 kpRect(kp, chr = 'rDNA_locus', x0 = 2001, x1 = 5657 , y0 = 0, y1 = 1, col = "red", data.panel = "ideogram", borders= NA) #marks 5'ETS
 kpRect(kp, chr = 'rDNA_locus', x0 = 5658, x1 = 7526, y0 = 0, y1 = 1, col = "green", data.panel = "ideogram", borders= NA) #marks 18S
 kpRect(kp, chr = 'rDNA_locus', x0 = 7527, x1 = 8596, y0 = 0, y1 = 1, col = "orange", data.panel = "ideogram", borders= NA) #marks ITS1S
 kpRect(kp, chr = 'rDNA_locus', x0 = 8597, x1 = 8753, y0 = 0, y1 = 1, col = "brown", data.panel = "ideogram", borders= NA) #marks 5.8S
 kpRect(kp, chr = 'rDNA_locus', x0 = 8754, x1 = 9920, y0 = 0, y1 = 1, col = "pink", data.panel = "ideogram", borders= NA) #marks ITS2
 kpRect(kp, chr = 'rDNA_locus', x0 = 9921, x1 = 14971, y0 = 0, y1 = 1, col = "lightyellow", data.panel = "ideogram", borders= NA) #marks 28S
 kpRect(kp, chr = 'rDNA_locus', x0 = 14972, x1 = 15332, y0 = 0, y1 = 1, col = "salmon", data.panel = "ideogram", borders= NA) #marks 3'ETS
 kpRect(kp, chr = 'rDNA_locus', x0 = 15333, x1 = 46838, y0 = 0, y1 = 1, col = "peachpuff", data.panel = "ideogram", borders= NA) #marks IGS
 
 
 
 
 
 kpPlotRegions(kp, data=promoter_rdna6_template, col="#1414E1", r0= -0.5, r1= -1.9)
 kpPlotRegions(kp, data=promoter_rdna6_nontemplate, col="#E21515", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width

 
  
###Plotting with promoter and excluding IGS
 
 ###Plot3 including promoters
 #read the rlfs that overlapped with rdna locus
 promoter_rdna<- fread("KY962518_inclu_2kb_promoter_qmrlfs.out.bed", sep = "\t", header = FALSE) #208
 promoter_rdna$V1= "rDNA_locus"
 promoter_rdna6<- promoter_rdna %>% select(1:6)
 colnames(promoter_rdna6)<- c("chr", "start", "end", "name", "score", "strand")
 promoter_rdna6_no_igs<- promoter_rdna6 %>% filter(start <= 15332)
 
 
 ##separate as per strand
 promoter_rdna6_nontemplate_no_igs<- promoter_rdna6_no_igs %>% filter(strand=="+")
 #because in NCBI keep nontemplate sequence.
 promoter_rdna6_template_no_igs<- promoter_rdna6_no_igs %>% filter(strand=="-")
 
 
 ##plotting begins
 custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=15332))
 
 ##make separate data.table for each rdna components
 
 #"mintcream", "red", "green", "orange", "brown", "pink", "lightyellow", "salmon", "peachpuff"
 
 kp <- plotKaryotype(genome=custom_genome, plot.type = 2)
 kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =2000 , y0 = 0, y1 = 1, col = "mintcream", data.panel = "ideogram", borders= NA) #marks 2kb promoter
 kpRect(kp, chr = 'rDNA_locus', x0 = 2001, x1 = 5657 , y0 = 0, y1 = 1, col = "red", data.panel = "ideogram", borders= NA) #marks 5'ETS
 kpRect(kp, chr = 'rDNA_locus', x0 = 5658, x1 = 7526, y0 = 0, y1 = 1, col = "green", data.panel = "ideogram", borders= NA) #marks 18S
 kpRect(kp, chr = 'rDNA_locus', x0 = 7527, x1 = 8596, y0 = 0, y1 = 1, col = "orange", data.panel = "ideogram", borders= NA) #marks ITS1S
 kpRect(kp, chr = 'rDNA_locus', x0 = 8597, x1 = 8753, y0 = 0, y1 = 1, col = "brown", data.panel = "ideogram", borders= NA) #marks 5.8S
 kpRect(kp, chr = 'rDNA_locus', x0 = 8754, x1 = 9920, y0 = 0, y1 = 1, col = "pink", data.panel = "ideogram", borders= NA) #marks ITS2
 kpRect(kp, chr = 'rDNA_locus', x0 = 9921, x1 = 14971, y0 = 0, y1 = 1, col = "lightyellow", data.panel = "ideogram", borders= NA) #marks 28S
 kpRect(kp, chr = 'rDNA_locus', x0 = 14972, x1 = 15332, y0 = 0, y1 = 1, col = "salmon", data.panel = "ideogram", borders= NA) #marks 3'ETS
 
 
 
 
 kpPlotRegions(kp, data=promoter_rdna6_template_no_igs, col="#1414E1", r0= -0.5, r1= -1.9)
 kpPlotRegions(kp, data=promoter_rdna6_nontemplate_no_igs, col="#E21515", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width
 
 
 
 ##plotting with IGS on start and promoter at end.
 # i manually added 3500bp from IGS at the front. Out of which, first 1298 belongs to IGS and remaining 2202bp will be considered as promoter. 
 # rationale is that I took the sequences form hg38 and checked how many nucleotides match with IGS. I realized only 2200 match from hg38 to IGS. 
 # so the total length of KY962518_added_3500nt_IGS_upstream_nontemplate.fasta is 44838+3500 = 48338
 
 ###Plot4 showing entire rDNA and begining of next
 
 #open terminal
 #(base) jyotiadala@Jyotis-MacBook-Pro Downloads % conda activate python2.7
 #(python2.7) jyotiadala@Jyotis-MacBook-Pro Downloads % python QmRLFS-finder.py -bed -i KY962518_added_3500nt_IGS_upstream_nontemplate.fasta -o KY962518_added_3500nt_IGS_upstream_qmrlfs
 #QmRLFS-finder.py (version v1.5)
 #run on Fri Feb 21 2025 05:16:52 
 #command line: python QmRLFS-finder.py -bed -i KY962518_added_3500nt_IGS_upstream_nontemplate.fasta -o KY962518_added_3500nt_IGS_upstream_qmrlfs
 
 
 #read the rlfs that overlapped with rdna locus
 entire_rdna<- fread("KY962518_added_3500nt_IGS_upstream_qmrlfs.out.bed", sep = "\t", header = FALSE) #208, this contain double entry for promoter and IGS. 
  
 
 entire_rdna$V1= "rDNA_locus"
 entire_rdna6<- entire_rdna %>% select(1:6)
 colnames(entire_rdna6)<- c("chr", "start", "end", "name", "score", "strand")
 
 ##separate as per strand
 entire_rdna6_nontemplate<- entire_rdna6 %>% filter(strand=="+") #94
 #because in NCBI keep nontemplate sequence.
 
 
 entire_rdna6_template<- entire_rdna6 %>% filter(strand=="-")#114
 

 png("rdna_both_strand_rlfs.png", width = 15, height= 10, units= "in", res = 600) #png is supported by biorender
 
 ##plotting begins
 custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=48338))
 #end is 48338 because I added 3500 to 44838
 #changing  colors 

 kp <- plotKaryotype(genome=custom_genome, plot.type = 2)
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =1298 , y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks last 1298bp from IGS representing previous rdna 
 kpRect(kp, chr = 'rDNA_locus', x0 = 1299, x1 =3500 , y0 = 0, y1 = 1, col = "#B6FFF4", data.panel = "ideogram", borders= NA) #marks 2202 bp of  promoter
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 3501, x1 = 7157 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram", borders= NA) #marks 5'ETS (3501+(3657-1))
 #3501+(3657-1) = 7157
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 7158, x1 = 9026, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram", borders= NA) #marks 18S
 #7158+ (1869-1) = 9026
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 9027, x1 = 10096, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram", borders= NA) #marks ITS1S
 #9027+ (1070-1) = 10096 
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 10097, x1 = 10253, y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram", borders= NA) #marks 5.8S
 #10097+ (157-1) = 10253
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 10254, x1 = 11420, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram", borders= NA) #marks ITS2
 #10254+(1167-1) = 11420
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 11421, x1 = 16471, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram", borders= NA) #marks 28S
 #11421+(5051-1) = 16471
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 16472, x1 = 16832, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram", borders= NA) #marks 3'ETS
 #16472+(361-1) = 16832
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 16833, x1 = 48338, y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks IGS
 #16833+ (31506-1)= 48338
 
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 46137, x1 = 48338, y0 = 0, y1 = 1, col = "#B6FFF4", data.panel = "ideogram", borders= NA) #marks 2202bp of promoter
 #48338-2201
 
 kpPlotRegions(kp, data=entire_rdna6_template, col="#1414E1", r0= -0.5, r1= -1.9)
 kpPlotRegions(kp, data=entire_rdna6_nontemplate, col="#E21515", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width
 dev.off()

 #use zoom option, took screenshot and edited in powerpoint
 
 
 #for checking how many RLFS are formed after rule count please see junction_rloop_2018.R file
 
 
 
 

 #need to plot only from 5'ETS to 3'ETS
 png("rdna_nontemplate_rlfs_coverage.png", width = 15, height= 10, units= "in", res = 600)
 
 custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=19000))
 kp <- plotKaryotype(genome=custom_genome, plot.type = 2)
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =1298 , y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks last 1298bp from IGS representing previous rdna 
 kpRect(kp, chr = 'rDNA_locus', x0 = 1299, x1 =3500 , y0 = 0, y1 = 1, col = "#B6FFF4", data.panel = "ideogram", borders= NA) #marks 2202 bp of  promoter
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 3501, x1 = 7157 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram", borders= NA) #marks 5'ETS (3501+(3657-1))
 #3501+(3657-1) = 7157
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 7158, x1 = 9026, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram", borders= NA) #marks 18S
 #7158+ (1869-1) = 9026
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 9027, x1 = 10096, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram", borders= NA) #marks ITS1S
 #9027+ (1070-1) = 10096 
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 10097, x1 = 10253, y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram", borders= NA) #marks 5.8S
 #10097+ (157-1) = 10253
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 10254, x1 = 11420, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram", borders= NA) #marks ITS2
 #10254+(1167-1) = 11420
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 11421, x1 = 16471, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram", borders= NA) #marks 28S
 #11421+(5051-1) = 16471
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 16472, x1 = 16832, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram", borders= NA) #marks 3'ETS
 #16472+(361-1) = 16832
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 16833, x1 = 19000, y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks IGS
 
 #16472+(361-1) = 16832
 
 kpPlotCoverage(kp, data=entire_rdna6_nontemplate, col = "#E21515")
 kpPlotRegions(kp, data=entire_rdna6_nontemplate, data.panel=2, col = "#E21515")
 dev.off()
 
 
 #need to plot only from 5'ETS to 3'ETS
 png("rdna_template_rlfs_coverage.png", width = 15, height= 10, units= "in", res = 600)
 
 custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=19000))
 kp <- plotKaryotype(genome=custom_genome, plot.type = 2)
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =1298 , y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks last 1298bp from IGS representing previous rdna 
 kpRect(kp, chr = 'rDNA_locus', x0 = 1299, x1 =3500 , y0 = 0, y1 = 1, col = "#B6FFF4", data.panel = "ideogram", borders= NA) #marks 2202 bp of  promoter
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 3501, x1 = 7157 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram", borders= NA) #marks 5'ETS (3501+(3657-1))
 #3501+(3657-1) = 7157
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 7158, x1 = 9026, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram", borders= NA) #marks 18S
 #7158+ (1869-1) = 9026
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 9027, x1 = 10096, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram", borders= NA) #marks ITS1
 #9027+ (1070-1) = 10096 
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 10097, x1 = 10253, y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram", borders= NA) #marks 5.8S
 #10097+ (157-1) = 10253
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 10254, x1 = 11420, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram", borders= NA) #marks ITS2
 #10254+(1167-1) = 11420
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 11421, x1 = 16471, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram", borders= NA) #marks 28S
 #11421+(5051-1) = 16471
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 16472, x1 = 16832, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram", borders= NA) #marks 3'ETS
 #16472+(361-1) = 16832
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 16833, x1 = 19000, y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks IGS
 #16472+(361-1) = 16832
 
 kpPlotCoverage(kp, data=entire_rdna6_template, col = "#1414E1")
 kpPlotRegions(kp, data=entire_rdna6_template, data.panel=2, col = "#1414E1")
 dev.off()
 
 #visually i dot think one can differentiate where ITS2 RLFS is beginning, it seems we have alot of ITS1 RLFS.
 #to see that 
 
 #i made three seperate garphs
 its2_test <- entire_rdna6_nontemplate %>%
               filter(start > 10254 & start < 11421)
 
 its1_test<- entire_rdna6_nontemplate %>%
              filter(start > 9027 & start < 10097)
 
 its1and2_test<- entire_rdna6_nontemplate %>%
                  filter(start > 9027 & start < 11421)
 
 custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=19000))
 kp <- plotKaryotype(genome=custom_genome, plot.type = 2)
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =1298 , y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks last 1298bp from IGS representing previous rdna 
 kpRect(kp, chr = 'rDNA_locus', x0 = 1299, x1 =3500 , y0 = 0, y1 = 1, col = "#B6FFF4", data.panel = "ideogram", borders= NA) #marks 2202 bp of  promoter
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 3501, x1 = 7157 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram", borders= NA) #marks 5'ETS (3501+(3657-1))
 #3501+(3657-1) = 7157
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 7158, x1 = 9026, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram", borders= NA) #marks 18S
 #7158+ (1869-1) = 9026
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 9027, x1 = 10096, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram", borders= NA) #marks ITS1S
 #9027+ (1070-1) = 10096 
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 10097, x1 = 10253, y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram", borders= NA) #marks 5.8S
 #10097+ (157-1) = 10253
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 10254, x1 = 11420, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram", borders= NA) #marks ITS2
 #10254+(1167-1) = 11420
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 11421, x1 = 16471, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram", borders= NA) #marks 28S
 #11421+(5051-1) = 16471
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 16472, x1 = 16832, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram", borders= NA) #marks 3'ETS
 #16472+(361-1) = 16832
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 16833, x1 = 19000, y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks IGS
 
 #16472+(361-1) = 16832
 
 
 kpPlotCoverage(kp, data=its1_test, col = "#E21515")
 kpPlotRegions(kp, data=its1_test, data.panel=2, col = "#E21515")
 
 kpPlotCoverage(kp, data=its1and2_test, col = "#E21515")
 kpPlotRegions(kp, data=its1and2_test, data.panel=2, col = "#E21515")
 
 kpPlotCoverage(kp, data=its2_test, col = "#E21515")
 kpPlotRegions(kp, data=its2_test, data.panel=2, col = "#E21515")
 
 # i just took teh screenshot not saved as pdf. 
 
 #after looking at it i can see that in the final image, you can differentiate the the one that starts from ITS2 are actually present at the bottom
 #so last 7 RLFS between ITS1 to IT2 actually belong to ITS2
 #bruce told table will have cooridnates so we dont need to show separate ITS RLFS
 
 
 
 
 #I am trying to make bedgraphs directly in karyoplote
 ##plotting begins for entire region till IGS
 
 entire_rdna<- fread("KY962518_added_3500nt_IGS_upstream_qmrlfs.out.bed", sep = "\t", header = FALSE) #208, this contain double entry for promoter and IGS. 
 entire_rdna$V1= "rDNA_locus"
 entire_rdna6<- entire_rdna %>% select(1:6)
 colnames(entire_rdna6)<- c("chr", "start", "end", "name", "score", "strand")
 
 ##separate as per strand
 entire_rdna6_nontemplate<- entire_rdna6 %>% filter(strand=="+") #94
 #because in NCBI keep nontemplate sequence.
 
 
 entire_rdna6_template<- entire_rdna6 %>% filter(strand=="-")#114
 
 
 ##non template
 png("entire_rdna_nontemplate_rlfs_coverage.png", width = 30, height= 30, units= "in", res = 150)
 custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=48338))
 kp <- plotKaryotype(genome=custom_genome, plot.type = 2)
 
kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =1298 , y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks last 1298bp from IGS representing previous rdna 
 kpRect(kp, chr = 'rDNA_locus', x0 = 1299, x1 =3500 , y0 = 0, y1 = 1, col = "#B6FFF4", data.panel = "ideogram", borders= NA) #marks 2202 bp of  promoter
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 3501, x1 = 7157 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram", borders= NA) #marks 5'ETS (3501+(3657-1))
 #3501+(3657-1) = 7157
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 7158, x1 = 9026, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram", borders= NA) #marks 18S
 #7158+ (1869-1) = 9026
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 9027, x1 = 10096, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram", borders= NA) #marks ITS1S
 #9027+ (1070-1) = 10096 
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 10097, x1 = 10253, y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram", borders= NA) #marks 5.8S
 #10097+ (157-1) = 10253
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 10254, x1 = 11420, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram", borders= NA) #marks ITS2
 #10254+(1167-1) = 11420
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 11421, x1 = 16471, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram", borders= NA) #marks 28S
 #11421+(5051-1) = 16471
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 16472, x1 = 16832, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram", borders= NA) #marks 3'ETS
 #16472+(361-1) = 16832
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 16833, x1 = 48338, y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks IGS
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 46137, x1 = 48338, y0 = 0, y1 = 1, col = "#B6FFF4", data.panel = "ideogram", borders= NA) #marks 2202bp of promoter
 #48338-2201
 #16472+(361-1) = 16832
 kpPlotCoverage(kp, data=entire_rdna6_nontemplate, col = "#E21515")
 kpPlotRegions(kp, data=entire_rdna6_nontemplate, data.panel=2, col = "#E21515")
 dev.off()
 
 
 #template
 png("entire_rdna_template_rlfs_coverage.png", width = 30, height= 30, units= "in", res = 150)
 custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=48338))
 kp <- plotKaryotype(genome=custom_genome, plot.type = 2)
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =1298 , y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks last 1298bp from IGS representing previous rdna 
 kpRect(kp, chr = 'rDNA_locus', x0 = 1299, x1 =3500 , y0 = 0, y1 = 1, col = "#B6FFF4", data.panel = "ideogram", borders= NA) #marks 2202 bp of  promoter
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 3501, x1 = 7157 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram", borders= NA) #marks 5'ETS (3501+(3657-1))
 #3501+(3657-1) = 7157
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 7158, x1 = 9026, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram", borders= NA) #marks 18S
 #7158+ (1869-1) = 9026
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 9027, x1 = 10096, y0 = 0, y1 = 1, col = "#EF9B20", data.panel = "ideogram", borders= NA) #marks ITS1S
 #9027+ (1070-1) = 10096 
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 10097, x1 = 10253, y0 = 0, y1 = 1, col = "#A0322B", data.panel = "ideogram", borders= NA) #marks 5.8S
 #10097+ (157-1) = 10253
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 10254, x1 = 11420, y0 = 0, y1 = 1, col = "#FFCC17", data.panel = "ideogram", borders= NA) #marks ITS2
 #10254+(1167-1) = 11420
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 11421, x1 = 16471, y0 = 0, y1 = 1, col = "#E5FFB6", data.panel = "ideogram", borders= NA) #marks 28S
 #11421+(5051-1) = 16471
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 16472, x1 = 16832, y0 = 0, y1 = 1, col = "#3B8CC4", data.panel = "ideogram", borders= NA) #marks 3'ETS
 #16472+(361-1) = 16832
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 16833, x1 = 48338, y0 = 0, y1 = 1, col = "#A4A2A8", data.panel = "ideogram", borders= NA) #marks IGS
 
 kpRect(kp, chr = 'rDNA_locus', x0 = 46137, x1 = 48338, y0 = 0, y1 = 1, col = "#B6FFF4", data.panel = "ideogram", borders= NA) #marks 2202bp of promoter
 #48338-2201
 #16472+(361-1) = 16832
 
 kpPlotCoverage(kp, data=entire_rdna6_template, col = "#1414E1")
 kpPlotRegions(kp, data=entire_rdna6_template, data.panel=2, col = "#1414E1")
 dev.off()
 
 

 
 
 
 
 
 
 
 
 
 
 