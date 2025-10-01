##To make figure for canonical G-quadruplex forming sequences in rDNA locus and see how it looks

library(data.table)
library(tidyverse)
library(BiocManager)
BiocManager::install("karyoploteR")
library(karyoploteR)




setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/g4s_and_rdna/human")


##Two ways to do it:
#1) using karyoploteR

 

#read the pG4CS that overlapped with rdna locus
rdna<- fread("output_pG4CS_KY962518_1_humanrDNA_2018_nontemplate.txt", sep = "\t", header = FALSE)
#210

rdna$V1= "rDNA_locus"
colnames(rdna)<- c("chr", "start", "end", "sequence", "name", "strand")

##separate as per strand
rdna_nontemplate<- rdna %>% filter(strand=="+") #60
#because in NCBI keep nontemplate sequence.
rdna_template<- rdna %>% filter(strand=="-") #150



#Plus Strand R-loops (Non-template Strand):
#R-loops annotated as being on the plus strand (non-template strand) indicate RNA is hybridized to the minus strand (template strand).
#These R-loops are formed during sense transcription, where the RNA transcript is complementary to the DNA template strand.

#Minus Strand R-loops (Template Strand):

#R-loops annotated as being on the minus strand (template strand) indicate RNA is hybridized to the plus strand (non-template strand).
#These R-loops are formed during antisense transcription, where RNA transcripts are complementary to the non-template strand.





##Initiate the steps to save the plot
##search code

#png("canonicalg4s_rdna.png", width = 30, height= 30, units= "in", res = 1000) # at this moment I have not saved the file because we might have to increase the resolution as paper writing begins 


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

kpPlotRegions(kp, data=rdna_template, col="maroon3", r0= -0.5, r1= -1.9)
kpPlotRegions(kp, data=rdna_nontemplate, col="royalblue", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width

#dev.off()




## if you want to add labels 
#kpAddLabels(kp, labels = "rDNA components", side = "left", r0 = 0.1, r1 = 0.3, cex = 1.2)
#kpAddLabels(kp, labels = "pG4CS Regions", side = "left", r0 = -0.5, r1 = -0.8, cex = 1.2)

# I open this plots in full screen in a monitor and then took the snip using snipping tool and paste in a powerpoint.

###Plot2 including promoters
#read the pG4CS that overlapped with rdna locus
promoter_rdna<- fread("output_pG4CS_KY962518_inclu_2kb_promoter_nontemplate.txt", sep = "\t", header = FALSE)
promoter_rdna$V1= "rDNA_locus"
colnames(promoter_rdna)<- c("chr", "start", "end", "sequence", "name", "strand")

##separate as per strand
promoter_rdna_nontemplate<- promoter_rdna %>% filter(strand=="+") #60
#because in NCBI keep nontemplate sequence.


promoter_rdna_template<- promoter_rdna %>% filter(strand=="-") #162

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



kpPlotRegions(kp, data=promoter_rdna_template, col="maroon3", r0= -0.5, r1= -1.9)
kpPlotRegions(kp, data=promoter_rdna_nontemplate, col="royalblue", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width



###Plotting with promoter and excluding IGS

###Plot3 including promoters
#read the pG4CS that overlapped with rdna locus
promoter_rdna<- fread("output_pG4CS_KY962518_inclu_2kb_promoter_nontemplate.txt", sep = "\t", header = FALSE)
promoter_rdna$V1= "rDNA_locus"
promoter_rdna<- promoter_rdna %>% select(1:6)
colnames(promoter_rdna)<- c("chr", "start", "end", "sequence", "name", "strand")
promoter_rdna_no_igs<- promoter_rdna %>% filter(start <= 15332)


##separate as per strand
promoter_rdna_nontemplate_no_igs<- promoter_rdna_no_igs %>% filter(strand=="+") #40
#because in NCBI keep nontemplate sequence.
promoter_rdna_template_no_igs<- promoter_rdna_no_igs %>% filter(strand=="-") #97


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





kpPlotRegions(kp, data=promoter_rdna_template_no_igs, col="maroon3", r0= -0.5, r1= -1.9)
kpPlotRegions(kp, data=promoter_rdna_nontemplate_no_igs, col="royalblue", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width


## ##plotting with IGS on start and promoter at end.
###Plot4 showing entire rDNA and begining of next
#read the pG4CS that overlapped with rdna locus

#open terminal
#(python2.7) jyotiadala@Jyotis-MacBook-Pro Downloads % conda activate python3.11
#(python3.11) jyotiadala@Jyotis-MacBook-Pro Downloads %  python g4_canonical_finder_3.11python.py KY962518_added_3500nt_IGS_upstream_nontemplate.fasta >output_pG4CS_KY962518_added_3500nt_IGS_upstream_humanrDNA.txt   

#read the pG4CS that overlapped with rdna locus
entire_rdna<- fread("output_pG4CS_KY962518_added_3500nt_IGS_upstream_humanrDNA.txt", sep = "\t", header = FALSE) #222, this contain double entry for promoter and IGS.

entire_rdna$V1= "rDNA_locus"
colnames(entire_rdna)<- c("chr", "start", "end", "sequence", "name", "strand")

##separate as per strand
entire_rdna_nontemplate<- entire_rdna %>% filter(strand=="+") #60
#because in NCBI keep nontemplate sequence.


entire_rdna_template<- entire_rdna %>% filter(strand=="-") #162


png("rdna_both_strand_pG4CS.png", width = 15, height= 10, units= "in", res = 1000)
##plotting begins
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=48338))
#end is 48338 because I added 3500 to 44838

##make separate data.table for each rdna components

#"yellowgreen", "red", "green", "orange", "brown", "pink", "lightyellow", "salmon", "peachpuff"

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


kpPlotRegions(kp, data=entire_rdna_template, col="#1414E1", r0= -0.5, r1= -1.9, lwd =5)
kpPlotRegions(kp, data=entire_rdna_nontemplate, col="#E21515", r0= -0.5, r1= -1.3, lwd =5) #-1.5 to make blue with more width
dev.off()
#use zoom option, took screenshot and edited in powerpoint
#for checking how many pG4CS are formed after rule count please see junction_rloop_2018.R file

#need to plot only from 5'ETS to 3'ETS
png("rdna_nontemplate_pG4CS_coverage.png", width = 10, height= 10, units= "in", res = 1000)

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

kpPlotCoverage(kp, data=entire_rdna_nontemplate, col = "#E21515", lwd =5)
kpPlotRegions(kp, data=entire_rdna_nontemplate, data.panel=2, col = "#E21515",lwd=5)
dev.off()


#need to plot only from 5'ETS to 3'ETS
png("rdna_template_pG4CS_coverage.png", width = 10, height= 10, units= "in", res = 1000)

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

kpPlotCoverage(kp, data=entire_rdna_template, col = "#1414E1", lwd =5)
kpPlotRegions(kp, data=entire_rdna_template, data.panel=2, col = "#1414E1", lwd =5)
dev.off()


#I am trying to make bedgraphs directly in karyoplote
entire_rdna<- fread("output_pG4CS_KY962518_added_3500nt_IGS_upstream_humanrDNA.txt", sep = "\t", header = FALSE) #222, this contain double entry for promoter and IGS.

entire_rdna$V1= "rDNA_locus"
colnames(entire_rdna)<- c("chr", "start", "end", "sequence", "name", "strand")

##separate as per strand
entire_rdna_nontemplate<- entire_rdna %>% filter(strand=="+") #60
#because in NCBI keep nontemplate sequence.


entire_rdna_template<- entire_rdna %>% filter(strand=="-") #162

##plotting begins
png("entire_rdna_nontemplate_pG4CS_coverage.png", width = 30, height= 30, units= "in", res = 150)
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=48338))


#wanted to plot bedgraph
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


kpPlotCoverage(kp, data=entire_rdna_nontemplate, col = "#E21515")
kpPlotRegions(kp, data=entire_rdna_nontemplate, data.panel=2, col = "#E21515")
dev.off()

#template
png("entire_rdna_template_pG4CS_coverage.png", width = 30, height= 30, units= "in", res = 150)
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


kpPlotCoverage(kp, data=entire_rdna_template, col = "#1414E1")
kpPlotRegions(kp, data=entire_rdna_template, data.panel=2, col = "#1414E1")
dev.off()

#need to chnage colors based on final comments from bruce!!! refer color.R in downloads


{
###good to know but was failed attempt
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rloop_and_rdna/human/one_rDNA_seq/output/QmpG4CS_results_2018")

promoter_rdna_pG4CS<- fread("KY962518_inclu_2kb_promoter_qmpG4CS.out.bed", sep = "\t", header = FALSE)
promoter_rdna_pG4CS$V1= "rDNA_locus"
promoter_rdna6_pG4CS<- promoter_rdna_pG4CS %>% select(1:6)
colnames(promoter_rdna6_pG4CS)<- c("chr", "start", "end", "name", "score", "strand")
promoter_rdna6_pG4CS_no_igs<- promoter_rdna6_pG4CS %>% filter(start <= 15332)


##separate as per strand
promoter_rdna6_pG4CS_nontemplate_no_igs<- promoter_rdna6_pG4CS_no_igs %>% filter(strand=="+") #84
#because in NCBI keep nontemplate sequence.
promoter_rdna6_pG4CS_template_no_igs<- promoter_rdna6_pG4CS_no_igs %>% filter(strand=="-") #95

kpPlotRegions(kp, data=promoter_rdna6_pG4CS_template_no_igs, col="maroon3", r0= -0.5, r1= -1.9)
kpPlotRegions(kp, data=promoter_rdna6_pG4CS_nontemplate_no_igs, col="royalblue", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width


#looks very bad...
}
