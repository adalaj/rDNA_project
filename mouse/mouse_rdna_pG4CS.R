## ##plotting with IGS on start and promoter at end.
###Plot4 showing entire rDNA and begining of next
#read the pG4CS that overlapped with rdna locus

#open terminal
#(python2.7) jyotiadala@Jyotis-MacBook-Pro Downloads % conda activate python3.11
#(python3.11) jyotiadala@Jyotis-MacBook-Pro Downloads %  python g4_canonical_finder_3.11python.py BK000964_added_3500nt_IGS_upstream_nontemplate.fasta >output_pG4CS_BK000964_added_3500nt_IGS_upstream_humanrDNA.txt

library(data.table)
library(tidyverse)
library(BiocManager)
BiocManager::install("karyoploteR")
library(karyoploteR)


#read the pG4CS that overlapped with rdna locus
entire_rdna<- fread("output_pG4CS_BK000964_added_3500nt_IGS_upstream_humanrDNA.txt", sep = "\t", header = FALSE) #91
entire_rdna$V1= "rDNA_locus"
colnames(entire_rdna)<- c("chr", "start", "end", "sequence", "name", "strand")

##separate as per strand
entire_rdna_nontemplate<- entire_rdna %>% filter(strand=="+") #51
#because in NCBI keep nontemplate sequence.


entire_rdna_template<- entire_rdna %>% filter(strand=="-") #40

##plotting begins
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=48806))
#end is 48806 because I added 3500 to 45306

##make separate data.table for each rdna components

#"yellowgreen", "red", "green", "orange", "brown", "pink", "lightyellow", "salmon", "peachpuff"

kp <- plotKaryotype(genome=custom_genome, plot.type = 2)

kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =1299 , y0 = 0, y1 = 1, col = "#DE9A22", data.panel = "ideogram", borders= NA) #marks last 1298bp from IGS representing previous rdna 
kpRect(kp, chr = 'rDNA_locus', x0 = 1300, x1 =3500 , y0 = 0, y1 = 1, col = "#F5FEFB", data.panel = "ideogram", borders= NA) #marks 2200 bp of  promoter

kpRect(kp, chr = 'rDNA_locus', x0 = 3501, x1 = 7507 , y0 = 0, y1 = 1, col = "#E21515", data.panel = "ideogram", borders= NA) #marks 5'ETS (3501+(3657-1))
#3501+(4007-1) = 7507

kpRect(kp, chr = 'rDNA_locus', x0 = 7508, x1 = 9377, y0 = 0, y1 = 1, col = "#5AAA46", data.panel = "ideogram", borders= NA) #marks 18S
#7508+(1870-1) = 9377

kpRect(kp, chr = 'rDNA_locus', x0 = 9378, x1 = 10377, y0 = 0, y1 = 1, col = "#F36017", data.panel = "ideogram", borders= NA) #marks ITS1
#9378+(1000-1) = 10377

kpRect(kp, chr = 'rDNA_locus', x0 = 10378, x1 = 10534, y0 = 0, y1 = 1, col = "#6B1519", data.panel = "ideogram", borders= NA) #marks 5.8S
#10378+(157-1) = 10534

kpRect(kp, chr = 'rDNA_locus', x0 = 10535, x1 = 11622, y0 = 0, y1 = 1, col = "#818689", data.panel = "ideogram", borders= NA) #marks ITS2
#10535+(1088-1) = 11622

kpRect(kp, chr = 'rDNA_locus', x0 = 11623, x1 = 16352, y0 = 0, y1 = 1, col = "#ECE612", data.panel = "ideogram", borders= NA) #marks 28S
#11623+(4730-1)= 16352

kpRect(kp, chr = 'rDNA_locus', x0 = 16353, x1 = 16903, y0 = 0, y1 = 1, col = "#E07F80", data.panel = "ideogram", borders= NA) #marks 3'ETS
#16353+(551-1) = 16903

kpRect(kp, chr = 'rDNA_locus', x0 = 16904, x1 = 48806, y0 = 0, y1 = 1, col = "#DE9A22", data.panel = "ideogram", borders= NA) #marks IGS
#16904+(31903-1)= 48338


kpRect(kp, chr = 'rDNA_locus', x0 = 46606, x1 = 48806, y0 = 0, y1 = 1, col = "#F5FEFB", data.panel = "ideogram", borders= NA) #marks 2202bp of promoter
#48806-2200


kpPlotRegions(kp, data=entire_rdna_template, col="maroon3", r0= -0.5, r1= -1.3)
kpPlotRegions(kp, data=entire_rdna_nontemplate, col="royalblue", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width

#use zoom option, took screenshot and edited in powerpoint

##i want plot bar graph after the rule for G4s that has 3500 bp added to 5ETS and also has promoter at the end. 
##basically input file will be output_pG4CS_BK000964_added_3500nt_IGS_upstream_humanrDNA.txt

##teh rule: Counting the presence of G4s where it is first detected. For example, if G4s start towards the
# end of 5'ETS but stretches to 18S then it will counted under 5'ETS. 
# if we would have defined boundary in the first place then we would have lost this junction in our counting. 

library(tidyverse)
library(data.table)


entire_g4s_rdna<- fread("output_pG4CS_BK000964_added_3500nt_IGS_upstream_humanrDNA.txt", sep = "\t", header = FALSE) #91
entire_g4s_rdna<- entire_g4s_rdna %>% filter(V2>1299)#87
entire_g4s_rdna$v7 <- "junction"
entire_g4s_rdna$v8<- entire_g4s_rdna$V3-entire_g4s_rdna$V2


entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 1300 & entire_g4s_rdna$V2 < 3500] <- "Promoter"
sum(entire_g4s_rdna$v7=="Promoter")
#0

entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 3500] <- "Promoter and 5'ETS junction"


entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 3501] <- "5'ETS"

entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 7507 ] <- "5'ETS and 18S junction"

entire_g4s_rdna$v7[entire_g4s_rdna$V2 >= 7508] <- "18S"
entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 9377 ] <- "18S and ITS1 junction"


entire_g4s_rdna$v7[entire_g4s_rdna$V2 >= 9378   ] <- "ITS1"
entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 10377] <- "ITS1 and 5.8S junction"

entire_g4s_rdna$v7[entire_g4s_rdna$V2 >= 10378 ] <- "5.8S"
entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 10534 ] <- "5.8S and ITS2 junction"

entire_g4s_rdna$v7[entire_g4s_rdna$V2 >= 10535] <- "ITS2"
entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 11622] <- "ITS2 and 28S junction"

entire_g4s_rdna$v7[entire_g4s_rdna$V2 >= 11623 ] <- "28S"
entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 16352 ] <- "28S and 3'ETS junction"

entire_g4s_rdna$v7[entire_g4s_rdna$V2 >= 16353] <- "3'ETS"
entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 16903 ] <- "3'ETS and IGS junction"

entire_g4s_rdna$v7[entire_g4s_rdna$V2 > 16904 & entire_g4s_rdna$V2 < 48806] <- "IGS"


colnames(entire_g4s_rdna) <- c("GenBank_Accession", "g4s_start", "g4s_end", "g4s_sequence", "g4s_name", "strand", "rDNA_region", "g4s_length")

fwrite(entire_g4s_rdna, "pG4CS_BK000964_added_3500nt_IGS_upstream_at_junctn_details.csv", sep = ",")

entire_g4s_rdna_summary<- entire_g4s_rdna %>% group_by(rDNA_region) %>% count()
#rows order is lost 

sum(entire_g4s_rdna_summary$n)
#[1] 87

entire_g4s_rdna_summary[,1]

#1 28S         3'ETS       5'ETS       
#4 IGS          ITS1       ITS2       



new_rows<- data.table(rDNA_region = c("Promoter", "Promoter and 5'ETS junction", "5'ETS and 18S junction", "18S", 
                                      "18S and ITS1 junction","ITS1 and 5.8S junction","5.8S",
                                      "5.8S and ITS2 junction","ITS2 and 28S junction",
                                      "28S and 3'ETS junction", "3'ETS and IGS junction"))

new_rows$n<- 0



entire_g4s_rdna_summary<- rbind(entire_g4s_rdna_summary, new_rows)
entire_g4s_rdna_summary$rDNA_region
#[1] "28S"                         "3'ETS"                      
#[3] "5'ETS"                       "IGS"                        
#[5] "ITS1"                        "ITS2"                       
#[7] "Promoter"                    "Promoter and 5'ETS junction"
#[9] "5'ETS and 18S junction"      "18S"                        
#[11] "18S and ITS1 junction"       "ITS1 and 5.8S junction"     
#[13] "5.8S"                        "5.8S and ITS2 junction"     
#[15] "ITS2 and 28S junction"       "28S and 3'ETS junction"     
#[17] "3'ETS and IGS junction"

entire_g4s_rdna_summary <- entire_g4s_rdna_summary[c(7,8,3,9,10,11,5,12,13,14,6,15,1,16,2,17,4),]

names(entire_g4s_rdna_summary)[2] <- "pG4CS_count"


entire_g4s_rdna_summary<- entire_g4s_rdna_summary %>% mutate(norm_pG4CS_count = pG4CS_count/sum(entire_g4s_rdna_summary$pG4CS_count)) %>% 
  mutate(norm_pG4CS_count= round(norm_pG4CS_count, 2))

fwrite(entire_g4s_rdna_summary, "pG4CS_BK000964_added_3500nt_IGS_upstream_at_junctn_graphinput.csv", sep = ",")

entire_g4s_rdna_summary$rDNA_region <- factor(entire_g4s_rdna_summary$rDNA_region, 
                                                levels = c("Promoter","Promoter and 5'ETS junction", "5'ETS", "5'ETS and 18S junction", 
                                                           "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                           "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                           "28S and 3'ETS junction", "3'ETS", "3'ETS and IGS junction", "IGS" ))


g4s_norm_3500igs<- ggplot(entire_g4s_rdna_summary, aes(x= rDNA_region, y = norm_pG4CS_count, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized pG4CS distribution in the Mouse rDNA locus", 
       x= "Mouse rDNA region", 
       y= "Normalized pG4CS count")+
  scale_y_continuous(breaks= seq(0, 0.60, by = 0.1), limits =c(0,0.60))+
  geom_text(aes(label= pG4CS_count, vjust= -0.5, size= 50))+
  scale_fill_manual(values= c( "#F5FEFB","maroon", "#E21515", "steelblue", "#5AAA46","darkviolet", "#F36017","burlywood2", "#6B1519", 
                               "pink", "#818689","aquamarine", "#ECE612","greenyellow", "#E07F80","turquoise2", "#DE9A22"))+
  #scale_fill_manual(values = combined_colors)+
  theme_minimal()+
  theme(axis.text.x = element_blank())+ 
  theme(panel.grid = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        theme(panel.grid = element_blank()),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +  # Center Y-axis title
  theme(axis.ticks.y = element_line(color = "black"))
#coord_flip()

ggsave( "Normalized_pG4CS_distribution_in_mouse_rDNA_subcomponents_after_rule.tiff", 
        plot = g4s_norm_3500igs, width=18,height=10, dpi=150)



