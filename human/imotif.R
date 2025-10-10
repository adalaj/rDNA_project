#Below code used for generation of Fig.4B to G
##plotting imotif in human rDNA locus
# go to https://im-seeker.org/

#paper: iM-Seeker: a webserver for DNA i-motifs prediction and scoring via automated machine learning

#install libraries
library(data.table)
library(tidyverse)
library(Biostrings)
library(BiocManager)
BiocManager::install("karyoploteR")
library(karyoploteR)


# I added KY962518_added_3500nt_IGS_upstream_nontemplate.fasta to predict>upload DNA FILE (FASTA FORMAT) (https://im-seeker.org/predict)>end to end with default model>predict
# I got a csv file which i saved as imotif_prediction_default_setting_nontemplate.csv

human_rDNA<- readDNAStringSet(file = "KY962518_added_3500nt_IGS_upstream_nontemplate.fasta") #belongs to package BIOSTRINGS
human_rDNA_seq<- human_rDNA[[1]]
nchar(human_rDNA_seq)
# 48338

template_human_rdna<- reverseComplement(human_rDNA_seq)
nchar(template_human_rdna)
#48338

template_human_rdna_set <- DNAStringSet(template_human_rdna)
                                        
names(template_human_rdna_set) <- "template_KY962518"
writeXStringSet(template_human_rdna_set, "KY962518_added_3500nt_IGS_upstream_template.fasta")

#go to website
# I added KY962518_added_3500nt_IGS_upstream_nontemplate.fasta to predict>upload DNA FILE (FASTA FORMAT) (https://im-seeker.org/predict)>continue>end to end with default model>predict
# I got a csv file which i renamed as human_imotif_end_to_end_prediction_default_setting_nontemplate.csv


template<- fread("human_imotif_prediction_end_to_end_prediction_default_setting_template.csv", sep = ",", header = TRUE)

template<- template %>%
  unite(imotif_sequence, c(ctract1, loop1, ctract2, loop2, ctract3, loop3, ctract4), sep = "", remove = FALSE)
template$strand<- "-"
template$chr = "KY962518"

nontemplate<- fread("human_imotif_prediction_end_to_end_prediction_default_setting_nontemplate.csv", sep = ",", header = TRUE)

nontemplate<- nontemplate %>%
  unite(imotif_sequence, c(ctract1, loop1, ctract2, loop2, ctract3, loop3, ctract4), sep = "", remove = FALSE)
nontemplate$strand <- "+"


fwrite(template, "human_imotif_prediction_end_to_end_prediction_default_setting_template.csv")
fwrite(nontemplate, "human_imotif_prediction_end_to_end_prediction_default_setting_nontemplate.csv")

master<- rbind(nontemplate, template)
master<- master %>% mutate(length= end-beg)

fwrite(master, "human_imotif_prediction_end_to_end_prediction_default_setting_master.csv")


#plotting begins
#Fig 4B

#read the imotif that overlapped with rdna locus
entire_rdna<- fread("human_imotif_prediction_end_to_end_prediction_default_setting_master.csv", sep = ",", header = TRUE) #91
entire_rdna<- entire_rdna %>% select(chr, beg, end, imotif_sequence, predict_score, predict_tranPH,strand) #bcoz this has to be in bed format
entire_rdna$chr = "rDNA_locus"


##separate as per strand
entire_rdna6_nontemplate<- entire_rdna %>% filter(strand=="+") #63
#because in NCBI keep nontemplate sequence.


entire_rdna6_template<- entire_rdna %>% filter(strand=="-") #30

##plotting begins
custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=48338))
#end is 50306 because I added 3500 to 45306

##make separate data.table for each rdna components

#"yellowgreen", "red", "green", "orange", "brown", "pink", "lightyellow", "salmon", "peachpuff"

png("human_rdna_both_strand_imotif.png", width = 15, height= 10, units= "in", res = 1000)

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

kpPlotRegions(kp, data=entire_rdna6_template, col="#1414E1", r0= -0.5, r1= -0.6,lwd=5)
kpPlotRegions(kp, data=entire_rdna6_nontemplate, col="#E21515", r0= -0.5, r1= -0.6,lwd =5 ) #-1.5 to make blue with more width
dev.off()



##wanted to plot only till 3'ets and nontemplate
##plotting begins
#Fig 4F
png("human_rdna_nontemplate_imotif_coverage.png", width = 10, height= 10, units= "in", res = 1000)
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

kpPlotRegions(kp, data=entire_rdna6_nontemplate,data.panel = 1, col = "#E21515", r0= -0.5, r1= -0.6,lwd=5)
kpPlotCoverage(kp, data=entire_rdna6_nontemplate, data.panel = 2, col = "#E21515",  r0= -0.5, r1= -0.6, lwd=5)
dev.off()


#Fig 4G
png("human_rdna_template_imotif_coverage.png", width = 10, height= 10, units= "in", res = 1000)
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

kpPlotRegions(kp, data=entire_rdna6_template,data.panel = 1, col = "#1414E1", r0= -0.5, r1= -0.6, lwd = 5)
kpPlotCoverage(kp, data=entire_rdna6_template, data.panel = 2, col = "#1414E1",  r0= -0.5, r1= -0.6, lwd =5)
dev.off()


##i want plot bar graph after the rule for imotif that has 3500 bp added to 5ETS and also has promoter at the end. 
##basically input file will be imotif_prediction_end_to_end_prediction_default_setting_master.csv

##the rule: Counting the presence of imotif where it is first detected. For example, if imotif start towards the
# end of 5'ETS but stretches to 18S then it will counted under 5'ETS. 
# if we would have defined boundary in the first place then we would have lost this junction in our counting. 

library(tidyverse)
library(data.table)


#Fig 4C to 4E
entire_imotif_rdna<- fread("imotif_prediction_end_to_end_prediction_default_setting_master.csv", sep = ",", header = TRUE) #93

# strand specificity doesnt matter here because we are cropping region of our interest
entire_imotif_rdna<- entire_imotif_rdna[entire_imotif_rdna$beg >1299 & entire_imotif_rdna$beg <46137] #total length is 48338
entire_imotif_rdna$rDNA_region <- "junction"



#strand specificity will matter here now we want to allocate imotif to rdna region sub components  based on their direction
# to do that, i will be creating a column that will have actual imotif start based on strand specificity
entire_imotif_rdna<- entire_imotif_rdna %>% mutate(actual_imotif_start = ifelse(entire_imotif_rdna$strand == "+", beg, end))
entire_imotif_rdna<- entire_imotif_rdna %>% mutate(actual_imotif_end = ifelse(entire_imotif_rdna$strand=="+", end, beg))



entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start > 1299 & entire_imotif_rdna$actual_imotif_start < 3500] <- "Promoter"
sum(entire_imotif_rdna$rDNA_region=="Promoter")
#5

entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start > 3500] <- "Promoter and 5'ETS junction"


entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start >= 3501] <- "5'ETS"

entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start > 7157 ] <- "5'ETS and 18S junction"

entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start >= 7158 ] <- "18S"
entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start > 9026 ] <- "18S and ITS1 junction"


entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start >= 9027  ] <- "ITS1"
entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start > 10096] <- "ITS1 and 5.8S junction"

entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start >= 10097 ] <- "5.8S"
entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start > 10253 ] <- "5.8S and ITS2 junction"

entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start >= 10254] <- "ITS2"
entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start > 11420 ] <- "ITS2 and 28S junction"

entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start >= 11421 ] <- "28S"
entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start > 16471 ] <- "28S and 3'ETS junction"

entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start >= 16472 ] <- "3'ETS"
entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start > 16832 ] <- "3'ETS and IGS junction"

entire_imotif_rdna$rDNA_region[entire_imotif_rdna$actual_imotif_start > 16833 & entire_imotif_rdna$actual_imotif_start < 46137] <- "IGS" # 48338-2201


entire_imotif_rdna <- entire_imotif_rdna %>%
  mutate(rDNA_region_length = case_when(
    rDNA_region == "Promoter" ~ 2202,
    rDNA_region == "5'ETS"    ~ 3657,
    rDNA_region == "18S"      ~ 1869,
    rDNA_region == "ITS1"     ~ 1070,
    rDNA_region == "5.8S"     ~ 157,
    rDNA_region =="ITS2"      ~ 1167,
    rDNA_region =="28S"       ~ 5051,
    rDNA_region == "3'ETS"    ~ 361,
    rDNA_region == "IGS"      ~ 29305, 
    TRUE ~ 0)) #unmatched make it zero


fwrite(entire_imotif_rdna, "imotif_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv", sep = ",")

entire_imotif_rdna_summary<- entire_imotif_rdna %>% group_by(rDNA_region, rDNA_region_length) %>% count()
#rows order is lost 

sum(entire_imotif_rdna_summary$n)
#[1] 85

entire_imotif_rdna_summary[,1]

#1 18S     28S     3'ETS       5'ETS       
#4 IGS    ITS1       ITS2   Promoter      



new_rows<- data.table(rDNA_region = c("Promoter and 5'ETS junction", "5'ETS and 18S junction", 
                                      "18S and ITS1 junction","ITS1 and 5.8S junction","5.8S",
                                      "5.8S and ITS2 junction","ITS2 and 28S junction",
                                      "28S and 3'ETS junction", "3'ETS and IGS junction"),
                      n= 0,
                      rDNA_region_length=0)



entire_imotif_rdna_summary<- rbind(entire_imotif_rdna_summary, new_rows)
entire_imotif_rdna_summary$rDNA_region
#[1] "18S"                         "28S"                        
#[3] "3'ETS"                       "5'ETS"                      
#[5] "IGS"                         "ITS1"                       
#[7] "ITS2"                        "Promoter"                   
#[9] "Promoter and 5'ETS junction" "5'ETS and 18S junction"     
#[11] "18S and ITS1 junction"       "ITS1 and 5.8S junction"     
#[13] "5.8S"                        "5.8S and ITS2 junction"     
#[15] "ITS2 and 28S junction"       "28S and 3'ETS junction"     
#[17] "3'ETS and IGS junction"     

row.names(entire_imotif_rdna_summary)<- entire_imotif_rdna_summary$rDNA_region
entire_imotif_rdna_summary <- entire_imotif_rdna_summary[c("Promoter","Promoter and 5'ETS junction", "5'ETS", "5'ETS and 18S junction", 
                                                     "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                     "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                     "28S and 3'ETS junction", "3'ETS", "3'ETS and IGS junction", "IGS" ),]

names(entire_imotif_rdna_summary)[2] <- "rDNA_region_length"
names(entire_imotif_rdna_summary)[3] <- "iMFS_count"
entire_imotif_rdna_summary <- entire_imotif_rdna_summary %>%
  mutate(rDNA_region_length = case_when(
    rDNA_region == "Promoter" ~ 2202,
    rDNA_region == "5'ETS"    ~ 3657,
    rDNA_region == "18S"      ~ 1869,
    rDNA_region == "ITS1"     ~ 1070,
    rDNA_region == "5.8S"     ~ 157,
    rDNA_region =="ITS2"      ~ 1167,
    rDNA_region =="28S"       ~ 5051,
    rDNA_region == "3'ETS"    ~ 361,
    rDNA_region == "IGS"      ~ 29305, 
    TRUE ~ 0)) #unmatched make it zero




entire_imotif_rdna_summary<- entire_imotif_rdna_summary %>% mutate(iMFS_density = iMFS_count/rDNA_region_length) %>% 
  mutate(iMFS_density= round(iMFS_density, 4))

fwrite(entire_imotif_rdna_summary, "imotif_KY962518_added_3500nt_IGS_upstream_at_junctn_graphinput.csv", sep = ",")



entire_imotif_rdna_summary<- fread("imotif_KY962518_added_3500nt_IGS_upstream_at_junctn_graphinput.csv", sep = ",", header = TRUE)

colnames(entire_imotif_rdna_summary)<- c("rDNA_region", "rDNA_region_length", "iMFS_count", "iMFS_density")

#Out of all iMFS across the rDNA locus, what fraction comes from each region?
entire_imotif_rdna_summary<- entire_imotif_rdna_summary %>% mutate(iMFS_proportion_perc = round((iMFS_count/sum(iMFS_count)*100),2))
#Which region contributes the largest share of iMFS overall.
#Biased toward longer regions (they naturally accumulate more iMFS simply because they have more bases)
fwrite(entire_imotif_rdna_summary, "iMFS_KY962518_added_3500nt_IGS_upstream_at_junctn_after_rule_graphinput.csv", sep = ",")


#Fig 4D
entire_imotif_rdna_summary<- fread("iMFS_KY962518_added_3500nt_IGS_upstream_at_junctn_after_rule_graphinput.csv", sep=",", header = TRUE)

imotif_rdna_summary<- entire_imotif_rdna_summary[!grepl("junction", entire_imotif_rdna_summary$rDNA_region),]


#When you use coord_flip(), the order of the factor levels in rDNA_region determines how the categories are displayed along the flipped y-axis.
#By default, the first factor level appears at the bottom when flipped, and the last factor level appears at the top.

imotif_rdna_summary$rDNA_region <- factor(imotif_rdna_summary$rDNA_region, 
                                       levels = rev(c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                      "ITS2","28S", "3'ETS", "IGS" )))



#To reverse the order so that "Promoter" appears at the top when flipped, modify the levels of the factor like this

max_value<- max(entire_imotif_rdna_summary$iMFS_density, na.rm = TRUE)

imotif_norm_3500igs_nojuntn<- ggplot(imotif_rdna_summary, aes(x= rDNA_region, y = iMFS_density, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(#title= "Normalized imotif distribution in the human rDNA locus", 
       x= "Human rDNA region", 
       y= "iMFS density", 
       fill = "rDNA")+
  scale_y_continuous(breaks= seq(0, max_value, by = 0.001), limits =c(0,max_value), expand = expansion(mult = c(0, 0.05)))+
  #geom_text(aes(label= iMFS_count, hjust= -0.2, vjust= 0.5), size= 30)+
  scale_fill_manual(values= rev(c("#B6FFF4", "#FDCCE5","#D0B6FF", "#EF9B20", "#A0322B", 
                                  "#FFCC17", "#E5FFB6", "#3B8CC4", "#A4A2A8")))+
  #guides(fill = guide_legend(reverse = TRUE))
  theme_minimal()+
  theme(axis.ticks = element_line(color = "black", linewidth = 4),
        axis.ticks.length = unit(50, "pt"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5), #face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 100),
        axis.line = element_line(color = "black", linewidth = 4),
        axis.title.x = element_text(vjust = 0.5, hjust = 0.5, color = "black"), 
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, margin = margin(r=20)),  # Center Y-axis title
        axis.text.x  = element_text(color = "black"),
        axis.text.y  = element_text(color = "black"),
        legend.position = "none")+
  coord_flip()

ggsave("Normalized_imotif_distribution_in_human_rDNA_subcomponents_after_rule.tiff", 
       plot = imotif_norm_3500igs_nojuntn, width=30,height=18, dpi=600)


#Fig 4C
max_value<- max(entire_imotif_rdna_summary$iMFS_proportion_perc, na.rm = TRUE)

iMFS_prop_3500igs_nojuntn<- ggplot(imotif_rdna_summary, aes(x= rDNA_region, y = iMFS_proportion_perc, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(#title= "Normalized iMFS distribution in the Human rDNA locus", 
    x= "Human rDNA region", 
    y= "iMFS proportion (%)", 
    fill = "rDNA")+
  scale_y_continuous(breaks= seq(0, max_value, by = 10), limits =c(0,max_value), expand = expansion(mult = c(0, 0.08)))+
  #geom_text(aes(label= iMFS_count, hjust= -0.2, vjust= 0.5), size= 30)+
  scale_fill_manual(values= rev(c("#B6FFF4", "#FDCCE5","#D0B6FF", "#EF9B20", "#A0322B", 
                                  "#FFCC17", "#E5FFB6", "#3B8CC4", "#A4A2A8")))+
  #guides(fill = guide_legend(reverse = TRUE))
  theme_minimal()+
  theme(axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.ticks.length = unit(50, "pt"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5), #face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 100),
        axis.line = element_line(color = "black", linewidth = 4),
        axis.title.x = element_text(vjust = 0.6, hjust = 0.5, colour = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, margin = margin(r=20)),   # Center Y-axis title
        axis.text.x  = element_text(color = "black"),
        axis.text.y  = element_text(color = "black"),
        legend.position = "none")+
  coord_flip()

ggsave("iMFS_proportion_distribution_in_human_rDNA_subcomponents_after_rule.tiff", 
       plot = iMFS_prop_3500igs_nojuntn, width=30,height=18, dpi=600)



#FIG 4E
#to make template and non-template
entire_imotif_rdna_summary2<- entire_imotif_rdna %>% group_by(rDNA_region, rDNA_region_length, strand) %>% count()
names(entire_imotif_rdna_summary2)[4] <- "iMFS_count"

new_rows<- data.table(rDNA_region = c("18S", "28S", "3'ETS", "5'ETS", "ITS1","Promoter", "5.8S", "5.8S"),
                      rDNA_region_length = c(0, 0,0,0,0,0, 0, 0),
                      strand= c("-", "-","-","-", "-","-", "+", "-"),
                      iMFS_count = c(0, 0,0,0,0,0, 0, 0))


entire_imotif_rdna_summary2<- rbind(entire_imotif_rdna_summary2, new_rows)

entire_imotif_rdna_summary2 <- entire_imotif_rdna_summary2 %>%
  mutate(rDNA_region_length = case_when(
    rDNA_region == "Promoter" ~ 2202,
    rDNA_region == "5'ETS"    ~ 3657,
    rDNA_region == "18S"      ~ 1869,
    rDNA_region == "ITS1"     ~ 1070,
    rDNA_region == "5.8S"     ~ 157,
    rDNA_region =="ITS2"      ~ 1167,
    rDNA_region =="28S"       ~ 5051,
    rDNA_region == "3'ETS"    ~ 361,
    rDNA_region == "IGS"      ~ 29305, 
    TRUE ~ 0)) #unmatched make it zero



entire_imotif_rdna_summary2$rDNA_region <- factor(entire_imotif_rdna_summary2$rDNA_region, 
                                               levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                          "ITS2","28S", "3'ETS", "IGS" ))

entire_imotif_rdna_summary2<- entire_imotif_rdna_summary2 %>% mutate(iMFS_density = iMFS_count/rDNA_region_length) %>% 
  mutate(iMFS_density= round(iMFS_density, 4))

colnames(entire_imotif_rdna_summary2)[4]<- "iMFS_count"
colnames(entire_imotif_rdna_summary2)[5]<- "iMFS_density"

fwrite(entire_imotif_rdna_summary2, "imotif_KY962518_added_3500nt_IGS_upstream_no_junctn_strandwise_AR_graphinput.csv")

entire_imotif_rdna_summary2<- fread("imotif_KY962518_added_3500nt_IGS_upstream_no_junctn_strandwise_AR_graphinput.csv", sep = ",", header = TRUE)


max_value<- max(entire_imotif_rdna_summary2$iMFS_density, na.rm = TRUE)

entire_imotif_rdna_summary2$rDNA_region <- factor(entire_imotif_rdna_summary2$rDNA_region, 
                                               levels = rev(c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                                              "ITS2","28S", "3'ETS", "IGS" )))


entire_imotif_rdna_summary2$strand <- factor(
  entire_imotif_rdna_summary2$strand,
  levels = rev(c("+", "-"))  # "+" = Non-template first, "-" = Template second
)

imotif_strandwise<- ggplot(entire_imotif_rdna_summary2, aes(x= rDNA_region, y = iMFS_density, fill= strand)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(#title= "Normalized imotif strandwise distribution in the human rDNA locus", 
       x= "Human rDNA region", 
       y= "iMFS density", 
       fill= NULL)+
  scale_y_continuous(breaks= seq(0, max_value, by = 0.001), limits =c(0,max_value), expand = expansion(mult = c(0, 0.08)))+
  #geom_text(aes(label= iMFS_count, hjust=-0.2, vjust=0.5), size=20, position = position_dodge(width = 0.9))+
  scale_fill_manual(values= c("+" = "#E21515", "-" = "#1414E1"), 
                    labels = c("+" = "Non-template strand", "-" = "Template strand"),
                    breaks = c("+", "-"))+
  #scale_fill_manual(values = combined_colors)+
  theme_minimal()+
  theme(axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.ticks.length = unit(50, "pt"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5), #face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 100),
        axis.line = element_line(color = "black", linewidth = 4),
        axis.title.x = element_text(vjust = 0.6, hjust = 0.5, color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, margin = margin(r=20)),   # Center Y-axis title
        axis.text.x  = element_text(color = "black"),
        axis.text.y  = element_text(color = "black"),
        legend.position = "top", 
        legend.text = element_text(size=80),
        legend.key.size = unit(3, "cm"))+
  coord_flip()

ggsave( "Normalized_strandwise_imotif_distribution_in_human_rDNA_subcomponents_AR.tiff", 
        plot = imotif_strandwise, width=30,height=18, dpi=600)



####################################
#EXTRA#

entire_imotif_rdna_summary$rDNA_region <- factor(entire_imotif_rdna_summary$rDNA_region, 
                                                 levels = c("Promoter","Promoter and 5'ETS junction", "5'ETS", "5'ETS and 18S junction", 
                                                            "18S", "18S and ITS1 junction", "ITS1", "ITS1 and 5.8S junction", "5.8S", 
                                                            "5.8S and ITS2 junction",  "ITS2", "ITS2 and 28S junction","28S", 
                                                            "28S and 3'ETS junction", "3'ETS", "3'ETS and IGS junction", "IGS" ))

max_value<- round(max(entire_imotif_rdna_summary$iMFS_density, na.rm = TRUE),4)


imotif_norm_3500igs<- ggplot(entire_imotif_rdna_summary, aes(x= rDNA_region, y = iMFS_density, fill= rDNA_region)) + 
  geom_bar(stat= 'identity', color= "black") +
  labs(title= "Normalized imotif distribution in the human rDNA locus", 
       x= "human rDNA region", 
       y= "i-motif density")+
  scale_y_continuous(breaks= seq(0, max_value, by = 0.001), limits =c(0,max_value))+
  geom_text(aes(label= iMFS_count, vjust= -0.5, size= 50))+
  scale_fill_manual(values= c( "#B6FFF4","maroon", "#FDCCE5", "steelblue", "#D0B6FF","darkviolet", "#EF9B20","burlywood2", "#A0322B", 
                               "pink", "#FFCC17","aquamarine", "#E5FFB6","greenyellow", "#3B8CC4","turquoise2", "#A4A2A8"))+
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

ggsave( "Normalized_imotif_distribution_in_human_rDNA_subcomponents_incld_junctn_AR.tiff", 
        plot = imotif_norm_3500igs, width=18,height=10, dpi=150)



nontemplate<- entire_imotif_rdna_summary2 %>% filter(strand == "+")
nontemplate$rDNA_region <- factor(nontemplate$rDNA_region, 
                                  levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                             "ITS2","28S", "3'ETS", "IGS" ))


max_value<- round(max(nontemplate$iMFS_density, na.rm = TRUE)+0.0009,4) #because previous max value was 0.0069 and our max is 0.006

imotif_nontemplate <- ggplot(nontemplate, aes(x= rDNA_region, y = iMFS_density, fill= rDNA_region)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized Non-template imotif distribution in the human rDNA locus", 
       x= "human rDNA region", 
       y= "Non-template imotif density", 
       fill= "rDNA")+
  scale_y_continuous(breaks= seq(0, max_value, by = 0.001), limits =c(0,max_value))+
  geom_text(aes(label= iMFS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  
  scale_fill_manual(values= c("#B6FFF4", "#FDCCE5", "#D0B6FF", "#EF9B20", "#A0322B", 
                              "#FFCC17", "#E5FFB6", "#3B8CC4", "#A4A2A8"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))


ggsave( "Normalized_nontemplate_imotif_distribution_in_human_rDNA_subcomponents_AR.tiff", 
        plot = imotif_nontemplate, width=18,height=10, dpi=150)


template<- entire_imotif_rdna_summary2 %>% filter(strand == "-")
template$rDNA_region <- factor(template$rDNA_region, 
                               levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", 
                                          "ITS2","28S", "3'ETS", "IGS" ))


max_value<- round(max(template$iMFS_density, na.rm = TRUE)+0.0001,4) #because max is 0.0009

imotif_template <- ggplot(template, aes(x= rDNA_region, y = iMFS_density, fill= rDNA_region)) + 
  geom_bar(stat= "identity", position ="dodge", color = "black") +
  labs(title= "Normalized Template imotif distribution in the human rDNA locus", 
       x= "human rDNA region", 
       y= "Template imotif density", 
       fill= "rDNA")+
  scale_y_continuous(breaks= seq(0, max_value, by = 0.0002), limits =c(0,max_value))+
  geom_text(aes(label= iMFS_count), vjust= -1.0, size= 6, position = position_dodge(width = 0.9))+
  
  scale_fill_manual(values= c("#B6FFF4", "#FDCCE5", "#D0B6FF", "#EF9B20", "#A0322B", 
                              "#FFCC17", "#E5FFB6", "#3B8CC4", "#A4A2A8"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 20), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Center Y-axis title
        axis.ticks.y = element_line(color = "black"))


ggsave( "Normalized_template_imotif_distribution_in_human_rDNA_subcomponents_AR.tiff", 
        plot = imotif_template, width=18,height=10, dpi=150)






