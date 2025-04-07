library(data.table)
library(tidyverse)


#From this paper 

#i found customised genome that has hg38 plus KY962518 
# However, they mentioned that they have introduced 
#"Generation of customized human and mouse genomes with rDNA.
#Reference human (NCBI Genbank KY962518.1) and mouse (NCBI Genbank BK000964.3) 
#rDNA sequences modified by cutting the original Genbank references at 35,500 and 36,000 nt positions, respectively, 
#and transposing the downstream sequence (containing the terminal part of IGS and promoters) upstream of the transcription start site".
# shared a git hub page https://github.com/vikramparalkar/rDNA-Mapping-Genomes/blob/main/README.md
#they have shared alot of custom bigwig files for polr1A, TBP, UBTF, CTCF and EP300 that can be viewed in IGV

#i wanted to add RLFS and pG4CS bed files to this session 
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/IGV/human")
rlfs<- fread("Human_KY962518.1_Modified_qmrlfs.out.bed", sep = "\t", header = FALSE) #195
pg4cs<- fread("output_pG4CS_Human_KY962518.1_Modified.txt", sep = "\t", header = FALSE) #210



# i have to add strand colur 
rlfs <- rlfs %>% mutate(V9 = ifelse(V6 == "+", "65,105,225", "205,50,120"))
rlfs$V1<- "chrR"
fwrite(rlfs, "Human_KY962518.1_Modified_qmrlfs_extended.bed", sep = "\t")
#removed header manually


pg4cs<- pg4cs %>% mutate(V4= 0) # chnage sequece information to zero.
pg4cs<- pg4cs %>% select(V1, V2, V3, V5, V4, V6) # change the order as per bed file
pg4cs<- pg4cs %>% mutate(V7 = V2)
pg4cs<- pg4cs %>% mutate(V8 =V3)
pg4cs <- pg4cs %>% mutate(V9 = ifelse(V6 == "+", "65,105,225", "205,50,120"))
pg4cs<- pg4cs %>% mutate(V10 = 2)
pg4cs<- pg4cs %>% mutate(V11 = V2)
pg4cs<-pg4cs %>% mutate(V12 = V3)
pg4cs$V1<- "chrR"


fwrite(pg4cs, "output_pG4CS_Human_KY962518.1_Modified_extended.bed", sep = "\t")
#removed header manually


##I wanted to do same for mouse 
#i wanted to add RLFS and pG4CS bed files to this session 
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/IGV/mouse")

rlfs<- fread("Mouse_BK000964.3_Modified_Snapgene_qmRLFS.out.bed", sep = "\t", header = FALSE) #68
pg4cs<- fread("Output_pG4CS_Mouse_BK000964.3_Modified_Snapgene.txt", sep = "\t", header = FALSE) #87



# i have to add strand colur 
rlfs <- rlfs %>% mutate(V9 = ifelse(V6 == "+", "65,105,225", "205,50,120"))
rlfs$V1<- "chrR"
fwrite(rlfs, "Mouse_BK000964.3_Modified_qmrlfs_extended.bed", sep = "\t")
#removed header manually


pg4cs<- pg4cs %>% mutate(V4= 0) # chnage sequece information to zero.
pg4cs<- pg4cs %>% select(V1, V2, V3, V5, V4, V6) # change the order as per bed file
pg4cs<- pg4cs %>% mutate(V7 = V2)
pg4cs<- pg4cs %>% mutate(V8 =V3)
pg4cs <- pg4cs %>% mutate(V9 = ifelse(V6 == "+", "65,105,225", "205,50,120"))
pg4cs<- pg4cs %>% mutate(V10 = 2)
pg4cs<- pg4cs %>% mutate(V11 = V2)
pg4cs<-pg4cs %>% mutate(V12 = V3)
pg4cs$V1<- "chrR"


fwrite(pg4cs, "output_pG4CS_Mouse_BK000964.3_Modified_extended.bed", sep = "\t")
#removed header manually

