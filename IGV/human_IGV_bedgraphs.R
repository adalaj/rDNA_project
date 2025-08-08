# to make bedgraphs for IGV 
library(tidyverse)
library(data.table)
library(Biostrings)


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/IGV/human/IGV_files/bedfiles")


#From this paper 

#i found customized genome that has hg38 plus KY962518 
# However, they mentioned that they have introduced 
#"Generation of customized human and mouse genomes with rDNA.
#Reference human (NCBI Genbank KY962518.1) and mouse (NCBI Genbank BK000964.3) 
#rDNA sequences modified by cutting the original Genbank references at 35,500 and 36,000 nt positions, respectively, 
#and transposing the downstream sequence (containing the terminal part of IGS and promoters) upstream of the transcription start site".
# shared a git hub page https://github.com/vikramparalkar/rDNA-Mapping-Genomes/blob/main/README.md
#they have shared alot of custom bigwig files for polr1A, TBP, UBTF, CTCF and EP300 that can be viewed in IGV


#i used python code to convert snapgene.dna file to FASTA file
#from Bio import SeqIO

# Convert SnapGene (.dna) to FASTA
#import os
#os.chdir("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/IGV/human/Human_hg38-rDNA_genome_v1.0_annotation")

#input_file = "Human_KY962518.1_Modified_Snapgene.dna"
#output_file = "output.fasta"

#with open(output_file, "w") as fasta_out:
  #for record in SeqIO.parse(input_file, "snapgene"):
  #SeqIO.write(record, fasta_out, "fasta")

#print("Conversion complete!")

#read snapgene_to_fasta.ipynb how do i complete complete conversion

#please note here the letters are in lowercase which in incase of RLFS and G4s, it doesnt matter but in case of imotif, the im seeker couldnt find any i-motifs.


#for imotif:
# I added Human_KY962518.1_Modified.fasta to predict>upload DNA FILE (FASTA FORMAT) (https://im-seeker.org/predict)>end to end with default model>predict
# it says no imotif 
# so may be lowercase it the issue so i changed the lowercase to uppercase
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/IGV/human")
human_rDNA<- readDNAStringSet(file = "Human_KY962518.1_Modified.fasta") #belongs to package BIOSTRINGS
#this modified fasta is from paralkar paper
human_rDNA_seq<- human_rDNA[[1]]
nchar(human_rDNA_seq)
# 44838


template_human_rdna<- reverseComplement(human_rDNA_seq)
nchar(template_human_rdna)
#44838



# Make uppercase and keep as DNAString
nontemplate_upper <- DNAString(toupper(human_rDNA_seq))
template_upper <- DNAString(toupper(template_human_rdna))


nontemplate_upper_set<- DNAStringSet(nontemplate_upper)
template_upper_set<- DNAStringSet(template_upper)

# Name the sequence (important for FASTA headers)
names(nontemplate_upper_set) <- "Nontemplate_KY962518"
names(template_upper_set) <- "Template_KY962518"


# Write to FASTA
writeXStringSet(nontemplate_upper_set, filepath = "Human_KY962518.1_Modified_nontemplate.fasta")
writeXStringSet(template_upper_set, filepath = "Human_KY962518.1_Modified_template.fasta")

#for i motif 
# I added Human_KY962518.1_Modified_nontemplate.fasta and Human_KY962518.1_Modified_template.fasta to predict>upload DNA FILE (FASTA FORMAT) (https://im-seeker.org/predict)>end to end with default model>predict
# output saved as output_imotif_human_default_end_to_end_prediction_nontemplate_modified.csv (csv is the default saving option)
# and output_imotif_human_default_end_to_end_prediction_template_modified.csv



imotif_nontemplate<- fread("output_imotif_human_default_end_to_end_prediction_nontemplate_modified.csv", header = TRUE, sep = ",") #57
imotif_template<- fread("output_imotif_human_default_end_to_end_prediction_template_modified.csv", header = TRUE, sep = ",") #28
imotif_nontemplate$strand<- "+"
imotif_template$strand<- "-"

datasets<- list(
  nontemplate= imotif_nontemplate,
  template = imotif_template
)

for (i in names(datasets)) {
  filt <- datasets[[i]] 
  filt<- filt %>% select(chr, beg, end, strand) 
  filt$chr<- "chrR"
  filt$name<- "imotif"
  filt$score<- 0
  filt<- filt %>% select(chr, beg, end, name, score, strand) 
  fwrite(filt, paste0("imotif_human_modified_", i, ".bed"), sep = "\t")
}


#i wanted to add RLFS and pG4CS bed files to this session 
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/IGV/human")
rlfs<- fread("Human_KY962518.1_Modified_qmrlfs.out.bed", sep = "\t", header = FALSE) #195
entire_g4s_rdna<- fread("output_pG4CS_Human_KY962518.1_Modified.txt", sep = "\t", header = FALSE) #210

 
# you can add color in strand specific manner but fro IGV is it not needed. As. you can change the setting icon for color
# you can change column V9 for color of your interest. previously, this column was zero by qmrlfs output and in G4s it dont exist.
#rlfs <- rlfs %>% mutate(V9 = ifelse(V6 == "+", "65,105,225", "205,50,120"))  
#rlfs$V1<- "chrR"
#fwrite(rlfs, "Human_KY962518.1_Modified_qmrlfs_extended.bed", sep = "\t")
#removed header manually


#pg4cs<- pg4cs %>% mutate(V4= 0) # chnage sequece information to zero.
#pg4cs<- pg4cs %>% select(V1, V2, V3, V5, V4, V6) # change the order as per bed file
#pg4cs<- pg4cs %>% mutate(V7 = V2)
#pg4cs<- pg4cs %>% mutate(V8 =V3)
#pg4cs <- pg4cs %>% mutate(V9 = ifelse(V6 == "+", "65,105,225", "205,50,120"))
#pg4cs<- pg4cs %>% mutate(V10 = 2)
#pg4cs<- pg4cs %>% mutate(V11 = V2)
#pg4cs<-pg4cs %>% mutate(V12 = V3)
#pg4cs$V1<- "chrR"
#fwrite(pg4cs, "output_pG4CS_Human_KY962518.1_Modified_extended.bed", sep = "\t")
#removed header manually
#entire_g4s_rdna<- fread("output_pG4CS_Human_KY962518.1_Modified_extended.bed", header = FALSE, sep = "\t") #210
#rlfs<- fread("Human_KY962518.1_Modified_qmrlfs_extended.bed", header = FALSE, sep = "\t") #195


#prepare datasets:
datasets<- list(pG4CS=entire_g4s_rdna,
                RLFS=rlfs)

for (i in names(datasets)) {
  filt <- datasets[[i]] 
  filt<- filt %>% select(1:6) 
  nontemplate<- filt %>% filter(V6=="+")
  #because in NCBI keep nontemplate sequence.
  template<- filt %>% filter(V6=="-")
  fwrite(nontemplate, paste0(i, "_human_modified_nontemplate.bed"), sep = "\t")
  fwrite(template, paste0(i, "_human_modified_template.bed"), sep = "\t")
}

# i removed manually all col names from this file which is V1 to V6
# and moved to different folder called Bedgraphs



#removed the headers manually

#to make bed graphs
#created a genoem bed that has chrR 1 and 44838 tab separated file
#(base) jyotiadala@Mac ~ % cd /Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/IGV/human/IGV_files/Bedgraphs/input_bedfiles/
#(base) jyotiadala@Jyotis-MacBook-Pro input_bedfiles % bedtools genomecov -i pG4CS_human_modified_nontemplate.bed -g genome.bed -bga -trackline> pG4CS_human_modified_nontemplate.bedgraph
#(base) jyotiadala@Jyotis-MacBook-Pro input_bedfiles % bedtools genomecov -i pG4CS_human_modified_template.bed -g genome.bed -bga -trackline> pG4CS_human_modified_template.bedgraph 
#(base) jyotiadala@Jyotis-MacBook-Pro input_bedfiles % bedtools genomecov -i imotif_human_modified_nontemplate.bed -g genome.bed -bga -trackline> imotif_human_modified_nontemplate.bedgraph
#(base) jyotiadala@Jyotis-MacBook-Pro input_bedfiles % bedtools genomecov -i imotif_human_modified_template.bed -g genome.bed -bga -trackline> imotif_human_modified_template.bedgraph 
#(base) jyotiadala@Jyotis-MacBook-Pro input_bedfiles % bedtools genomecov -i RLFS_human_modified_nontemplate.bed -g genome.bed -bga -trackline> RLFS_human_modified_nontemplate.bedgraph
#(base) jyotiadala@Jyotis-MacBook-Pro input_bedfiles % bedtools genomecov -i RLFS_human_modified_template.bed -g genome.bed -bga -trackline> RLFS_human_modified_template.bedgraph 



#this was fine with IGV but for quantitative comparison and correlation, deepTools is the better choice.
# i downloaded deeptools read deeptools_bigwig_signal_correlation.txt for more details

#one of the feature i need is multiBigwigSummary for deeptools
#Given typically two or more bigWig files, multiBigwigSummary computes the average scores for each of the files in every genomic region

#the input file for this is bigwig which worked for all the POL1 and UBTF but what about RLFS, PG4CS and imotif
# google says bed to bedgraph and bedgraph to bigwig is recommended step
# bedgraph to bigwig is done by UCSC bedgraphtobigwig 

#(base) jyotiadala@Jyotis-MacBook-Pro output_bedgraphs % nano genome.chrom.sizes which is same as genome.bed
#genome.chrom.sizes/genome.bed is chrR 1 to 44838


#(base) jyotiadala@Jyotis-MacBook-Pro input_bedfiles % bedgraphToBigWig RLFS_human_modified_template.bedgraph genome.bed RLFS_human_modified_template.bw 
#(base) jyotiadala@Jyotis-MacBook-Pro input_bedfiles % bedgraphToBigWig RLFS_human_modified_nontemplate.bedgraph genome.bed RLFS_human_modified_nontemplate.bw
#(base) jyotiadala@Jyotis-MacBook-Pro input_bedfiles % bedgraphToBigWig pG4CS_human_modified_template.bedgraph genome.bed pG4CS_human_modified_template.bw
#(base) jyotiadala@Jyotis-MacBook-Pro input_bedfiles % bedgraphToBigWig pG4CS_human_modified_nontemplate.bedgraph genome.bed pG4CS_human_modified_nontemplate.bw
#(base) jyotiadala@Jyotis-MacBook-Pro input_bedfiles % bedgraphToBigWig imotif_human_modified_nontemplate.bedgraph genome.bed imotif_human_modified_nontemplate.bw
#(base) jyotiadala@Jyotis-MacBook-Pro input_bedfiles % bedgraphToBigWig imotif_human_modified_template.bedgraph genome.bed imotif_human_modified_template.bw



#i also made bigiwg files irrespective of strandwise 
#(base) jyotiadala@Mac input_bedfiles % bedtools genomecov -i Human_KY962518.1_Modified_qmrlfs.out.bed -g genome.bed -bga -trackline > Entire_RLFS_human_modified.bedgraph
#(base) jyotiadala@Mac input_bedfiles % bedtools genomecov -i output_pG4CS_Human_KY962518.1_Modified.txt  -g genome.bed -bga -trackline > Entire_pG4CS_human_modified.bedgraph

imotif_nontemplate<- fread("imotif_human_modified_nontemplate.bed", header = FALSE, sep = "\t") #57
imotif_template<- fread("imotif_human_modified_template.bed", header = FALSE, sep = "\t") #28
imotif_entire<- rbind(imotif_nontemplate, imotif_template)
fwrite(imotif_entire, "output_entire_imotif_human_modified_default_end_to_end_prediction.bed", sep="\t")

#(base) jyotiadala@Mac input_bedfiles % bedtools genomecov -i output_entire_imotif_human_modified_default_end_to_end_prediction.bed  -g genome.bed -bga -trackline > Entire_imotif_human_modified.bedgraph

#perform the same convert to bigwig
#(base) jyotiadala@Jyotis-MacBook-Pro input_bedfiles % bedgraphToBigWig Entire_RLFS_human_modified.bedgraph genome.bed Entire_RLFS_human_modified.bw
#(base) jyotiadala@Mac input_bedfiles % bedgraphToBigWig Entire_imotif_human_modified.bedgraph genome.bed Entire_imotif_human_modified.bw
#(base) jyotiadala@Mac input_bedfiles % bedgraphToBigWig Entire_pG4CS_human_modified.bedgraph genome.bed Entire_pG4CS_human_modified.bw




##I wanted to do same for mouse, 
## extended file needs to be removed. Task is pending.
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








