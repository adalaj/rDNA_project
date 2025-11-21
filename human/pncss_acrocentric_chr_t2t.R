# ------------------------------------------------------------------------------
# This code is part of paper: In silico Mapping of Non-Canonical DNA Structures
# Across the Human Ribosomal DNA Locus.
#
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
#   This script maps non-canonical DNA structures (G4FS, RLFS, iMFS) predicted
#   on isolated rDNA repeat units back onto their genomic positions across
#   acrocentric chromosomes (13, 14, 15, 21, 22) in the T2T-CHM13 assembly.
#   It generates template/nontemplate BED files and corresponding coverage
#   BEDGRAPH tracks for visualization and downstream quantification in the
#   UCSC Genome Browser.
#
# Major Steps:
#   1. Define acrocentric chromosomes and rDNA start coordinates (T2T-CHM13).
#   2. Load predicted G4FS calls from g4_canonical_finder output files and:
#        - Replace local rDNA coordinates with genomic coordinates.
#        - Split features into template (−) and nontemplate (+) strands.
#        - Merge across all acrocentric chromosomes.
#   3. Load predicted RLFS calls from QmRLFS-finder BED output and:
#        - Shift local coordinates by chromosome-specific rDNA_start.
#        - Split into template (−) and nontemplate (+) BED tracks.
#   4. Load predicted iMFS calls from iM-seeker CSV output and:
#        - Convert to BED-like format with genomic coordinates.
#        - Separate template and nontemplate tracks.
#   5. Write consolidated BED files for each structure and strand
#        (G4FS, RLFS, iMFS; template vs nontemplate).
#   6. Use bedtools genomecov to generate base-resolution BEDGRAPH coverage
#        tracks, suitable for UCSC Genome Browser visualization.
# Input: fasta file for chr 13,14,15,21,22 template and nontemplate
# Output: bed and bedgraphs


#install libraries
library(data.table)
library(tidyverse)


#G4FS detection for all acrocentric chromosome at rdna locus in T2T genome assembly
#open terminal
#Downloads % conda activate python3.11
# Downloads % python g4_canonical_finder_3.11python.py chr15_rdna.fasta > chr15_rdna_G4FS_t2t.txt
# Downloads % python g4_canonical_finder_3.11python.py chr14_rdna.fasta > chr14_rdna_G4FS_t2t.txt
# Downloads % python g4_canonical_finder_3.11python.py chr13_rdna.fasta > chr13_rdna_G4FS_t2t.txt
# Downloads % python g4_canonical_finder_3.11python.py chr21_rdna.fasta > chr21_rdna_G4FS_t2t.txt
# Downloads % python g4_canonical_finder_3.11python.py chr22_rdna.fasta > chr22_rdna_G4FS_t2t.txt

# Define acrocentric chromosomes and rDNA start coordinates
acrocentric_chr<- list(13,14,15,21,22)
#start i got from UCSC genome browser refer refer T2T genome assembly analysis

rdna_start <- c(
  `13` = 5770549, #please note here '13' is character
  `14` = 2099538,
  `15`= 2506443,
  `21` = 3108299,
  `22` = 4793795)

all_template_G4FS <- data.frame(
  V1 = character(),
  start = numeric(),
  end = numeric(),
  name = character(),
  score = numeric(),
  strand = character()
)


all_nontemplate_G4FS <- data.frame(
  V1 = character(),
  start = numeric(),
  end = numeric(),
  name = character(),
  score = numeric(),
  strand = character()
)



for (i in acrocentric_chr){
  filename<- paste0("chr",i, "_rdna_G4FS_t2t.txt")
  G4FS<- fread(filename, header = FALSE, sep = "\t")
  G4FS$V1<- paste0("chr",i)
  G4FS<- G4FS %>% mutate(
    start=V2+rdna_start[as.character(i)],
      end=V3+rdna_start[as.character(i)])
  G4FS$name= "G4FS"
  G4FS$score= 0
  G4FS_bed<-G4FS %>% select(V1,start,end,name,score,V6) 
  G4FS_template<- G4FS_bed %>% filter(V6=="-")
  G4FS_nontemplate<- G4FS_bed %>% filter(V6=="+")
  all_template_G4FS<- rbind(all_template_G4FS, G4FS_template)
  all_nontemplate_G4FS <- rbind(all_nontemplate_G4FS, G4FS_nontemplate)
  
}

fwrite(all_nontemplate_G4FS,"acrocentric_chr_nontemplate_rdna_G4FS_t2t.bed", sep = "\t", col.names = FALSE)
fwrite(all_template_G4FS,"acrocentric_chr_template_rdna_G4FS_t2t.bed", sep = "\t", col.names = FALSE)


# Downloads % bedtools genomecov -i acrocentric_chr_template_rdna_G4FS_t2t.bed -g t2t_acrocentric_chr_length.bed -bga -trackline > acrocentric_chr_template_rdna_G4FS_t2t.bedgraph
# Downloads % bedtools genomecov -i acrocentric_chr_nontemplate_rdna_G4FS_t2t.bed -g t2t_acrocentric_chr_length.bed -bga -trackline > acrocentric_chr_nontemplate_rdna_G4FS_t2t.bedgraph

#Both bed and bedgraph uploaded in UCSC genome browser


###################################################################################################
#RLFS detection for all acrocentric chromosome at rdna locus in T2T genome assembly
# Downloads % conda activate python2.7 
# Downloads % python QmRLFS-finder.py -bed -i chr15_rdna.fasta -o chr15_rdna_qmrlfs_t2t
#QmRLFS-finder.py (version v1.5)
#run on Mon Nov 17 2025 06:37:10 
#command line: python QmRLFS-finder.py -bed -i chr15_rdna.fasta -o chr15_rdna_qmrlfs_t2t
#Time used: 217.83 mins


# Downloads % python QmRLFS-finder.py -bed -i chr14_rdna.fasta -o chr14_rdna_qmrlfs_t2t
#QmRLFS-finder.py (version v1.5)
#run on Mon Nov 17 2025 10:52:02 
#command line: python QmRLFS-finder.py -bed -i chr14_rdna.fasta -o chr14_rdna_qmrlfs_t2t
#Time used: 11.75 mins


# Downloads % python QmRLFS-finder.py -bed -i chr13_rdna.fasta -o chr13_rdna_qmrlfs_t2t
#QmRLFS-finder.py (version v1.5)
#run on Mon Nov 17 2025 11:05:44 
#command line: python QmRLFS-finder.py -bed -i chr13_rdna.fasta -o chr13_rdna_qmrlfs_t2t
#Time used: 111.05 mins


# Downloads % python QmRLFS-finder.py -bed -i chr21_rdna.fasta -o chr21_rdna_qmrlfs_t2t
#QmRLFS-finder.py (version v1.5)
#run on Mon Nov 17 2025 13:11:22 
#command line: python QmRLFS-finder.py -bed -i chr21_rdna.fasta -o chr21_rdna_qmrlfs_t2t
#Time used: 40.76 mins

# Downloads % python QmRLFS-finder.py -bed -i chr22_rdna.fasta -o chr22_rdna_qmrlfs_t2t
#QmRLFS-finder.py (version v1.5)
#run on Mon Nov 17 2025 13:53:25 
#command line: python QmRLFS-finder.py -bed -i chr22_rdna.fasta -o chr22_rdna_qmrlfs_t2t
#Time used: 15.24 mins


acrocentric_chr<- list(13,14,15,21,22)
#start i got from UCSC genome browser refer refer T2T genome assembly analysis

rdna_start <- c(
  `13` = 5770549, #please note here '13' is character 
  `14` = 2099538,
  `15`= 2506443,
  `21` = 3108299,
  `22` = 4793795)

all_template_RLFS <- data.frame(
  V1 = character(),
  start = numeric(),
  end = numeric(),
  name = character(),
  V5 = numeric(),
  V6 = character()
)


all_nontemplate_RLFS <- all_template_RLFS



for (i in acrocentric_chr){
  filename<- paste0("chr",i, "_rdna_qmrlfs_t2t.out.bed")
  RLFS<- fread(filename, header = FALSE, sep = "\t")
  RLFS$V1<- paste0("chr",i)
  RLFS<- RLFS %>% mutate(
    start=V2+rdna_start[as.character(i)],
    end=V3+rdna_start[as.character(i)])
  RLFS$name= "RLFS"
  RLFS_bed<-RLFS %>% select(V1,start,end,name,V5,V6) 
  RLFS_template<- RLFS_bed %>% filter(V6=="-")
  RLFS_nontemplate<- RLFS_bed %>% filter(V6=="+")
  all_template_RLFS<- rbind(all_template_RLFS, RLFS_template)
  all_nontemplate_RLFS <- rbind(all_nontemplate_RLFS, RLFS_nontemplate)
  
}

fwrite(all_nontemplate_RLFS,"acrocentric_chr_nontemplate_rdna_RLFS_t2t.bed", sep = "\t", col.names = FALSE)
fwrite(all_template_RLFS,"acrocentric_chr_template_rdna_RLFS_t2t.bed", sep = "\t", col.names = FALSE)


# bedtools genomecov -i acrocentric_chr_template_rdna_RLFS_t2t.bed -g t2t_acrocentric_chr_length.bed -bga -trackline > acrocentric_chr_template_rdna_RLFS_t2t.bedgraph
# bedtools genomecov -i acrocentric_chr_nontemplate_rdna_RLFS_t2t.bed -g t2t_acrocentric_chr_length.bed -bga -trackline > acrocentric_chr_nontemplate_rdna_RLFS_t2t.bedgraph
#Both bed and bedgraph uploaded in UCSC genome browser


###################################################################################################
#RLFS detection for all acrocentric chromosome at rdna locus in T2T genome assembly

#used iM-seeker ((https://im-seeker.org/) predict>upload fasta (both template and nontemplate) >default setting 
#used UCSC genome browser to get template and nontemaplte info refer T2T genome assembly analysis text file

acrocentric_chr<- list(13,14,15,21,22)
#start i got from UCSC genome browser refer refer T2T genome assembly analysis

rdna_start <- c(
  `13` = 5770549, #please note here '13' is character
  `14` = 2099538,
  `15`= 2506443,
  `21` = 3108299,
  `22` = 4793795)

all_template_iMFS <- data.frame(
  chr = character(),
  beg = numeric(),
  end = numeric(),
  name = character(),
  score = numeric(),
  strand = character()
)


all_nontemplate_iMFS <- all_template_iMFS



for (i in acrocentric_chr){
  filename<- paste0("chr",i, "_iMFS_template.csv")
  iMFS_template<- fread(filename, header = TRUE, sep = ",")
  iMFS_template<- iMFS_template %>% select(chr, beg, end)
  iMFS_template$chr<- paste0("chr",i)
  iMFS_template<- iMFS_template %>% mutate(
    beg=beg+rdna_start[as.character(i)],
    end=end+rdna_start[as.character(i)])
  iMFS_template$name= "iMFS"
  iMFS_template$score= 0
  iMFS_template$strand= "-"
  all_template_iMFS<- rbind(all_template_iMFS, iMFS_template)
  
  
  filename2<- paste0("chr",i, "_iMFS_nontemplate.csv")
  iMFS_nontemplate<- fread(filename2, header = TRUE, sep = ",")
  iMFS_nontemplate<- iMFS_nontemplate %>% select(chr, beg, end)
  iMFS_nontemplate$chr<- paste0("chr",i)
  iMFS_nontemplate<- iMFS_nontemplate %>% mutate(
    beg=beg+rdna_start[as.character(i)],
    end=end+rdna_start[as.character(i)])
  iMFS_nontemplate$name= "iMFS"
  iMFS_nontemplate$score= 0
  iMFS_nontemplate$strand= "-"
  all_nontemplate_iMFS <- rbind(all_nontemplate_iMFS, iMFS_nontemplate)
  
}

fwrite(all_nontemplate_iMFS,"acrocentric_chr_nontemplate_rdna_iMFS_t2t.bed", sep = "\t", col.names = FALSE)
fwrite(all_template_iMFS,"acrocentric_chr_template_rdna_iMFS_t2t.bed", sep = "\t", col.names = FALSE)



# Downloads bedtools genomecov -i acrocentric_chr_nontemplate_rdna_iMFS_t2t.bed -g t2t_acrocentric_chr_length.bed -bga -trackline > acrocentric_chr_nontemplate_rdna_iMFS_t2t.bedgraph
# Downloads bedtools genomecov -i acrocentric_chr_template_rdna_iMFS_t2t.bed -g t2t_acrocentric_chr_length.bed -bga -trackline > acrocentric_chr_template_rdna_iMFS_t2t.bedgraph









