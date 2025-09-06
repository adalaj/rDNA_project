library(Biostrings)

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/input")

hg38_rDNA<- readDNAStringSet(file = "KY962518_added_3500nt_IGS_upstream_nontemplate.fasta", format = "FASTA")
nchar(hg38_rDNA)

vmatchPattern(pattern_dna, hg38_rDNA[[1]])

#I took hairpin sequence from paper: The structure of NoRC‐associated RNA is crucial for targeting the chromatin remodelling complex NoRC to the nucleolus. 
# figure 1 and wrote teh sequence in the sublime text. They specified position but not sure which reference they used.
# to see which reference rdna they have used for position, i looked their online supplementary file embr2008109-sup-0001.txt
# it says the first line that they have used Homo sapiens rDNA promoter (X01547).
#if i google, then i coulse see this number in link:https://www.ncbi.nlm.nih.gov/nuccore/U13369 and search X01547
# it says that misc_feature    42486..42999
#note="corresponds to GenBank Accession Numbers J00304,
# X01547, K01105, M30906, M13001"


#The conservation landscape of the human ribosomal RNA gene repeats, this paper they have done comparison of genebank U13369.1 and KY962518.
#i downloaded this pdf took the S5 Appendix. Sequence alignment between human rDNA reference sequence and KY962518
#I aim to find the match for hairpin sequence in my KY962518 rdna and if there a match i will find where exactly it is. 
# if not then i will search Rloop in U13369 just in promoter.

# first i dont know which direction, if this hairpin is by Pol I or II.
# mostly it will pRNA so POL I, meaning template strand. 



#Paper says that from -137  to -50 wrt to TSS the reading teh hairpin would look like (fig1F)
hairpin_seq <- "CGAUGGUGGCGUUUUUGGGGACAGGUGUCCGUGUCGCGUGUCGCGCGUCGCCUGGGCCGGCGGCGUGGUCGGUGACGCGACCUCCCGGCCCCGGGGA"


#if you wanna read for -50 t0 -137 then simply reverse test sequnece
#this is how it would look like instead 
"AGGGGCCCCGGCCCUCCAGCGCAGUGGCUGGUGCGGCGGCCGGGUCCGCUGCGCGCUGUGCGCUGUGCCUGUGGACAGGGGUUUUUGCGGUGGUAGC"


#here we see alot of U meaning its RNA, I need DNA
rna_hairpin_seq <- DNAString(gsub("U", "T", hairpin_seq))

dna_hairpin_seq<- complement(rna_hairpin_seq)
dna_hairpin_seq
#97-letter DNAString object
#seq: GCTACCACCGCAAAAACCCCTGTCCACAGGCACAGCGCACAGCGCGCAGCGGACCCGGCCGCCGCACCAGCCACTGCGCTGGAGGGCCGGGGCCCCT  #this isi what rdna nontemplate woudl be made up of

#when i looked this sequence in KY962518_added_3500nt_IGS_upstream_nontemplate.fasta using sublime text, i didnt see a match.
# only these GCTACCAC match from left side (start) and GGCCCCT from right side (end)

#this means that may be we see this hairpin in template strand instead 

template_dna_hairpin_seq<- complement(dna_hairpin_seq)
template_dna_hairpin_seq
#97-letter DNAString object
#seq: CGATGGTGGCGTTTTTGGGGACAGGTGTCCGTGTCGCGTGTCGCGCGTCGCCTGGGCCGGCGGCGTGGTCGGTGACGCGACCTCCCGGCCCCGGGGA # this will change reading from -50 to -137 
     #CGATGGTGGCGTTTTTGGGGACAGGTGTCCGTGTCGCG  from left side (start)                and this GCCCCGGGGA from left side

# to be sure if this the right or rough sequence 
# i looked at fig1 A, they have provided reference, which looks identical to U13369.1 and some what similar to KY962518 from position +148

#so this is match sequence
    #CGATGGTGGCGTTTTTGGGGACAGGTGTCCGTGTCGCG  from left side (start) 
 #CTGCGATGGTGGCGTTTTTGGGGACAGGTGTCCGTGTCCGTGTCGCGCGTCGCCTGGGCC (this is from the reference of 1A, also matching identical to U13)
    #CGAUGGUGGCGUUUUUGGGGACAGGUGUCCGUGUCGCGUGUCGCGCGUCGCCUGGGCCGGCGGCGUGGUCGGUGACGCGACCUCCCGGCCCCGGGGA (this is hairpin, -137 to -50)
 #CTGCGATGGTGGCGTTTTTGGGGACAGGTGTCCGTGTCCGTGTCGCGCGTCGCCTGGGCCGGCGGCGTGGTCGGTGACGCGACCTCCCGGCCCCGGGGGA


#next i wanted see if U13369.1 and KY962518 from end (promoter ) are same or not.
#from this paper The conservation landscape of the human ribosomal RNA gene repeats
#i downloaded this pdf took the S5 Appendix. Sequence alignment between human rDNA reference sequence and KY962518
# they have done comparison of genebank U13369.1 and KY962518
# interesting observation is that towards teh end they say both are identical but 
# in reality they are not 
#https://www.ncbi.nlm.nih.gov/nuccore/U13369 thsi end wil gggttatat whereas 
#https://www.ncbi.nlm.nih.gov/nuccore/KY962518 this end with gggttatt

#this the last 140 nucleotides for KY
#CTGCGATGGTGGCGTTTTTGGGGACAGGTGTCCGTGTCG       CGCGTCGCCTGGGCCGGCGGCGTGGTCGGTGACGCGACCTCCCGGCCCCGGGGGAGGTATATCTTTCGCTCCGAGTCGGCA TTTTGGGCCGCCGGGTTATT

# and this the from U133.              ****
#CTGCGATGGTGGCGTTTTTGGGGACAGGTGTCCGTGTCCGT GTCGCGCGTCGCCTGGGCCGGCGGCGTGGTCGGTGACGCGACCTCCCGGCCCCGGGGGAGGTATATCTTTCGCTCCGAGTCGGCAATTTTGGGCCGCCGGGTTAT AT
#They are similar but not identical which was claimed by S5 Appendix. Sequence alignment between human rDNA reference sequence and KY962518.pdf

#so plan is to make visualization just for my satisfaction and if bruce agrees then i can add to discussion later.

#for visualization i will use U13369 dna sequence only 2KB promoter and 5'ETS region
#1..3656 is "5' external transcribed spacer"
# and 1..42999 is entire

# from paper it says -143 to -22 is tip5 binding site
# -137 to -50 predicted hairpin pRNA

#here is table for that
prna_details<- data.table(chr = c("rDNA_locus", "rDNA_locus"),
                          start = c(1857, 1863), #2000-137 = 1863, #2000-143 = 1857
                          end = c(1978, 1950), #2000-22 = 1978, #2000-50 = 1950
                          name = c("hairpin_seq", "pRNA"), 
                          score = c(0,0),
                          strand = c("+", "+"))


setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/input")
U13369_rDNA<- readDNAStringSet(file = "U13369_1_humanrDNA_2010_nontemplate.fasta", format = "fasta", use.names = FALSE)
U13369_rDNA_seq<- as.character(U13369_rDNA[[1]])
nchar(U13369_rDNA_seq)
#42999

ets5<- subseq(U13369_rDNA_seq, start= 1, end = 3656) 
promoter<- subseq(U13369_rDNA_seq, start= 40999, end = 42999)
s18<- subseq(U13369_rDNA_seq, start = 3657, end = 5527)


promoter_ets5_s18<- paste0(promoter, ets5, s18)
nchar(promoter_ets5_s18)
#7528


  write(">U13369_promoter_ets5_s18",                                            
        file = "pRNA_U13369_rdna.fasta",
        append = TRUE)

  write(promoter_ets5_s18,                                            
        file = "pRNA_U13369_rdna.fasta",
        append = TRUE)
  


  # i simply used wedsite and got the output
  
  entire_rdna<- fread("pRNA_U13369_rdna_qmrlfs.out.bed", sep = "\t", header = FALSE) #208, this contain double entry for promoter and IGS. 
  
  
  entire_rdna$V1= "rDNA_locus"
  entire_rdna6<- entire_rdna %>% select(1:6)
  colnames(entire_rdna6)<- c("chr", "start", "end", "name", "score", "strand")
  
  ##separate as per strand
  entire_rdna6_nontemplate<- entire_rdna6 %>% filter(strand=="+") #94
  #because in NCBI keep nontemplate sequence.
  
  
  entire_rdna6_template<- entire_rdna6 %>% filter(strand=="-")#114
  
  
  
  custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=7528))
  
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2)
  
  kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =2000 , y0 = 0, y1 = 1, col = "#B6FFF4", data.panel = "ideogram", borders= NA) #marks 2202 bp of  promoter
  
  kpRect(kp, chr = 'rDNA_locus', x0 = 2000, x1 = 5657 , y0 = 0, y1 = 1, col = "#FDCCE5", data.panel = "ideogram", borders= NA) #marks 5'ETS (3501+(3657-1))
  #2000+3657 = 5657
  
  kpRect(kp, chr = 'rDNA_locus', x0 = 5658, x1 = 7528, y0 = 0, y1 = 1, col = "#D0B6FF", data.panel = "ideogram", borders= NA) #marks 18S
  #5658+1870 = 9026
  
  kpPlotRegions(kp, data=test, col="black", r0= -0.5, r1= -0.6) # you will see two black lines
  kpPlotRegions(kp, data=entire_rdna6_template, col="#1414E1", r0= -0.5, r1= -1.9)
  kpPlotRegions(kp, data=entire_rdna6_nontemplate, col="#E21515", r0= -0.5, r1= -1.3) #-1.5 to make blue with more width
  
  
  


  entire_rdna6<- entire_rdna6 %>% setkey(chr, start, end)
  prna_details<-  prna_details %>% setkey(chr, start, end)
  
  prna_ovlp_rlfs <- foverlaps(prna_details,entire_rdna6) %>% na.omit()
fwrite(prna_ovlp_rlfs, "prna_haripin_ovlp_rlfs.csv")
