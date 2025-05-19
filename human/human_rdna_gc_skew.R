

#Task is to visually show GC skew for human RDNA locus. 

#first make GC for each region and then calculate GC skew using sliding window in entire rDNA length .


# for first analysis, I already defined a function to calculate GC skew  

#load GC skew function
source("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/rDNA_project/human/gc_skew_function.R")

library(stringr) #needed for GC skew function
library(tidyverse)
library(data.table)

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/output")
rdna_human<- fread ("rdna_hg38_chr21_2018_dataset_details_v2.csv", sep = ",", header = TRUE)

attempt1<- data.frame((matrix(nrow = 0, ncol=6)))

for ( i in 1:nrow(rdna_human)){
  gc_skew_data<- gc_skew(rdna_human$Sequences[i])
  attempt1<- rbind(attempt1, gc_skew_data)
}

selected_gc_skew_data <- attempt1 %>% select(GC_skew_value) 
rdna_human<- cbind(rdna_human, selected_gc_skew_data)


#calculated GC skew for template strand as well
attempt2<- data.frame((matrix(nrow = 0, ncol=6)))
for ( i in 1:nrow(rdna_human)){
  rev_seq<- as.character(reverseComplement(DNAString(rdna_human$Sequences[i])))
  temp_gc_skew<- gc_skew(rev_seq)
  attempt2<- rbind(attempt2, temp_gc_skew)
  
}

template_GC_skew<- attempt2 %>% select(window_seq, G_count, C_count, GC_skew_value)
colnames(template_GC_skew)<- c("template_seq", "temp_G_count", "temp_C_count", "temp_GC_skew_value")
rdna_human<- cbind(rdna_human, template_GC_skew)


# wanted to add cumulative sum of total nucleotides for future 
rdna_human<- rdna_human %>% mutate(norm_GC_skew = GC_skew_value/sum(GC_skew_value))
rdna_human<- rdna_human %>% mutate(temp_norm_GC_skew = temp_GC_skew_value/sum(GC_skew_value))

rdna_human$x_axis <- cumsum(rdna_human$Total_nucleotides)

#calculated GC perc
rdna_human<- rdna_human %>% mutate(GC_perc = round(((G+C)/(Total_nucleotides))*100,2))


#also wanted to calculate CpG or CG dinucleotide

rdna_human <- rdna_human %>%
    mutate(
      CpG_count = str_count(Sequences, "CG"),
      CpG_perc = (CpG_count/Total_nucleotides)*100,
      CpG_OE = (CpG_count * Total_nucleotides) / (G * C))



fwrite(rdna_human, "rdna_hg38_chr21_2018_dataset_details_v3.csv")




full_length_gc_skew<- ggplot(rdna_human, aes(x = x_axis, y = norm_GC_skew)) +
  geom_line(color = "#E21515") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "GC Skew Across rDNA Sequence KY962518",
       x = "rDNA",
       y = "Normalized GC Skew") +
  scale_x_continuous(breaks = c(2202, 5859, 7728,8798, 8955, 10122, 15173, 15534, 47040 ),
                     labels = c("promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS", "IGS"))+
                  
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=20))

ggsave("KY962518_full_length_normalised_gc_skew.tiff", 
       plot = full_length_gc_skew, width=18, height = 10, dpi = 150)


full_length_gc_skew_excld_igs<- ggplot(rdna_human[1:8,], aes(x = x_axis, y = norm_GC_skew)) +
  geom_line(color = "#E21515") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "GC Skew Across rDNA Sequence KY962518 (Exculding IGS)",
       x = "rDNA",
       y = "Normalised GC Skew") +
  scale_x_continuous(breaks = c(2202, 5859, 7728,8798, 8955, 10122, 15173, 15534 ),
                     labels = c("promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS"))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),#if you want to add a rectangle box or you can use theme_minimal()
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=20))

ggsave("KY962518_full_length_excld_IGS_normalised_gc_skew.tiff", 
       plot = full_length_gc_skew_excld_igs, width=18, height = 10, dpi = 150)

#plotting GC content
full_length_gc_content_excld_igs<- ggplot(rdna_human[1:8,], aes(x = x_axis, y = GC_perc)) +
  geom_line(color = "#E21515")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 40.89, linetype = "dashed", color="darkgreen")+
  labs(title = "GC Content Across rDNA Sequence KY962518 (Exculding IGS)",
       x = "rDNA",
       y = "GC Content") +
  scale_x_continuous(breaks = c(2202, 5859, 7728,8798, 8955, 10122, 15173, 15534 ),
                     labels = c("promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS"))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),#if you want to add a rectangle box or you can use theme_minimal()
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=20))

ggsave("KY962518_full_length_excld_IGS_GC_content.tiff", 
       plot = full_length_gc_content_excld_igs, width=18, height = 10, dpi = 150)



#plotting CpG dinucleotide content
full_length_CpG_excld_igs<- ggplot(rdna_human[1:8,], aes(x = x_axis, y = CpG_OE)) +
  geom_line(color = "#E21515")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0.6, linetype = "dashed", color="darkorange")+
  labs(title = "CpG dinucleotides O/E ratio Across rDNA Sequence KY962518 (Exculding IGS)",
       x = "rDNA",
       y = "CpG dinucleotides O/E ratio") +
  scale_x_continuous(breaks = c(2202, 5859, 7728,8798, 8955, 10122, 15173, 15534 ),
                     labels = c("promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS"))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),#if you want to add a rectangle box or you can use theme_minimal()
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=20))

ggsave("KY962518_full_length_excld_IGS_CpG_dinucleotides_OE_ratio.tiff", 
       plot = full_length_CpG_excld_igs, width=18, height = 10, dpi = 150)


full_length_CpG_perc_excld_igs<- ggplot(rdna_human[1:8,], aes(x = x_axis, y = CpG_perc)) +
  geom_line(color = "#E21515")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "CpG dinucleotides rDNA Sequence KY962518 (Exculding IGS)",
       x = "rDNA",
       y = "CpG dinucleotides percent") +
  scale_x_continuous(breaks = c(2202, 5859, 7728,8798, 8955, 10122, 15173, 15534 ),
                     labels = c("promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS"))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),#if you want to add a rectangle box or you can use theme_minimal()
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=20))

ggsave("KY962518_full_length_excld_IGS_CpG_dinucleotides_perc.tiff", 
       plot = full_length_CpG_perc_excld_igs, width=18, height = 10, dpi = 150)



# there is no need to make template ones because it will be mirror image of non -template

# Above result showed that there is a GC skew rDNA region. the maximum is seen in promoter region, followed by 28s region than 5'ETS. 
# lowest is seen in 18S.

# I made a loop that will calculate GC skew and GC content in 40, 70, 100 overlapping and non-overlapping window size for rdna locus (as a whole not compartmentalized) 
#starting from 5'ETS and ending at 3'ETS).


library(Biostrings) #needed to read fasta file
#load GC skew and GC content function
#calculate GC skew
source("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/rDNA_project/human/gc_skew_function.R")

#calculate GC content
source("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/rDNA_project/human/GC_content_function.R")

#calculate CpG content
source("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/rDNA_project/human/CpG_dinucleotides_function.R")

rdna_human_seq<- readDNAStringSet("KY962518_added_3500nt_IGS_upstream_nontemplate.fasta")
seq<- as.character(rdna_human_seq[[1]]) 
seq2<- str_sub(seq, 1299, 16832) #just decided to have 5'ETS to 3'ETS
nchar(seq2) 
#[1] 15534

#bin_size=c(30, 70, 100)
bin_size = c(100)
for ( i in bin_size){
  
  gc_skew_data<- gc_skew(seq2, window_size = i)
  sliding_skew<- gc_skew_data$sliding_window_results
  fwrite(sliding_skew, paste0("KY962518_5ETS_TO_3ETS_gc_skew_sliding_data_", i,"bp.csv"))
  
  
  #if you want to run GC skew separately: 
  #for ( i in bin_size){
  #filename1<- paste0("KY962518_5ETS_TO_3ETS_gc_skew_sliding_data_", i,"bp.csv")
  #sliding_skew<- fread(filename1, sep = ",", header = TRUE)
  
  #filename2<- paste0("KY962518_5ETS_TO_3ETS_gc_skew_non_sliding_data_", i,"bp.csv")
  #non_sliding_skew<- fread(filename2, sep = ",", header = TRUE)
  
  gc_content_data<- gc_content(seq2, window_size = i)
  sliding_content<- gc_content_data$sliding_window_results
  fwrite(sliding_content, paste0("KY962518_5ETS_TO_3ETS_gc_content_sliding_data_", i,"bp.csv"))
  
  
  
  
  sliding_graph<- ggplot(sliding_skew, aes(x = start, y = GC_skew_value)) +
    geom_line(color = "#E21515") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "GC Skew Across rDNA Sequence KY962518",
         subtitle= paste0(i, " bp Sliding window size"),
         x = "Position (bp)",
        y = "GC Skew (G−C / G+C)") +
    scale_x_continuous(breaks = c(0, 2202, 5859, 7728,8798, 8955, 10122, 15173, 15534), 
                       labels =c("Promoter","5'ETS", "18S","ITS1", "5.8S", "ITS2", "28S", "3'ETS", "IGS"))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20))
  
  ggsave(paste0("KY962518_5ETS_TO_3ETS_gc_skew_sliding_", i, "bp.tiff"), 
         plot = sliding_graph, width=18, height = 10, dpi = 150)
  
  
  #sliding<- paste0("KY962518_5ETS_TO_3ETS_gc_content_sliding_data_", i,"bp.csv")
  #sliding_content<- fread(sliding, sep = ",", header = TRUE)
  
  
  sliding_content_graph<- ggplot(sliding_content, aes(x = start, y = gc_content_perc)) +
    geom_line(color = "#E21515") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 40.89, linetype = "dashed", color="darkgreen")+
    labs(title = "GC base pairs Content Across rDNA Sequence KY962518",
         subtitle= paste0(i, " bp sliding window size"),
         x = "Position (bp)",
         y = "GC Fraction (G+C / A+T+G+C)") +
    scale_y_continuous(limits= c(0,100), breaks=seq(0,100, by=20))+
    scale_x_continuous(breaks = c(0, 2202, 5859, 7728,8798, 8955, 10122, 15173, 15534), 
                       labels =c("Promoter","5'ETS", "18S","ITS1", "5.8S", "ITS2", "28S", "3'ETS", "IGS"))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20))
  
  ggsave(paste0("KY962518_5ETS_TO_3ETS_gc_content_sliding_", i, "bp.tiff"), 
         plot =  sliding_content_graph, width=18, height = 10, dpi = 150)
  
} 
  #CpG dinucleotides 
# as per Gardiner-Garden and Frommer: paper title : CpG Islands in vertebrate genomes (1987) teh CpG has to have a length of 200, GC >50%, and CpG O/E >0.6
# but in 2002 a new paper published takai and Jones (Comprehensive analysis of CpG islands in human chromosomes 21 and 22)
# they said that earlier formula was made when sequencing was not performed, contain repetitive and alu elements. 
# This demanded new formula. as per which length of 500bp, GC% >55%, and CpG O/E >0.65 is needed

bin_size = c(200, 500)
for ( i in bin_size){

  CpG_dinucleotides_data<- CpG_content(seq2, window_size = i)
  sliding_CpG<- CpG_dinucleotides_data$sliding_window_results
  fwrite(sliding_CpG, paste0("KY962518_5ETS_TO_3ETS_CpG_dinucleotides_sliding_data_", i,"bp.csv"))
  
  
  sliding_CpG_graph<- ggplot(sliding_CpG, aes(x = start, y = CpG_OE)) +
    geom_line(color = "#E21515") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0.65, linetype = "dashed", color="orange")+
    labs(title = "CpG dinucleotides O/E ratio across rDNA Sequence KY962518",
         subtitle= paste0(i, " bp sliding window size"),
         x = "Position (bp)",
         y = "CpG dinucleotides O/E ratio") +
    scale_x_continuous(breaks = c(0, 2202, 5859, 7728,8798, 8955, 10122, 15173, 15534), 
                       labels =c("Promoter","5'ETS", "18S","ITS1", "5.8S", "ITS2", "28S", "3'ETS", "IGS"))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20))
  
  ggsave(paste0("KY962518_5ETS_TO_3ETS_sliding_CpG_graph_", i, "bp.tiff"), 
         plot =  sliding_CpG_graph, width=18, height = 10, dpi = 150)
  
  
  
}  






# if you want to run non-sliding skew
non_sliding_skew<- gc_skew_data$fixed_window_results
fwrite(non_sliding_skew, paste0("KY962518_5ETS_TO_3ETS_gc_skew_non_sliding_data_", i,"bp.csv"))

  non_sliding_graph<- ggplot(non_sliding_skew, aes(x = start, y = GC_skew_value)) +
    geom_line(color = "#E21515") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "GC Skew Across rDNA Sequence KY962518",
         subtitle= paste0(i, " bp non sliding window size"),
         x = "Position (bp)",
        y = "GC Skew (G−C / G+C)") +
    scale_x_continuous(breaks = c(0, 3657,5526, 6596, 6753, 7920, 12971, 13333), 
                       labels =c("5'ETS", "18S","ITS1", "5.8S", "ITS2", "28S", "3'ETS", "IGS"))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20))
  
  ggsave(paste0("KY962518_5ETS_TO_3ETS_gc_skew_non_sliding_", i, "bp.tiff"), 
         plot = non_sliding_graph, width=18, height = 10, dpi = 150)
  
  non_sliding<- paste0("KY962518_5ETS_TO_3ETS_gc_content_non_sliding_data_", i,"bp.csv")
  non_sliding_content<- fread(non_sliding, sep = ",", header = TRUE)
  
  
  #if you want to re-run only GC content thing
  
  #for ( i in bin_size){
    
    
  
  #if you want to run non sliding GC content
  non_sliding_content<- gc_content_data$fixed_window_results
  fwrite(non_sliding_content, paste0("KY962518_5ETS_TO_3ETS_gc_content_non_sliding_data_", i,"bp.csv"))
  
  non_sliding_content_graph<- ggplot(non_sliding_content, aes(x = start, y = gc_content_perc)) +
    geom_line(color = "#E21515") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 40.89, linetype = "dashed", color="darkgreen")+
    labs(title = "GC Content Across rDNA Sequence KY962518",
         subtitle= paste0(i, " bp non sliding window size"),
         x = "Position (bp)",
        y = "GC Fraction (G+C / A+T+G+C)") +
    scale_y_continuous(limits= c(0,100), breaks=seq(0,100, by=20))+
    scale_x_continuous(breaks = c(0, 3657,5526, 6596, 6753, 7920, 12971, 13333), 
                       labels =c("5'ETS", "18S","ITS1", "5.8S", "ITS2", "28S", "3'ETS", "IGS"))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
          axis.ticks.y = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=20))
  
  ggsave(paste0("KY962518_5ETS_TO_3ETS_gc_content_non_sliding_", i, "bp.tiff"), 
         plot = non_sliding_content_graph, width=18, height = 10, dpi = 150)
  


  #if you want to visualize 
  library(karyoploteR)
  
  custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=48338)) #entire region rDNA including IGS
  kp <- plotKaryotype(genome=custom_genome, plot.type = 2)
  
  kpRect(kp, chr = 'rDNA_locus', x0 = 1, x1 =1298 , y0 = 0, y1 = 1, col = "#DE9A22", data.panel = "ideogram", borders= NA) #marks last 1298bp from IGS representing previous rdna 
  kpRect(kp, chr = 'rDNA_locus', x0 = 1299, x1 =3500 , y0 = 0, y1 = 1, col = "#F5FEFB", data.panel = "ideogram", borders= NA) #marks 2202 bp of  promoter
  
  kpRect(kp, chr = 'rDNA_locus', x0 = 3501, x1 = 7157 , y0 = 0, y1 = 1, col = "#E21515", data.panel = "ideogram", borders= NA) #marks 5'ETS (3501+(3657-1))
  #3501+(3657-1) = 7157
  
  kpRect(kp, chr = 'rDNA_locus', x0 = 7158, x1 = 9026, y0 = 0, y1 = 1, col = "#5AAA46", data.panel = "ideogram", borders= NA) #marks 18S
  #7158+ (1869-1) = 9026
  
  kpRect(kp, chr = 'rDNA_locus', x0 = 9027, x1 = 10096, y0 = 0, y1 = 1, col = "#F36017", data.panel = "ideogram", borders= NA) #marks ITS1S
  #9027+ (1070-1) = 10096 
  
  kpRect(kp, chr = 'rDNA_locus', x0 = 10097, x1 = 10253, y0 = 0, y1 = 1, col = "#6B1519", data.panel = "ideogram", borders= NA) #marks 5.8S
  #10097+ (157-1) = 10253
  
  kpRect(kp, chr = 'rDNA_locus', x0 = 10254, x1 = 11420, y0 = 0, y1 = 1, col = "#818689", data.panel = "ideogram", borders= NA) #marks ITS2
  #10254+(1167-1) = 11420
  
  kpRect(kp, chr = 'rDNA_locus', x0 = 11421, x1 = 16471, y0 = 0, y1 = 1, col = "#ECE612", data.panel = "ideogram", borders= NA) #marks 28S
  #11421+(5051-1) = 16471
  
  kpRect(kp, chr = 'rDNA_locus', x0 = 16472, x1 = 16832, y0 = 0, y1 = 1, col = "#E07F80", data.panel = "ideogram", borders= NA) #marks 3'ETS
  #16472+(361-1) = 16832
  
  kpRect(kp, chr = 'rDNA_locus', x0 = 16833, x1 = 48338, y0 = 0, y1 = 1, col = "#DE9A22", data.panel = "ideogram", borders= NA) #marks IGS
  #16833+ (31506-1)= 48338
  
  
  kpRect(kp, chr = 'rDNA_locus', x0 = 46137, x1 = 48338, y0 = 0, y1 = 1, col = "#F5FEFB", data.panel = "ideogram", borders= NA) #marks 2202bp of promoter
  #48338-2201
  
  
  kpLines(kp, chr = "rDNA_locus", x=gc_skew_50bp_data$start , y = gc_skew_50bp_data$GC_skew_value, col= "blue", lwd=2.0, r0= -1.2, r1=-2.0)
  
  
  



