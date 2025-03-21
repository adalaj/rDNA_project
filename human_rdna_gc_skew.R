

#Task is to create GC skew for human RDNA locus 
# and make coorelation with RLFS.

#first make GC for each region and then calculate correlation with RLFS in that region.
#next, I wanted to take average RLFS length and then calculate GC skew using sliding window.


# for first analysis, I already defined a function to calculate GC skew  

#load GC skew function
source("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/script/gc_skew_function.R")
rdna_human<- fread ("rdna_hg38_chr21_2018_dataset_details_v2.csv", sep = ",", header = TRUE)


library(stringr)
library(tidyverse)
library(data.table)

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/rloop_and_rdna/human/one_rDNA_seq/output")

attempt1<- data.frame((matrix(nrow = 0, ncol=7)))
                      
for ( i in 1:nrow(rdna_human)){
  gc_skew_data<- gc_skew(rdna_human$Sequences[i])
  attempt1<- rbind(attempt1, gc_skew_data)
}



selected_gc_skew_data <- attempt1 %>% select(GC_skew_value, norm_GC_skew) 
rdna_human<- cbind(rdna_human, selected_gc_skew_data)
fwrite(rdna_human, "rdna_hg38_chr21_2018_dataset_details_v3.csv")

full_length_gc_skew<- ggplot(rdna_human, aes(x = x_axis, y = norm_GC_skew)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
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
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "GC Skew Across rDNA Sequence KY962518",
       x = "rDNA",
       y = "Normalised GC Skew") +
  scale_x_continuous(breaks = c(2202, 5859, 7728,8798, 8955, 10122, 15173, 15534 ),
                     labels = c("promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS"))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=20))

ggsave("KY962518_full_length_excld_IGS_normalised_gc_skew.tiff", 
       plot = full_length_gc_skew_excld_igs, width=18, height = 10, dpi = 150)


# The resulst are not that informative, so decided to to GC skew in 100 bp sliding window

library(Biostrings)

rdna_human_seq<- readDNAStringSet("KY962518_added_3500nt_IGS_upstream_nontemplate.fasta")
seq<- as.character(rdna_human_seq[[1]])

gc_skew_100bp_data <- gc_skew(seq, window_size = 100)
fwrite(gc_skew_100bp_data, "KY962518_added_3500nt_IGS_upstream_nontemplate_gc_skew_data_100bp.csv", sep=",")

gc_skew_50bp_data <- gc_skew(seq, window_size = 50)
fwrite(gc_skew_50bp_data, "KY962518_added_3500nt_IGS_upstream_nontemplate_gc_skew_data_50bp.csv", sep=",")

gc_skew_10bp_data <- gc_skew(seq, window_size = 10)
fwrite(gc_skew_10bp_data, "KY962518_added_3500nt_IGS_upstream_nontemplate_gc_skew_data_10bp.csv", sep=",")








gc100<- ggplot(gc_skew_100bp_data, aes(x = start, y = GC_skew_value)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "GC Skew Across rDNA Sequence KY962518",
       subtitle= "100 bp Sliding Window",
       x = "Position (bp)",
       y = "GC Skew ") +
  scale_x_continuous(breaks = seq(0, 50000, by = 5000))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=20))

ggsave("KY962518_added_3500nt_IGS_upstream_gc_skew_100bp.tiff", 
       plot = gc100, width=18, height = 10, dpi = 150)

gc50<- ggplot(gc_skew_50bp_data, aes(x = start, y = GC_skew_value)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "GC Skew Across rDNA Sequence KY962518",
       subtitle= "50 bp Sliding Window",
       x = "Position (bp)",
       y = "GC Skew ") +
  scale_x_continuous(breaks = seq(0, 50000, by = 5000))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=20))

ggsave("KY962518_added_3500nt_IGS_upstream_gc_skew_50bp.tiff", 
       plot = gc50, width=18, height = 10, dpi = 150)


gc10<- ggplot(gc_skew_10bp_data, aes(x = start, y = GC_skew_value)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "GC Skew Across rDNA Sequence KY962518",
       subtitle= "10 bp Sliding Window",
       x = "Position (bp)",
       y = "GC Skew ") +
  scale_x_continuous(breaks = seq(0, 50000, by = 5000))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=20))

ggsave("KY962518_added_3500nt_IGS_upstream_gc_skew_10bp.tiff", 
       plot = gc10, width=18, height = 10, dpi = 150)


library(karyoploteR)

custom_genome <- toGRanges(data.frame(chr="rDNA_locus", start=1, end=48338))
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



