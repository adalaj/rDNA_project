# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj

#I wanted to make bar graph for GC percent of different regions of genome. 
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/gc_content_human")


library(tidyverse)
library(data.table)

gc_ranges <- data.frame(
  Category = c("Human genome (Overall)", 
               "Protein-coding Exons", 
               "Introns", 
               "CDS (coding sequence)", 
               "5'UTR", 
               "3'UTR", 
               "lncRNA", 
               "rDNA (coding region)"),
  Min = c(40, 50, 35, 52, 60, 42,35, 70),
  Max = c(41, 55, 45, 60, 70, 53,47, 71)
)

gc_ranges$Category<- factor(gc_ranges$Category,
                             levels = c("Human genome (Overall)", "Protein-coding Exons",
                                        "CDS (coding sequence)", "5'UTR", "3'UTR",
                                        "Introns", "lncRNA","rDNA (coding region)"))

gc_rdna<- ggplot(gc_ranges, aes(x = Category, ymin = Min, ymax = Max)) +
  geom_errorbar(width = 0.2, color = "darkblue", size = 1.2) +
  geom_point(aes(y = (Min + Max)/2), size = 4, color = "red") +
  scale_y_continuous(limits =c(0,80), breaks=seq(0,80, by=20)) +
  labs(title = "GC content Comparison:rDNA vs Others",
       x = "Category",
       y = "GC content (%)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 25), # Center Y-axis title
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1, size=20))+
  coord_flip() #you can remove this if you dont want then flipped ones then



ggsave("GC_content_comparison_rdna_vs_others_flipped.tiff", 
       plot =  gc_rdna, width=18, height = 10, dpi = 150)

#bruce also asked to make x axis in reverse to represent error bar


#bruce has not looked into this
gc_bar<- ggplot(gc_ranges, aes(x = Category, y= (Min+Max)/2)) +
  geom_col(fill = "skyblue", color = "darkblue", width = 0.6, size =2) +  # Bar plot
  geom_errorbar(aes(ymin=Min, ymax=Max), width = 0.2, color = "darkblue", size = 1) +
  geom_point(aes(y = (Min + Max)/2), size = 4, color = "red") +
  scale_y_continuous(limits =c(0,80), breaks=seq(0,80, by=20), expand = expansion(mult = c(0,0.03))) +
  labs(#title = "GC content Comparison:rDNA vs Others",
       x = "Category",
       y = "GC content (%)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 100),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 6),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5,  margin = margin(r=20)), # Center Y-axis title
        axis.ticks = element_line(color = "black", linewidth = 4),
        axis.ticks.length = unit(30, "pt"),
        axis.text.x = element_text(hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"))+
  coord_flip() #you can remove this if you dont want then flipped ones then

  
ggsave("GC_content_comparison_rdna_vs_others_bar_flipped.png", 
       plot =  gc_bar, width=35, height = 20, dpi = 150)


