#bar graph representation of UBF and polr1A signal intensity in rdna sub regions.

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/rDNA/IGV/human/IGV_files/Bedgraphs/output_bigwig/chrR_regionwise_correlation")


library(tidyverse)
library(data.table)

#raw data from deeptools
entire<- fread("ChIP_regionwise_nontemplate_summary_raw_counts.tab", sep = "\t", header = TRUE) 

files_of_int <- list(entire = "ChIP_regionwise_nontemplate_summary_raw_counts.tab")#, i realised template and nontemplate is for RLFS, UBF and POLR1 will be same 
                     #template = "ChIP_regionwise_nontemplate_summary_raw_counts.tab")  # Add actual file names here

for (i in names(files_of_int)) {
  # Step 1: Read the file
  data <- fread(files_of_int[[i]], sep="\t", header = TRUE)
  
  # Step 2: Slice relevant rows and rename columns
  filt <- data[2:9, 1:5]  # first row is IGS, rows 2–9 are rDNA regions promoter to 3'ETS
  
  colnames(filt) <- c("chr","start", "end", "mean_UBF_intensity", "mean_POLR1A_intensity")
  
  # Step 3: Add region names
  filt$rdna_region <- c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS")
  
  # Step 4: Add region length
  filt <- filt %>%
    mutate(region_length = end - start + 1)
  
  # Step 5: Normalize columns by region length
  filt <- filt %>%
    mutate(across(c(mean_UBF_intensity, mean_POLR1A_intensity),
                  ~ . / region_length, #function(x) x / region_length
                  .names = "norm_{.col}"))
  
  # Optional: save or store output
  fwrite(filt, paste0("normalized_", i, ".csv"))  # Creates a new object in the environment
  
  
  filt$rdna_region <- factor(filt$rdna_region, 
                             levels = c("Promoter", "5'ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS"))
  
  plot_data <- filt %>%
    select(rdna_region, starts_with("norm_")) %>%
    pivot_longer(cols = starts_with("norm_"),
                 names_to = "Signal",
                 values_to = "Normalized_Intensity") 
  
  
 normalised_graph<- ggplot(plot_data, aes(x = rdna_region, y = Normalized_Intensity, fill = Signal)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("mean_UBF_intesity" = "#FF66CC", "mean_POLR1A_intesity" = "#39FF14"),
                      labels = c("mean_UBF_intesity" = "UBF", "mean_POLR1A_intesity"= "POLR1A")) + 
    labs(title = paste0("Normalized Signal across rDNA Regions_",i),
         x = "rDNA Region",
         y = "Normalized Signal (per bp)",
         fill = "Signal") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0(i, "_UBF_POL1_association_normalized_by_length.tiff"), normalised_graph, width = 30, height = 30, dpi = 150)
  
  plot_data_new <- filt %>%
    select(rdna_region, starts_with("mean_")) %>%
    pivot_longer(cols = starts_with("mean_"),
                 names_to = "Signal",
                 values_to = "Raw_Intensity") 
  
  
  raw_graph<- ggplot(plot_data_new, aes(x = rdna_region, y = Raw_Intensity, fill = Signal)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("mean_UBF_intesity" = "#FF66CC", "mean_POLR1A_intesity" = "#39FF14"),
                      labels = c("mean_UBF_intesity" = "UBF", "mean_POLR1A_intesity"= "POLR1A")) + 
    labs(title = paste0("Normalized Signal across rDNA Regions_",i),
         x = "rDNA Region",
         y = "Raw Signal Intensity (per bp)",
         fill = "Signal") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0(i,"_UBF_POL1_association_raw_intensity.tiff"), raw_graph, width = 30, height = 30, dpi = 150)
  
  
  
  
}

 
#bruce also wanted  
