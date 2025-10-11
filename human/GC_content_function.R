# ------------------------------------------------------------------------------
# This code is part of paper: In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus.
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
#   Compute GC content—the fraction/percent of 'G' and 'C' bases—in a DNA
#   sequence. If no window size is provided, GC is computed for the full
#   sequence. If a window size is provided, GC is computed for:
#     (a) fixed, non-overlapping windows (1–10, 11–20, …)
#     (b) sliding, overlapping windows (1–10, 2–11, …)
#
# Input:
#   - seq: A character string representing a DNA sequence.
#   - window_size: (optional) Integer; window length for GC calculation.
#
# Output:
#   - If window_size = NULL → returns a data.frame with overall GC content.
#   - If window_size provided → returns a list with two data.frames:
#       * fixed_window_results
#       * sliding_window_results
#
# Notes:
#   - Function operates in memory (no input/output files).
#   - Case-insensitive; 'g' and 'c' treated as 'G' and 'C'.
#   - Non-ACGT characters are ignored in GC counting.
# ------------------------------------------------------------------------------


#load library
library(stringr)

gc_content <- function(seq, window_size= NULL){
  seq_length<- nchar(seq)
  
  # If window_size is NULL, compute GC content for the entire sequence
  if (is.null(window_size)){
    g_count <- str_count(seq, "G")
    c_count <- str_count(seq, "C")
    
    if (g_count+c_count>0){
      gc_content_value <- (g_count+c_count)/ seq_length
    }else{
      gc_content_value <- 0
    }
    
    
    
    return(data.frame(
      start = 1, 
      end = seq_length,
      window_seq = seq,
      G_count = g_count,
      C_count = c_count,
      gc_content_value = gc_content_value,
      gc_content_perc = gc_content_value*100, # Add normalized GC content
      stringsAsFactors = FALSE
    ))
  }
  
  # if window_size is provided, compute GC content using a fix window size slide
  #such as 1 to 10 then 11 to 20 etc 
  
  positions<- seq(1, seq_length-window_size+1, by=window_size)
  fixed_window_results <- data.frame(
    start = numeric(), 
    end = numeric(), 
    window_seq = character(), 
    G_count = numeric(), 
    C_count = numeric(), 
    gc_content_value = numeric(),
    gc_content_perc = numeric(),
    stringsAsFactors = FALSE  # Avoid automatic conversion to factors
  )
  
  for (i in seq_along(positions)){
    window_seq<- str_sub(seq, positions[i], positions[i]+window_size-1)
    g_count<- str_count(window_seq, "G")
    c_count<- str_count(window_seq, "C")
    
    if (g_count+c_count>0){
      gc_content_value <- (g_count+c_count)/ window_size
    }else{
      gc_content_value <- 0
    }
    
    
    fixed_window_results<- rbind(fixed_window_results, data.frame(
      start = positions[i], 
      end= positions[i]+window_size-1,
      window_seq = window_seq,
      G_count= g_count,
      C_count= c_count,
      gc_content_value= gc_content_value,
      gc_content_perc = gc_content_value*100))
    
  }
  
  #
  sliding_window_positions <- seq(1, seq_length-window_size+1)
  sliding_window_results <- data.frame(
    start = numeric(), 
    end = numeric(), 
    window_seq = character(), 
    G_count = numeric(), 
    C_count = numeric(), 
    gc_content_value = numeric(),
    gc_content_perc = numeric(),
    stringsAsFactors = FALSE  # Avoid automatic conversion to factors
  )
  
  for (i in seq_along(sliding_window_positions)){
    window_seq<- str_sub(seq, sliding_window_positions[i], sliding_window_positions[i]+window_size-1)
    g_count<- str_count(window_seq, "G")
    c_count<- str_count(window_seq, "C")
    
    if (g_count+c_count>0){
      gc_content_value <- (g_count+c_count)/ window_size
    }else{
      gc_content_value <- 0
    }
    
    
    sliding_window_results<- rbind(sliding_window_results, data.frame(
      start = sliding_window_positions[i], 
      end= sliding_window_positions[i]+window_size-1,
      window_seq = window_seq,
      G_count= g_count,
      C_count= c_count,
      gc_content_value = gc_content_value,
      gc_content_perc = gc_content_value*100))
    
  }
  
  return(list(sliding_window_results = sliding_window_results,fixed_window_results=fixed_window_results))
  
  
  
}


