# ------------------------------------------------------------------------------
# This code is part of paper: In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus.
# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj
#
# Purpose:
#   Compute GC skew across a DNA sequence, defined as:
#       GC_skew = (G - C) / (G + C)
#   GC skew indicates strand compositional bias and can reflect replication
#   origin or transcriptional asymmetry.
#
#   If no window size is provided, GC skew is computed for the entire sequence.
#   If a window size is provided, GC skew is computed for:
#     (a) fixed, non-overlapping windows (1–10, 11–20, …)
#     (b) sliding, overlapping windows (1–10, 2–11, …)
#
# Input:
#   - seq: A character string representing a DNA sequence.
#   - window_size: (optional) Integer; window length for skew calculation.
#
# Output:
#   - If window_size = NULL → returns a data.frame with overall GC skew value.
#   - If window_size provided → returns a list with two data.frames:
#       * fixed_window_results
#       * sliding_window_results
#
# Notes:
#   - Case-insensitive ('g' and 'c' treated as 'G' and 'C').
#   - Non-ACGT characters are allowed; they count toward window length but not
#     toward GC skew numerator/denominator.
#   - Function operates in memory (no input/output files).
# ------------------------------------------------------------------------------

#load library
library(stringr)

gc_skew <- function(seq, window_size= NULL){
  seq_length<- nchar(seq)
  
  # If window_size is NULL, compute GC skew for the entire sequence
  if (is.null(window_size)){
    g_count <- str_count(seq, "G")
    c_count <- str_count(seq, "C")
    
    if (g_count+c_count>0){
      skew_value <- (g_count-c_count)/ (g_count+c_count)
    }else{
      skew_value <- 0
    }
    

    return(data.frame(
      start = 1, 
      end = seq_length,
      window_seq = seq,
      G_count = g_count,
      C_count = c_count,
      GC_skew_value = skew_value,
      stringsAsFactors = FALSE
    ))
  }
  
  # if window_size is provided, compute GC skew using a fix window size slide
  #such as 1 to 10 then 11 to 20 etc 
  
  positions<- seq(1, seq_length-window_size+1, by=window_size)
  fixed_window_results <- data.frame(
    start = numeric(), 
    end = numeric(), 
    window_seq = character(), 
    G_count = numeric(), 
    C_count = numeric(), 
    GC_skew_value = numeric(),
    stringsAsFactors = FALSE  # Avoid automatic conversion to factors
  )
  
  for (i in seq_along(positions)){
    window_seq<- str_sub(seq, positions[i], positions[i]+window_size-1)
    g_count<- str_count(window_seq, "G")
    c_count<- str_count(window_seq, "C")
    
    if (g_count+c_count>0){
      skew_value <- (g_count-c_count)/ (g_count+c_count)
    }else{
      skew_value <- 0
    }
  
  
    fixed_window_results<- rbind(fixed_window_results, data.frame(
                    start = positions[i], 
                    end= positions[i]+window_size-1,
                    window_seq = window_seq,
                    G_count= g_count,
                    C_count= c_count,
                    GC_skew_value= skew_value))
    
  }
  
  #
  sliding_window_positions <- seq(1, seq_length-window_size+1)
  sliding_window_results <- data.frame(
    start = numeric(), 
    end = numeric(), 
    window_seq = character(), 
    G_count = numeric(), 
    C_count = numeric(), 
    GC_skew_value = numeric(),
    stringsAsFactors = FALSE  # Avoid automatic conversion to factors
  )
  
  for (i in seq_along(sliding_window_positions)){
    window_seq<- str_sub(seq, sliding_window_positions[i], sliding_window_positions[i]+window_size-1)
    g_count<- str_count(window_seq, "G")
    c_count<- str_count(window_seq, "C")
    
    if (g_count+c_count>0){
      skew_value <- (g_count-c_count)/ (g_count+c_count)
    }else{
      skew_value <- 0
    }
    
  
  sliding_window_results<- rbind(sliding_window_results, data.frame(
    start = sliding_window_positions[i], 
    end= sliding_window_positions[i]+window_size-1,
    window_seq = window_seq,
    G_count= g_count,
    C_count= c_count,
    GC_skew_value= skew_value))
  
  }
  
  return(list(sliding_window_results = sliding_window_results,fixed_window_results=fixed_window_results))
  
  

}


