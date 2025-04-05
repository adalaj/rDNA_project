#Calculating GC content function

#this function simply calculates GC content that is the percentage of G and C in a given DNA sequence.

#1) if window size is not provided, it calculated GC content for the entire sequence
#2) If window size is provided, it calculated GC content using fixed (Non overlapping) and sliding windows (overlapping window).

#non overlapping meaning bin size could be 1 to 10 then 11 to 20 
# overlapping meaning bin size could be 1 to 10, 2 to 11 and so on. 


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


