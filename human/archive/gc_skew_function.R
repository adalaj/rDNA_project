
#GC skew function
##define GC skew function based on a window size. If window size is not defined then
## this function will compute the function based on GC skew of the entire given sequence. 
## Additionally, for some analysis, we are calculating GC skew for multiple sequence of varying length
## As entire length sequence can sometime vary, so this function will also calculate the 
## normalized GC skew based on sequence length. 



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
    
    # Normalized GC skew
    norm_GC_skew = skew_value/ seq_length
    
    return(data.frame(
      start = 1, 
      end = seq_length,
      window_seq = seq,
      G_count = g_count,
      C_count = c_count,
      GC_skew_value = skew_value,
      norm_GC_skew = norm_GC_skew, # Add normalized GC skew
      stringsAsFactors = FALSE
    ))
  }
  
  # if window_size is provided, compute GC skew using a sliding window
  
  positions<- seq(1, seq_length-window_size+1, by=window_size)
  results <- data.frame(
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
  
  
  results<- rbind(results, data.frame(
                    start = positions[i], 
                    end= positions[i]+window_size-1,
                    window_seq = window_seq,
                    G_count= g_count,
                    C_count= c_count,
                    GC_skew_value= skew_value))
  }
  
  return(results)
  
}




