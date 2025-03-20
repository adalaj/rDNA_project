##define GC skew function

library(stringr)

gc_skew <- function(seq, window_size=10){
  seq_length<- nchar(seq)
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
      skew_values <- (g_count-c_count)/ (g_count+c_count)
    }else{
      skew_values <- 0
    }
  
  
  results<- rbind(results, data.frame(
                    start = positions[i], 
                    end= positions[i]+window_size-1,
                    window_seq = window_seq,
                    G_count= g_count,
                    C_count= c_count,
                    GC_skew_value= skew_values))
  }
  
  return(results)
  
}




