# Author: Jyoti Devendra Adala under supervision of Dr. Bruce Knutson
# For updates and contributions, visit : https://github.com/adalaj


#Calculating CpG content function

#this function simply calculates CpG nuclotides, their percent and their o/e ratio 
# below info is taken form UCSC genome browser https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=cpgIslandSuper
#The CpG count is the number of CG dinucleotides in the island. 
#The Percentage CpG is the ratio of CpG nucleotide bases (twice the CpG count) to the length. 
#The ratio of observed to expected CpG is calculated according to the formula (cited in Gardiner-Garden et al. (1987)):
  
  #Obs/Exp CpG = Number of CpG * N / (Number of C * Number of G)
  # where N = length of sequence.
  

#1) if window size is not provided, it calculated CpG content for the entire sequence
#2) If window size is provided, it calculated CpG content using fixed (Non overlapping) and sliding windows (overlapping window).

#non overlapping meaning bin size could be 1 to 10 then 11 to 20 
# overlapping meaning bin size could be 1 to 10, 2 to 11 and so on. 


#CpG dinucleotides 
# as per Gardiner-Garden and Frommer: paper title : CpG Islands in vertebrate genomes (1987) teh CpG has to have a length of 200, GC >50%, and CpG O/E >0.6
# but in 2002 a new paper published takai and Jones (Comprehensive analysis of CpG islands in human chromosomes 21 and 22)
# they said that earlier formula was made when sequencing was not performed, contain repetitive and alu elements. 
# This demanded new formula. as per which length of 500bp, GC% >55%, and CpG O/E >0.6is needed.


library(stringr)

CpG_content <- function(seq, window_size= NULL){
  seq_length<- nchar(seq)
  
  # If window_size is NULL, compute CpG content for the entire sequence
  if (is.null(window_size)){
    g_count <- str_count(seq, "G")
    c_count <- str_count(seq, "C")
    CpG_count <- str_count(seq, "CG")
    

    if (CpG_count>0){
      CpG_perc = (CpG_count/seq_length)*100
      CpG_OE = (CpG_count * seq_length) / (g_count*c_count)

    }else{
      CpG_count <- 0
    }
    
    return(data.frame(
      start = 1, 
      end = seq_length,
      window_seq = seq,
      G_count = g_count,
      C_count = c_count,
      CpG_count = CpG_count,
      CpG_perc= CpG_perc,
      CpG_OE = CpG_OE,
      stringsAsFactors = FALSE
    ))
  }
  
  # if window_size is provided, compute CpG content using a fix window size slide
  #such as 1 to 10 then 11 to 20 etc 
  
  positions<- seq(1, seq_length-window_size+1, by=window_size)
  fixed_window_results <- data.frame(
    start = numeric(), 
    end = numeric(), 
    window_seq = character(), 
    G_count = numeric(), 
    C_count = numeric(), 
    CpG_count = numeric(),
    CpG_perc= numeric(),
    CpG_OE = numeric(),
    stringsAsFactors = FALSE  # Avoid automatic conversion to factors
  )
  
  for (i in seq_along(positions)){
    window_seq<- str_sub(seq, positions[i], positions[i]+window_size-1)
    g_count<- str_count(window_seq, "G")
    c_count<- str_count(window_seq, "C")
    CpG_count <- str_count(window_seq, "CG")
    
    
    if (CpG_count>0){
      CpG_perc = (CpG_count/window_size)*100
      CpG_OE = (CpG_count * window_size) / (g_count*c_count)

    }else{
      CpG_count <- 0
    }
    
    
    fixed_window_results<- rbind(fixed_window_results, data.frame(
      start = positions[i], 
      end= positions[i]+window_size-1,
      window_seq = window_seq,
      G_count= g_count,
      C_count= c_count,
      CpG_count = CpG_count,
      CpG_perc= CpG_perc,
      CpG_OE = CpG_OE))
    
  }
  
  #
  sliding_window_positions <- seq(1, seq_length-window_size+1)
  sliding_window_results <- data.frame(
    start = numeric(), 
    end = numeric(), 
    window_seq = character(), 
    G_count = numeric(), 
    C_count = numeric(), 
    CpG_count = numeric(),
    CpG_perc= numeric(),
    CpG_OE = numeric(),
    stringsAsFactors = FALSE  # Avoid automatic conversion to factors
  )
  
  for (i in seq_along(sliding_window_positions)){
    window_seq<- str_sub(seq, sliding_window_positions[i], sliding_window_positions[i]+window_size-1)
    g_count<- str_count(window_seq, "G")
    c_count<- str_count(window_seq, "C")
    CpG_count <- str_count(window_seq, "CG")
    
    
    if (CpG_count>0){
      CpG_perc = (CpG_count/window_size)*100
      CpG_OE = (CpG_count * window_size) / (g_count*c_count)

    }else{
      CpG_count <- 0
    }
    
    
    sliding_window_results<- rbind(sliding_window_results, data.frame(
      start = sliding_window_positions[i], 
      end= sliding_window_positions[i]+window_size-1,
      window_seq = window_seq,
      G_count= g_count,
      C_count= c_count,
      CpG_count = CpG_count,
      CpG_perc= CpG_perc,
      CpG_OE = CpG_OE))
    
  }
  
  return(list(sliding_window_results = sliding_window_results,fixed_window_results=fixed_window_results))
  
}


