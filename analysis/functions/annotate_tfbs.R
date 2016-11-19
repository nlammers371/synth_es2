function(seq,tfbs){
  seqbp <- strsplit(seq, NULL)[[1]]
  seqL <- length(seqbp)
  
  #Remove known binding sites and create adjustment index to track gaps
  seqTruncbp <- seqbp
  for (i in 1:nrow(tfbs)){
    s = tfbs[i,1]
    e = tfbs[i,2]
    seqTruncbp[s:e] <- 'F'
  }
  
  seqTrunc <- ""
  for (i in 1:seqL){
    seqTrunc = paste0(seqTrunc,seqTruncbp[i])
  }
  return(seqTrunc)
}