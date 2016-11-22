function(seq,tfbs){
  seqbp <- strsplit(seq, NULL)[[1]]
  seqL <- length(seqbp)
  
  #Remove known binding sites
  seqTruncbp <- seqbp
  for (i in 1:nrow(tfbs)){
    s = tfbs[i,1]
    e = tfbs[i,2]
    seqTruncbp[s:e] <- 'F'
  }
  
  seqTrunc <- paste0(seqTruncbp,collapse = '')
  return(seqTrunc)
}