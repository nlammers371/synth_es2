function(seq,ref,subVec, nsub,tol){
  seqL <- length(subVec)
  dist <- stringDist(c(paste0(seq),paste0(ref)), method = "hamming",
             ignoreCase = FALSE, diag = FALSE, upper = FALSE)[1]
  
  frac <- 1-(dist)/(seqL)
  wt <- max((tol-frac)/(4*tol),0)
  probs <- append(wt,rep(1,3)/(1-tol))
  subInd <- sample(subVec,nsub,replace=FALSE)
  subChar <- vector()
  NewSeq <- seq
  for (i in 1:nsub){
    NewSeq[i] <- sample(sampMat[subInd[i],],1,prob=probs)
  }
  return(NewSeq)
}