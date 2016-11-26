function(alignment, wtVec, seqAnno, seq, lb, scores){

  seqL <- nchar(seqAnno)
  nSeq <- length(alignment)-1
  
  AlgnMat <- matrix(unlist(strsplit(paste0(alignment),NULL)), ncol = seqL, byrow = TRUE)
  AlgnMat <- AlgnMat[2:nrow(AlgnMat),]
  
  IndVec <- as.vector(gregexpr('-',seqAnno)[[1]])
  CharVec <- unlist(strsplit(seq,NULL)) 
  CharMat <- t(matrix(rep(CharVec,nSeq), nrow=seqL))
  
  wtMat <- matrix(rep(wtVec,seqL), nrow=length(wtVec))
  wtDenom <- colSums(wtMat)
  
  wtScore <- colSums((AlgnMat==CharMat)*wtMat)/wtDenom
  #wtScore[IndVec] <- 1.5
  es2 <- DNAString(seq)
  es2cons <- es2[wtScore>=lb]
  if (scores==1){
    return(wtScore)
  }
  else{
    return(es2cons)
  }
}