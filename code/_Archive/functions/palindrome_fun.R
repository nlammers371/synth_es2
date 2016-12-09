function(seq, len){ 
  #Function to find non-overlapping palindrome pairs of a specified length
  #Variables:
  #   seq (CharVec):  DNA sequence to be searched
  #   length (Int): Length of palindrome seq
  
  
  seqTruncbp = strsplit(seq,NULL)[[1]]
  seqL = nchar(seq)
  #reverse adjusted sequence
  seqRev <- paste0(rev(seqTruncbp), collapse = '')
  seqRev <- gsub('F','R',seqRev)
  seqRev <- substr(seqRev,1,floor(seqL/2))
  
  palF <- vector()
  palR <- vector()
  
  for (i in 1:(seqL-nchar(seqRev) -len + 1)){
    lead <- i + len - 1
    f <- substr(seq,i,lead)
    
    Rmatch <- as.vector(gregexpr(f,seqRev)[[1]])
    if (Rmatch[1] > 0){
      Rmatch <- seqL - Rmatch - len + 2
      Fmatch <- i  
      palR <- append(palR,Rmatch)
      palF <- append(palF,Fmatch)
    }
  }
  
  palRseq <- vector()
  palFseq <- vector()
  for (i in 1:length(palR)){
    Rev <- palR[i] 
    For <- palF[i]
    palFseq[i] <- substr(seq,For,For + len - 1)
    palRseq[i] <- substr(seq,Rev,Rev + len - 1)
  }
  palindromes <- cbind(palF,palFseq,palR,palRseq)
  return(palindromes)
}
