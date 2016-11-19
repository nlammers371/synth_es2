function(seq, length){ 
  base <- length #set relevant "word" size
  seqTruncbp = strsplit(seq,NULL)[[1]]
  seqL = nchar(seq)
  #reverse adjusted sequence
  seqRevbp <- rev(seqTruncbp)
  seqRev = ""
  for (i in 1:seqL){
    seqRev = paste0(seqRev,seqRevbp[i])
  }
  seqRev <- gsub('F','R',seqRev)
  seqRev <- substr(seqRev,1,floor((seqL-base + 1)/2))
  
  #change to apply family once I have basic method down
  palF <- vector()
  palR <- vector()
  
  for (i in 1:floor((seqL-base + 1)/2)){
    lead <- i + base - 1
    f <- substr(seq,i,lead)
    
    Rmatch <- as.vector(gregexpr(f,seqRev)[[1]])
    if (Rmatch[1] > 0){
      Rmatch <- seqL - Rmatch - base + 2
      Fmatch <- (Rmatch==Rmatch)*i  
      palR <- append(palR,Rmatch)
      palF <- append(palF,Fmatch)
    }
  }
  
  palRseq <- vector()
  palFseq <- vector()
  for (i in 1:length(palR)){
    RevInd <- palR[i] 
    ForInd <- palF[i]
    palFseq[i] <- substr(seq,ForInd,ForInd + base - 1)
    palRseq[i] <- substr(seq,RevInd,RevInd + base - 1)
  }
  palindromes <- cbind(palF,palFseq,palR,palRseq)
  return(palindromes)
}
