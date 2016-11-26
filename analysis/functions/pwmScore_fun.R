function(pwmMat=MasterPWM, tf, id, seq=es2, wt=0){  
  seqL <- length(seq)
  Rseq <- rev(es2)
  Cseq <- complement(seq)
  RCseq <- reverseComplement(seq)
  
  PWM <- filter(pwmMat,
                TF == as.character(tf),
                ID == as.integer(id))%>%
    select(A,C,G,T)
  
  pwm <- as.matrix(t(PWM), row.names = c("A","C","G","T"))
  matL <- ncol(pwm)
  
  if (wt == 0){
    norm <- maxScore(pwm)
    lag <- seqL-matL+1
    Zpad <-as.vector(rep(0,matL-1))
    
    a = append(PWMscoreStartingAt(pwm, seq, 1:lag),Zpad)
    b = append(PWMscoreStartingAt(pwm, Rseq, 1:lag),Zpad)
    c = append(PWMscoreStartingAt(pwm, Cseq, 1:lag),Zpad)
    d = append(PWMscoreStartingAt(pwm, RCseq, 1:lag),Zpad)
    
    scoreF <- rowSums(cbind(a,c))/norm
    scoreR <- rowSums(cbind(b,d))/norm
    cF <- cumsum(scoreF)
    cR <- cumsum(scoreR)
    avgF <- vector()
    avgR <- vector()
    
    for (i in 1:(matL)){
      avgF[i] <- cF[i]/i
      avgR[i] <- cR[i]/i 
    }
    
    avgF[(matL+1):seqL] <- (cF[(matL+1):seqL] - cF[1:(seqL-matL)]) / matL
    avgR[(matL+1):seqL] <- (cR[(matL+1):seqL] - cR[1:(seqL-matL)]) / matL
    avgR <- rev(avgR)
    out = colSums(rbind(avgF,avgR))    
  }  
  if (wt == 1){
    wt <- vector()
    wt[1:seqL] <- matL
    wt[1:(matL-1)] <- 1:(matL-1)
    wt[(seqL-matL+2):seqL] <- (matL-1):1
    wtTot <- wt/matL
  
    out <- wtTot
  }
  return(out) 
}