function(pwmMat=MasterPWM, tf, id, seq=es2, WT=0){  
  
  seqL <- length(seq)
  Rseq <- rev(es2)
  Cseq <- complement(seq)
  RCseq <- reverseComplement(seq)
  Zseq <- rep(0,seqL)
  
  PWM <- filter(pwmMat,
                TF == tf,
                ID == id)%>%
    select(A,C,G,T)
  
  pwm <- as.matrix(t(PWM), row.names = c("A","C","G","T"))
  matL <- ncol(pwm)
  med <- median(1:matL)
  f1 <- rep(0,floor(med)-1)
  f2 <- rep(0,ceiling(med)-1)
  t1 <- f2
  t2 <- f1
  
  mScore <- minScore(pwm)
  norm <- maxScore(pwm) - mScore
  lag <- seqL-matL+1
  #sum pwm score contributions from 4 "views"
  RawScore <- (PWMscoreStartingAt(pwm, seq, 1:lag)/2  + 
              PWMscoreStartingAt(pwm, Cseq, 1:lag)/2  +
          rev(PWMscoreStartingAt(pwm, Rseq, 1:lag)/2) +  
          rev(PWMscoreStartingAt(pwm, RCseq, 1:lag)/2))

  #define min score as zero point. Normalize by max value 
  NormScore <- (RawScore-4*mScore)/(4*norm)
  
  #weight vector to track share of pwm score received by each cell
  wt <- rep(.5, length(RawScore)) 
  
  ShiftScore <- Zseq + c(f1,NormScore,t1) + c(f2,NormScore,t2)
  weights <- Zseq + c(f1,wt,t1) + c(f2,wt,t2)
  wtNorm <- weights/max(weights)
  
  if (WT == 1){
    return(wtNorm)
  }
  else{
    return(ShiftScore)
  }
}