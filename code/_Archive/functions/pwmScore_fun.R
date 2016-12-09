function(pwmMat=MasterPWM, tf, id, seq=es2,complement){  
  
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
  #sum pwm score contributions from 2 "views"
  if (complement==0){
  RawScore <- (PWMscoreStartingAt(pwm, seq, 1:lag)  + 
              rev(PWMscoreStartingAt(pwm, Rseq, 1:lag)))
  }              
  if (complement==1){
    RawScore <- (PWMscoreStartingAt(pwm, Cseq, 1:lag)+  
                rev(PWMscoreStartingAt(pwm, RCseq, 1:lag)))
  }
  #define min score as zero point. Normalize by max value 
  NormScore <- (RawScore-2*mScore)/(2*norm)
  
  #weight vector to track share of pwm score received by each cell
  wt <- rep(.5, length(RawScore)) 
  
  weights <- Zseq + c(f1,wt,t1) + c(f2,wt,t2)
  ShiftScore <- Zseq + c(f1,NormScore,t1)/2 + c(f2,NormScore,t2)/2
  ShiftScore[weights!=1] = NA
  
  return(ShiftScore)
  
}