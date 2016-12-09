function(MAT,SEQ,DRIVER,COMP){
  
  scores = mapply(pwmScore_fun, id = DRIVER[,2], tf = DRIVER[,1], 
                  MoreArgs = list(pwmMat = MAT, seq=SEQ,complement=COMP))
  

  scMean <- as.data.frame(cbind(DRIVER[,1],t(scores)))%>%
    group_by(V1)%>%
    summarize_each(funs(mean))%>%
    select(-V1)
  
  return(scMean)
  
}
