function(MAT,SEQ,DRIVER,WEIGHT){
  
  scores = mapply(pwmScore_fun, id = DRIVER[,2], tf = DRIVER[,1], 
                  MoreArgs = list(pwmMat = MasterPWM, seq=SEQ, WT = 0))
  
  weights = mapply(pwmScore_fun, id = DRIVER[,2], tf = DRIVER[,1], 
                  MoreArgs = list(pwmMat = MasterPWM, seq=SEQ, WT = 1))
  

  wtScores = weights*scores
  
  wtSum <- as.data.frame(cbind(DRIVER[,1],t(weights)))%>%
    group_by(V1)%>%
    summarize_each(funs(sum))%>%
    ungroup()%>%
    select(-V1) 
  
  scSum <- as.data.frame(cbind(DRIVER[,1],t(wtScores)))%>%
    group_by(V1)%>%
    summarize_each(funs(sum))%>%
    ungroup()%>%
    select(-V1)
  
  scMean = scSum/wtSum
  
  if (WEIGHT==1){
    return(wtSum)
  }
  else{
    return(scMean)
  }
}
