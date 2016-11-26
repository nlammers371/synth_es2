function(SEQ,DRIVER){
  scores <- vector()

  #I'm not sure why this function is incompatible with mapply format. revisit
  for (i in 1:nrow(DRIVER)){
    scores <- rbind(scores,pwmScore_fun(pwmMat=MasterPWM,seq=SEQ,wt=0,id=DRIVER[i,2],tf=DRIVER[i,1]))
  }
  scMean <- as.data.frame(cbind(DRIVER[,1],scores))%>%
    group_by(V1)%>%
    summarize_each(funs(mean))%>%
    ungroup()%>%
    select(-V1)
  
  return(scMean)
}
