function(SEQ,DRIVER){
  weights <- vector()
  #I'm not sure why this function is incompatible with mapply format. revisit
  for (i in 1:nrow(DRIVER)){
    weights <- rbind(weights,pwmScore_fun(pwmMat=MasterPWM,seq=SEQ,wt=1,id=DRIVER[i,2],tf=DRIVER[i,1]))  
  }
  wtMean <- as.data.frame(cbind(DRIVER[,1],weights))%>%
    group_by(V1)%>%
    summarize_each(funs(mean))%>%
    ungroup()%>%
    select(-V1)
  return(wtMean)
}
