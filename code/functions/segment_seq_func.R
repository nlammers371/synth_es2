function(seqDNA,seqCASE,bpUnit,minStretch){
  
  seqCharVec <- as.vector(unlist(strsplit(seqCASE,NULL)))
  
  boundsLU <- as.vector(unlist(gregexpr(
    paste0('([a,t,g,c]{',bpUnit,'})([A,T,G,C]{',bpUnit,'})'),seqCASE)))+bpUnit-1
  boundsUL <- as.vector(unlist(gregexpr(
    paste0('([A,T,G,C]{',bpUnit,'})([a,t,g,c]{',bpUnit,'})'),seqCASE)))+bpUnit-1
  boundsLU <- cbind(rep(1,length(boundsLU)),boundsLU)
  boundsUL <- cbind(rep(2,length(boundsUL)),boundsUL)
  
  bounds <- data.frame(rbind(boundsUL,boundsLU)) %>%
            arrange(boundsUL)%>%
            mutate(len = boundsUL-lag(boundsUL)) %>%
            filter(len >= bpUnit) %>% 
             #enforce min segment length
            mutate(nextID = ifelse(is.na(lead(V1)),0,lead(V1))) 
  
  oldID <- bounds[1,1]
  L <- bounds[,3]
  newInd <- vector()
  idCheck <- vector()
  ind <- 1
  for (i in 2:nrow(bounds)){
    ID <- bounds[i,1]
    if (ID != oldID && L >= minStretch){
      newInd[ind] <- bounds[i-1,2]
      idCheck[ind] <- oldID
      L <- bounds[i,3]
      ind <- ind + 1
      oldID <- ID
    }
    else {L <- L + bounds[i,3]}
  }
  
  newInd <- newInd[newInd<=(nchar(seqCASE)-bpUnit)] #ensure that end piece >= bpUnit
  
  stop <- append(newInd,nchar(seqCASE))
  start <- append(0,newInd)+1
  seq <- vector()
  
  for (i in 1:length(stop)){
    seq[i] <- paste0(seqCharVec[start[i]:stop[i]],collapse='')
  }
  ES2Segments <- data.frame(start,stop,seq)%>%
    mutate(len = nchar(as.character(seq)),
           left = substr(seq,1,bpUnit),
           right = substr(seq,len-bpUnit+1,len),
           L = lag(right),
           R = lead(left))%>%
    select(start,stop,seq,left,right,L,R)
  
  #ES2Segments <- ES2Segments[2:(nrow(ES2Segments)-1),]
  return(ES2Segments)
}