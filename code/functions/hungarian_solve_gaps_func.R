function(ref,IDvec){
  r <- nrow(ref)
  if(r%%2==0){ScrambledVec <- append(1,IDvec)}
  if(r%%2==1){ScrambledVec <- append(append(1,IDvec),r)}
  
  footprint <- ref[ScrambledVec,] %>%
               mutate(fp= paste0(lag(right),left,lag(R),L))%>% #fp = new neighbors
               select(fp)
  fpl <-nrow(footprint)
  fpVec <- as.vector(footprint[2:(fpl),])
  
  filler <- ref[!(1:nrow(ref)%in%IDvec),]
  filler <- filler[2:(nrow(filler)-1),] %>%
                mutate(fp = paste0(L,R,left,right))%>% #fp = original neighbors
                select(fp)
  flVec <- as.vector(unlist(filler))
  
  compVec <- append(fpVec,flVec)
  
  distMat <- adist(compVec,ignore.case=TRUE) #insertions,deltions, and substitutions are weighted equally
  distMat <- distMat[(fpl):(2*(fpl-1)),1:(fpl-1)] #rows=filler, col=borders  
  
  gapInd <- sort((ScrambledVec[2:(fpl)]-1))
  find <- function(x) {return(which(ScrambledVec==x))}
  #identify indices that each gap is forbidden from selecting 
  index <-  cbind(1:length(gapInd),as.vector(unlist(mapply(find,gapInd-1)))+1-1
                  ,as.vector(unlist(mapply(find,gapInd+1)))-1)
  index[index[,2]>length(IDvec),2] <- index[index[,2]>length(IDvec),3]
    
  #diag(distMat) <- 10e6 #particle may not return to original location
  distMat[index[,1:2]] <- 10e6 #prevent segment from rejoining with original neighbors
  distMat[index[,c(1,3)]] <- 10e6 #prevent segment from rejoining with original neighbors
  
  out <- solve_LSAP(distMat, maximum = FALSE)
  
  solution <- data.frame(cbind(gapInd,out))%>%arrange(out)
  
  order <- append(rep(NA,nrow(ref)-1),nrow(ref))
  order[1:nrow(ref)%in%ScrambledVec] <- ScrambledVec
  order[1:nrow(ref)%in%gapInd] <- solution[,1]
  
  return(order)
}