function(iter,segDF,bpUnit){
  #Function to Search for novel segment confguration that minimizes edge conflicts
  
  #iter (int): number of independent configurations to test
  #segDF (DataFrame): frame containing segment boundary index and edge seq info
  #bpUnit: minimum threshold to qualify as "conserved" seq. Used to weight boundary dif scores
  
  #set diff weights. The closer a nucleotide is to the boundary of a 
  #segment the more independent potential motifs it will disrupt
  lwt <- bpUnit:1
  rwt <- 1:bpUnit
  Allwt <- c(lwt,rwt,rwt,lwt)
  csvSEG <- data.frame(lapply(csvSEG, as.character), stringsAsFactors=FALSE)
  
  csvSEG[is.na(csvSEG)] = paste0(rep('F',bpUnit),collapse='')
  
  idSEG <- csvSEG[,c(5:8)]
  for (i in 1:nrow(idSEG)){
    split <- paste0(csvSEG[i,5:8],collapse='')
    #amplify each nc based upon proximity to edge
    #ATGCGA -> AAAAAATTTTTGGGGCCCGGA
    str <- strrep(unlist(strsplit(split,NULL)),Allwt)
    idSEG[i,1] <- paste0(str[rwt],collapse = '')
    idSEG[i,2] <- paste0(str[(bpUnit+rwt)],collapse = '')
    idSEG[i,3] <- paste0(str[(2*bpUnit+rwt)],collapse = '')
    idSEG[i,4] <- paste0(str[(3*bpUnit+rwt)],collapse = '')
  }
  
  idVec <- seq(3,(nrow(idSEG)-1),2) 
  idVec <- idVec[idVec%%2==1] #randomly select islands, solve for optimal configuration within gaps

  #sampling matrix
  selectMat <- t(matrix(rep(1:length(idVec),length(idVec)),nrow=length(idVec)))
  selectMat <- selectMat - diag(diag(selectMat))
  filterMat <- selectMat[,1:(ncol(selectMat)-1)]
  for (i in 1:nrow(selectMat)){
    filterMat[i,] <- selectMat[i,selectMat[i,]!=0]
  }
  #random samples
  set.seed(123)
  r <- nrow(filterMat)
  row <- sample(1:r, 1)#randomly choose row to start
  shiftInd <- shifter(1:r,-(1+r-row))#shift index to keep track of order
  
  max <- 10*iter
  m <- 0
  control <- 0
  input <- vector()
  output <- vector()
  scores <- vector()
  
  while (control < iter && m <= max){
    sampVec <- as.vector(rep(NA,r))
    opt <- 1:r
    m <- m+1
    for (i in 1:r){
      shift <- shiftInd[i]#identify first segment
      sampleVec <- filterMat[shift,filterMat[shift,]%in%opt]
      if (length(sampleVec)>0){
        rand <- sample(sampleVec,1)#select index to insert new
        sampVec[rand] <- idVec[shift]    
        #remove selected val from opt
        opt <- opt[opt!=rand]
      }
      else{break}
    }
    if (max(is.na(sampVec))==0){
      vec <-hungarian_solve_gaps_func(IDvec=sampVec,ref=idSEG)
      output <- rbind(output,vec)
      input <- rbind(input,sampVec)
      scores <- append(scores,score(ref=idSEG,vec=vec,bpUnit=bpUnit,df=segDF))
      control<-control+1
    }
  }
  out <- list(input,output,scores)
  return(out)
}  
