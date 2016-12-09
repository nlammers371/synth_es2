#Variational Approach to minimize edge conflicts
rm(list = ls())
setwd(getwd())
library(Biostrings)
library(clue)
source('../utilities/header.R')
project <- '/ToFIMO/Results_11.30.16/'

csvSEG <- read.csv(paste0(AnalyzePath,project,'SementedES2_p0.01ms12.csv'),sep=',',header=TRUE)
bpUnit <- 6

lwt <- bpUnit:1
rwt <- 1:bpUnit
Allwt <- c(lwt,rwt,lwt,rwt)
csvSEG <- data.frame(lapply(csvSEG, as.character), stringsAsFactors=FALSE)

csvSEG[is.na(csvSEG)] = paste0(rep('F',bpUnit),collapse='')

DiffChar <- matrix(unlist(strsplit(paste0(unlist(csvSEG[,5:8]),collapse=''),NULL)),ncol=4*bpUnit)
idSEG <- csvSEG[,c(1,5:8)]
for (i in 1:nrow(idSEG)){
  str <- strrep(unlist(strsplit(paste0(DiffChar[i,],collapse=''),NULL)),Allwt)
  idSEG[i,2] <- paste0(str[rwt],collapse = '')
  idSEG[i,3] <- paste0(str[(bpUnit+rwt)],collapse = '')
  idSEG[i,4] <- paste0(str[(2*bpUnit+rwt)],collapse = '')
  idSEG[i,5] <- paste0(str[(3*bpUnit+rwt)],collapse = '')
}

idVec <- as.numeric(idSeg[2:(nrow(idSEG)-1),1])
shiftVec <- idVec
shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}
shiftVec[shiftVec%%2==0] <- shifter(shiftVec[shiftVec%%2==0],floor(length(shiftVec)/8))

idVec <- cbind(idVec,shiftVec)

score <- function(ref,vec){
  L <- length(vec)
  svec <- as.vector(append(1,vec,L))
  ldist <- diag(adist(c(ref[vec,4],ref[svec[1:(L)],2]))[1:L,(L+1):(2*L)])
  rdist <- diag(adist(c(ref[vec,5],ref[svec[1:(L)],3]))[1:L,(L+1):(2*L)])
  return(sum(ldist+rdist))
}
