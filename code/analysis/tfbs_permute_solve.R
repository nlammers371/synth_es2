#Implementation of Hungarian Algorithm to solve for best rearrangement or segments
rm(list = ls())
setwd(getwd())
library(Biostrings)
library(clue)
source('../utilities/header.R')
project <- 'Results_12.01.16'


files = c("c0.8ms6CONS","c0.8ms18CONS","c0.8ms36CONS",
          "p0.005ms6TFBS", "p0.005ms18", "p0.005ms24TFBS")

for (file in files){
  csvSEG <- read.csv(paste0(AnalyzePath,'/ToFIMO/',project,'/','SegmentedES2_',file,'.csv'),sep=',',header=TRUE)
  
  footprint <- data.frame(csvSEG)%>%
                    select(X,seq)
  
  distMat <- adist(footprint[2:(nrow(footprint)-1),2],ignore.case=TRUE) #insertions,deltions, and substitutions are weighted equally
  
  indVec <- 1:nrow(distMat)
  distMat[indVec%%2==0,indVec%%2==1] <- 10e6
  distMat[indVec%%2==1,indVec%%2==0] <- 10e6
  diag(distMat) <- 10e6
  
  
  out <- solve_LSAP(distMat, maximum = FALSE)
  
  solution <- cbind(1:nrow(footprint),append(append(1,as.vector(out)+1,),nrow(footprint)))
  
  solOrdered <- footprint[out,2]
  
  seqChar <- paste0(unlist(solOrdered),collapse='')
  seqDNA <- DNAStringSet(seqChar)
  seqOut <- csvSEG[out,][,1:4]
  
  write(seqChar, file = paste0(OutPath,'/ResultsTxt/',file,'_iter_','_',project,'SWAP.txt'))
  writeXStringSet(seqDNA,file=paste0(OutPath,'/ResultsFasta/',file,'_iter_','_',project,'SWAP.fa'))
  write.csv(seqOut,file=paste0(OutPath,'/ResultsCsv/',file,'_',project,'SWAP.csv'))
  
}