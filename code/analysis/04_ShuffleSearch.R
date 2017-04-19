#Variational Approach to minimize edge conflicts
rm(list = ls())
setwd("C:/Users/Nicholas/Documents/GitHub/synth_es2/code/analysis")
library(Biostrings)
library(clue)
source('../utilities/header.R')
source(paste0(FunPath,'utilities.R'))
project <- 'Results_12.01.16'
hungarian_solve_gaps_func <- dget(paste0(FunPath,"hungarian_solve_gaps_func.R")) 
hungarian_shuffle_search <- dget(paste0(FunPath,"hung_shuffle_search.R")) 

files = c("ms6TFBS", "ms36TFBS","c0.8ms6CONS","c0.8ms36CONS")
         # "p0.005ms6TFBS", "p0.005ms18TFBS", "p0.005ms24TFBS")

#Define vars
bpUnit = 6
ITER = 100
ptm <- proc.time()
for (file in files){
  csvSEG <- read.csv(paste0(AnalyzePath,'/ToFIMO/',project,'/','SegmentedES2_',file,'.csv'),sep=',',header=TRUE)
  
  
  results <- hungarian_shuffle_search(iter=ITER,segDF=csvSEG,bpUnit=6)
  
  minScore <- min(results[[3]])
  minID <- min(which(results[[3]]==minScore))
  BestSeq <- results[[2]][minID,]
  seqChar <- paste0(csvSEG[BestSeq,4],collapse='')
  seqDNA <- DNAStringSet(seqChar)
  
  seqOut <- csvSEG[BestSeq,][,1:4]
  
  write(seqChar, file = paste0(OutPath,'/ResultsTxt/',file,'_iter_',ITER,'_',project,'.txt'))
  writeXStringSet(seqDNA,file=paste0(OutPath,'/ResultsFasta/',file,'_iter_',ITER,'_',project,'.fa'))
  write.csv(seqOut,file=paste0(OutPath,'/ResultsCsv/',file,'_iter_',ITER,'_',project,'.csv'))
  write.csv(as.vector(results[[3]]),file=paste0(OutPath,'/ResultsCsv/',file,'_iter_',ITER,'_',project,'scores.csv'))
  
} 
proc.time() - ptm