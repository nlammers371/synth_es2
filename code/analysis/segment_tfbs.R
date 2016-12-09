#import and FIMO results and generate annotated ES2 seq
rm(list = ls())
setwd(getwd())
library(Biostrings)
library(dplyr)
source('../utilities/header.R')

project <- '/ToFIMO/Results_12.01.16/'


bpUnit <- 6 #minimum "meaningful" length scale
minStretch <- bpUnit #nothing shorter than this counts as an independent sequence
                       #most tfbs are < 10 bp
dfFIMO <- read.csv(paste0(AnalyzePath,project,'TFBSfiltered.csv'),sep=',',header=TRUE)
seqDNA <- unlist(readDNAStringSet(paste0(AnalyzePath,project,'dmelES2.fa'))) 

seqCharVec <- as.vector(strsplit(tolower(seqDNA),NULL)[[1]])

for(i in 1:nrow(dfFIMO)){
  seqCharVec[dfFIMO[i,3]:dfFIMO[i,4]] <- toupper(seqCharVec[dfFIMO[i,3]:dfFIMO[i,4]])
}

seqCharVec[(length(seqCharVec)-19):length(seqCharVec)] <- toupper(seqCharVec[(length(seqCharVec)-19):length(seqCharVec)])
seqCharVec[1:20] <- toupper(seqCharVec[1:20])
seqCASE <- paste0(seqCharVec, collapse='')


segment_seq_func <- dget(paste0(FunPath,"segment_seq_func.R"))
out <-segment_seq_func(seqDNA=seqDNA,seqCASE=seqCASE
                       ,bpUnit=bpUnit,minStretch=minStretch)

write(seqCASE, file = paste0(AnalyzePath,project,'dmel_ES2Case.txt'))
write.csv(out,file=paste0(AnalyzePath,project,'SegmentedES2_ms',minStretch,'TFBS.csv'))


