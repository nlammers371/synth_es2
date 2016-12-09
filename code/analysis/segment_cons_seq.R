#import and FIMO results and generate annotated ES2 seq
rm(list = ls())
setwd(getwd())
library(Biostrings)
library(dplyr)
source('../utilities/header.R')

#parameters used to generate ref files
consLim <- .8
RHO <- .3
project <- '/ToFIMO/Results_12.01.16/'

seqCASEraw <- as.character(unlist(read.csv(paste0(AnalyzePath,project,
                                  'dmel_es2CASECons_rho',RHO,'cl',consLim,'.txt'),header=FALSE)))

seqCharVec <- unlist(strsplit(seqCASEraw,NULL))
seqCharVec[(length(seqCharVec)-19):length(seqCharVec)] <- toupper(seqCharVec[(length(seqCharVec)-19):length(seqCharVec)])
seqCharVec[1:20] <- toupper(seqCharVec[1:20])
seqCASE <- paste0(seqCharVec, collapse='')

dmeles2 <- DNAString(seqCASE)
#####CALL SEGEMENTATION FUNC#####
segment_seq_func <- dget(paste0(FunPath,"segment_seq_func.R"))
bpUnit <- 6
minStretch <- bpUnit #not very sensitive to param val
out <-segment_seq_func(seqDNA=dmeles2,seqCASE=seqCASE
                       ,bpUnit=bpUnit,minStretch=minStretch)


write.csv(out,file=paste0(AnalyzePath,project,'SegmentedES2_c',consLim,'ms',minStretch,'CONS.csv'))
