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

#Import Eve2 String with Conserved Regions in Caps
seqCASEraw <- as.character(unlist(read.csv(paste0(AnalyzePath,project,
                                  'dmel_es2CASECons_rho',RHO,'cl',consLim,'.txt'),header=FALSE)))


#-----------------------------Initial Cleaning------------------------------------#
#Split into letter vector
seqCharVec <- unlist(strsplit(seqCASEraw,NULL))
#Ensure that leading and trailing 20 characters are conserved
seqCharVec[(length(seqCharVec)-19):length(seqCharVec)] <- toupper(seqCharVec[(length(seqCharVec)-19):length(seqCharVec)])
seqCharVec[1:20] <- toupper(seqCharVec[1:20])
#Recombine 
seqCASE <- paste0(seqCharVec, collapse='')
#Convert to DNA String Object
dmeles2 <- DNAString(seqCASE)

#------------------------------CALL SEGEMENTATION FUNC-----------------------------#
segment_seq_func <- dget(paste0(FunPath,"segment_seq_func.R"))
bpUnit <- 6
#nothing shorter than this counts as an independent sequence
#output for cons sequences is not very sensitive to minstretch val
minStretch <- bpUnit 

#call function that breaks sequence into sub sections while preserving conserved regions
out <-segment_seq_func(seqDNA=dmeles2,seqCASE=seqCASE
                       ,bpUnit=bpUnit,minStretch=minStretch)


#Export
write.csv(out,file=paste0(AnalyzePath,project,'SegmentedES2_c',consLim,'ms',minStretch,'CONS.csv'))
