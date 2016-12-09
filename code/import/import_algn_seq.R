#import+align minimal dmel es2 and comparison sequences
rm(list = ls())
setwd(getwd())
source('../utilities/header.R')
library(Biostrings)


#####Read in ES2 Sequences Compiled by Ciera####
CombES2 <- readDNAStringSet(paste0(SeqPath,'forContructTarget_eve-striped-2_with_Montium_and_melanogaster.fa'))                                                  
dmelES2 <- CombES2[18]
names(CombES2[18])
writeXStringSet(dmelES2,file=paste0(AnalyzePath,'dmelES2.fa'))


