#import+align minimal dmel es2 and comparison sequences
rm(list = ls())

library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(DECIPHER)

setwd('C:/Users/Nicholas/Documents/GitHub/synth_es2/data/')
InPath <- './in/'
WritePath <- './intermediate/'

DepEs2 <- readDNAStringSet("./intermediate/dmel_es2.fa") #minimal Dmel ES2 from Depace Paper
CompEs2 <- readDNAStringSet(paste0(InPath,'es2_comparison_seq.fa'))                                                  
CombEs2 <- c(DepEs2,LudEs2)

Es2Aligned <- AlignSeqs(CombEs2)
#index for non-missing sites in min dmel ES2
dmelInd <- as.vector(gregexpr('[A,T,G,C]',Es2Aligned[1])[[1]])
dmelLen <- width(DepEs2)

#create set containing only relevant stretches of each comparison sequence
Es2Min <- DNAStringSet(Es2Aligned, start = 1, end = dmelLen)
for (i in 1:length(Es2Min)){
  Fstring <- Es2Aligned[[i]]
  Tstring <- Fstring[dmelInd]
  Es2Min[[i]] <- Tstring 
}


writeXStringSet(Es2Min,file=paste0(WritePath,'es2_cons_align_lud.fa'))
writeXStringSet(Es2Aligned,file=paste0(WritePath,'es2_cons_full_align_lud.fa'))



