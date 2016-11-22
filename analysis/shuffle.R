#Code to generate "shuffled" ES2 Enhancers
rm(list = ls())
library(dplyr)
library(Biostrings)
library(iterpc)

setwd('C:/Users/Nicholas/Documents/GitHub/synth_es2/data/')

WritePath <- './out/'
FunPath <- 'C:/Users/Nicholas/Documents/GitHub/synth_es2/analysis/functions/'

es2csv <- read.csv('./in/redfly/ES2MinSeq.csv', header = TRUE, sep = ",")
TFranges <- read.csv('./intermediate/TFranges.csv', header = TRUE, sep = ",")[,2:3]

#create es2 seq char
es2Char <- as.character(es2csv[1,6])

#call tfbs annotation function
annotate_tfbs_fun <- dget(paste0(FunPath,"annotate_tfbs_fun.R")) 
es2AnnoChar <- gsub('F','-',annotate_tfbs_fun(seq=es2Char,tfbs=TFranges))

es2Anno <- DNAString(es2AnnoChar)
es2 <- DNAString(es2Char)

#The goal of this exercise is to shuffle the positions of uncharacterized
#sequences that separate known TFBS while preserving the local structure
#of those sequences, as well as the global ordering of TFBS

#############################################################################
##Identify candidate sequence segments
#############################################################################

#Any uncharacterized sequence of length >= 12bp is eligible. This threshold
#derives from my decision to take 6bp the basic structural for es2 
bound = 12
base = 6

Fbound <- append(as.vector(gregexpr('[A,T,C,G]-',es2AnnoChar)[[1]]),nchar(es2Char))
Lbound <- append(1,as.vector(gregexpr('-[A,T,C,G]',es2AnnoChar)[[1]] + 1))
boundaries <- as.data.frame(cbind(Lbound,Fbound))
seqID <- mutate(boundaries, 
                dist = Fbound - Lbound + 1) %>%
         filter(dist >= bound) %>%
         select(Fbound,Lbound) %>%
         mutate(tail = Fbound - base + 1,
                head = Lbound + base - 1)

#I find 8 uncharacterized stretches eligible for shuffling using a 12 bp bound
seqMatch <- vector()
for (i in 1:nrow(seqID)){
    head <- substr(es2AnnoChar,(seqID[i,2]),(seqID[i,4]))
    tail <- substr(es2AnnoChar,(seqID[i,3]),(seqID[i,1]))
    fringe <- paste0(head,tail)
    weights <- cbind(6:1,1:6) #place more emphasis on bp near tfbs boundary
    seqVec<-strrep(strsplit(fringe,NULL)[[1]],weights)
    seqMatch[i] <- paste0(seqVec,collapse = '')
}

#find permutation with that minimizes overall edit distance
#is there a better approach than searching the full space?
SeqDist <- adist(seqMatch)
l = length(seqMatch)

I <- getall(iterpc(l, ordered=TRUE)) #matrix containing all possible permutations
#Find shuffle scheme that is least disruptive to tfbs borders 
#(as defined by str edit distance)

#exclude all permutations with self loops
for (i in 1:l){
  I = I[I[,i]!=i,]
}

scores <- vector()
for (i in 1:nrow(I)){
  path <- I[i,]
  ind <- 1:l
  scores[i] <- sum(diag(SeqDist[path,ind]))
} 
min <- min(scores)
degeneracy <- length(scores[scores==min]) #check for multiple optimums
ID <- which(scores==min)
OptPath <- I[ID,]

#Generate scrambled ES2 using optimal insertion scheme
insertions <- as.data.frame(matrix(nrow=8,ncol=5))
es2shuff <- DNAString("")
start <- 1 
seqID <- arrange(seqID,Lbound)#This should be redundant 
for (i in 1:l){
  left <- seqID[i,2] 
  orig <- es2[start:(left-1)]
  New <- OptPath[i]
  insert <- es2[seqID[New,2]:seqID[New,1]]
  es2shuff <- paste0(es2shuff,orig,insert)
  i2 <- nchar(es2shuff)
  i1 <- i2-nchar(insert)+1
  insertions[i,] <- cbind(i1,i2,seqID[New,2]
                          ,seqID[New,1],as.character(insert))
  start <- seqID[i,1]+1  
}

es2shuffChar <- paste0(es2shuff,es2[start:length(es2)])
es2shuff <- DNAString(es2shuffChar)

#checks
nchar(es2shuffChar)
nchar(es2Char)
alphabetFrequency(es2shuff)
alphabetFrequency(es2)

#export
write.table(es2shuffChar, file = paste0(WritePath,"/csv/es2_shuff_variant.csv"),row.names=FALSE, na="",col.names=TRUE, sep=",")
write.table(insertions, file = paste0(WritePath,"/csves2_shuff_variant_insertions.csv"),row.names=FALSE, na="",col.names=TRUE, sep=",")

writeXStringSet(DNAStringSet(es2shuff), file=paste0(WritePath,"/fasta/es2_shuff_variant.fa"), width=nchar(es2shuffChar))
