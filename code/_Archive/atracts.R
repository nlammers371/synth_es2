#Code to generate "shuffled" ES2 Enhancers
rm(list = ls())
library(dplyr)
library(Biostrings)
library(iterpc)

setwd('C:/Users/Nicholas/Documents/GitHub/synth_es2/data/')

WritePath <- './out/'
FunPath <- 'C:/Users/Nicholas/Documents/GitHub/synth_es2/analysis/functions/'
annotate_tfbs_fun <- dget(paste0(FunPath,"annotate_tfbs_fun.R")) 

#####Import#####

es2csv <- read.csv('./in/redfly/ES2MinSeq.csv', header = TRUE, sep = ",")
TFranges <- read.csv('./intermediate/TFranges.csv', header = TRUE, sep = ",")[,2:3]


#####Initial Cleaning Steps#####

#create es2 seq char
es2Char <- as.character(es2csv[1,6])

#call tfbs annotation function

es2AnnoChar <- gsub('F','-',annotate_tfbs_fun(seq=es2Char,tfbs=TFranges))

es2Anno <- DNAString(es2AnnoChar)
es2 <- DNAString(es2Char)

#############################################################################
##Search for A-Tracts (Defined as any continuous A or T series >= 3bp)
#############################################################################
dict0=PDict(c("AAA","TTT"))
mm <- matchPDict(dict0, es2Anno)

Atracts <- rbind(as.data.frame(mm[[1]]),as.data.frame(mm[[2]]))[,1:2]

KeepRegions <- as.data.frame(rbind(as.matrix(Atracts),as.matrix(TFranges))) %>%
               arrange(start) 
              
#Create table of contiguous tfbs regions
endvec <- vector()
startvec <- vector()
start <- 0
end <- 0
e <- 1
s <- 1

for (i in 1:nrow(KeepRegions)){
  if (KeepRegions[i,1] > end){
    startvec[s] <- KeepRegions[i,1]
    endvec[e] <- KeepRegions[i,2]
    start <- KeepRegions[i,1]
    end <- KeepRegions[i,2]
    s <- s + 1
    e <- e + 1}
  
  if (KeepRegions[i,1] <= end && KeepRegions[i,2] > (end+1)){
    startvec[s] <- end + 1
    endvec[e] <- KeepRegions[i,2]
    start <- end
    end <- KeepRegions[i,2]
    s <- s + 1
    e <- e + 1      
  }
}
st <- as.vector(c(0,0))
keepSeg <- as.data.frame(rbind(st,cbind(startvec,endvec)))%>%
           mutate(newEnd = ifelse(is.na(lead(startvec)), 485,lead(startvec))-1,
                  newStart = endvec + 1,
                  gap = newEnd - newStart + 1) %>%
           filter(gap > 0) %>%
           select(newStart, newEnd, gap)

es2KeepChar <- gsub('F*F','-',annotate_tfbs_fun(seq=es2Char,tfbs=keepSeg[,1:2]))
es2KeepVec <- strsplit(es2KeepChar,'-')[[1]]
es2KeepVec <- es2KeepVec[es2KeepVec!='']
wVec <- keepSeg[,3]

toSiteOut <- c()
iter <- length(wVec) + length(es2KeepVec)

for (i in 1:iter){
  if (i%%2){
    toSiteOut[i] <- as.character(wVec[ceiling(i/2)])
  }
  else{
    toSiteOut[i] <- es2KeepVec[i/2]
  }
}

write(toSiteOut[],'./intermediate/a_Tracts_to_site_out.txt')
