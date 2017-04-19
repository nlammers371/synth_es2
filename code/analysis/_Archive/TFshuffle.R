rm(list = ls())
setwd(getwd())
library(Biostrings)
library(dplyr)
library(clue)
source('../utilities/header.R')

project <- '/ToFIMO/Results_12.01.16/'
pFimo <- .01
pFilt <- .005 #minTol for analysis

bpUnit <- 6 #minimum "meaningful" length scale
minStretch <- 3*bpUnit #nothing shorter than this counts as an independent sequence
#most tfbs are < 10 bp
csvFIMO <- read.csv(paste0(AnalyzePath,project,'dmel_fimo.csv'),sep=',',header=TRUE)
seqDNA <- unlist(readDNAStringSet(paste0(AnalyzePath,project,'dmelES2.fa'))) 

dfFIMO <- as.data.frame(csvFIMO)%>%filter(pValue <= pFilt) ##this filter step is important
seqCharVec <- as.vector(strsplit(tolower(seqDNA),NULL)[[1]])

for(i in 1:nrow(dfFIMO)){
  seqCharVec[dfFIMO[i,3]:dfFIMO[i,4]] <- toupper(seqCharVec[dfFIMO[i,3]:dfFIMO[i,4]])
}

seqCASE <- paste0(seqCharVec, collapse='')

sqL <- max(dfFIMO$stop)
joinSeq <- data.frame(1:seqL)

dfFIMOLong <-dfFIMO %>% 
             mutate(dummy=TRUE, ID = row_number()) %>%
             left_join(joinSeq %>% mutate(dummy=TRUE)) %>%
             filter(start<= X1.seqL, stop>= X1.seqL) %>%
             mutate(seq = sequence.1, bp = X1.seqL) %>%
             select(-dummy, -X1.seqL, -sequence.1)


dfBcd <- filter(dfFIMOLong,X.pattern == 'bcd') %>%
         select(ID, X.pattern, bp, start, stop, strand, score, pValue, seq) 

otherTF <- filter(dfFIMOLong, X.pattern != 'bcd') %>%
           mutate(flag = 1) %>%
           select(bp,flag)

BcdTagged <- left_join(x=dfBcd, y=otherTF, by=c("bp")) %>%
             distinct(bp,ID,strand,score,seq,flag) 

dfBcdFilt <- group_by(BcdTagged, ID,strand,seq)%>%
             mutate(sp = is.na(flag))%>%
             summarize(ind = min(sp)) %>%
             filter(ind==1) %>%
             sapply(as.character)

strVec <- vector()
RstrVec <- vector()
for(i in 1:nrow(dfBcdFilt)){
  dna <- DNAString(paste0(toString(dfBcdFilt[i,3])))
  if (dfBcdFilt[i,2] =="-"){dna = reverseComplement(dna)}
  strVec[i] <- paste0(dna)
  #RstrVec[i] <- paste0(rev(dna))
}

distMat <- adist(strVec)
  
distMat[distMat==0] <- 10e6

out <- solve_LSAP(distMat, maximum = FALSE)

dfBcdOrdered <- data.frame(cbind(dfBcdFilt[out,],dfBcdFilt)) %>%
                mutate(newID = as.numeric(ID.1),
                       oldID = as.numeric(ID))%>%
                select(-strand.1,-seq.1,-ind.1,-ID, -ID.1) 

bcdOrig <- select(dfFIMO,start,stop,score) %>%
           mutate(ID = row_number())

newCoord <- left_join(dfBcdOrdered, bcdOrig, c("newID"= "ID"))%>%
            mutate(newScore = score,
                   newStart = start,
                   newStop = stop) %>%
            select(-score, -start,-stop, -ind, -strand) %>%
            left_join(bcdOrig, c("oldID" = "ID" )) %>%
            mutate(oldScore = score,
                   oldStart = start,
                   oldStop = stop) 
              