#import and FIMO results and generate annotated ES2 seq
rm(list = ls())
setwd(getwd())
library(Biostrings)
library(dplyr)
source('../utilities/header.R')

#####Project#######
project <- '/ToFIMO/Results_12.01.16/'
#####IMPORTS######
es2Min <- as.numeric(substr(read.csv(paste0(ReadPath,'/ES2MinSeq.csv'), header = TRUE, sep = ",")[1,5],4,10))
TFranges <- read.csv(paste0(AnalyzePath,'/RFtfbs.csv'), header = TRUE, sep = ",")
csvFIMO <- read.csv(paste0(AnalyzePath,project,'fimo0.1.csv'),sep=',',header=TRUE)
FullTFRaw <- sapply(read.csv(paste0(ReadPath,'/redfly_tfbs.csv'), header = TRUE, sep = ",")[1:27,c(1,5,6)], as.character)

#####DEFINE VARIABLES, CLEAN SETS#####
offset <- 49 #offset between start of full ES2 and minimal
TFminNorm <- TFranges
TFminNorm[,3:4] <- TFminNorm[,3:4] - es2Min + 1 + offset
MinMax <- max(TFminNorm[,4])
seqL <- 898 #length of full dmel ES2
origin <- es2Min - offset

TFClean <- data.frame(FullTFRaw) %>%
           mutate(start = as.numeric(substr(coordinates,4,10))-origin + 1,
                  stop = as.numeric(substr(coordinates,13,19))-origin + 1,
                  name = tolower(substr(name,1,unlist(gregexpr('_',name))-1)),
                  pValue = 0,
                  pattern=name) %>%
           select(pattern, start, stop, sequence, pValue)

#####confrim that min and full coordinates match#####
check <- full_join(TFClean,
                   select(TFminNorm,-X), c("start" = "X2","stop" = "X3")) %>%
         arrange(start,stop)
#word

frame <- data.frame(1:seqL)

FIMOLong <- csvFIMO %>% 
              mutate(dummy=TRUE) %>%
              left_join(frame %>% mutate(dummy=TRUE)) %>%
              filter(start<= X1.seqL, stop>= X1.seqL) %>%
              mutate(bp = X1.seqL, tf = as.character(pattern)) %>%
              filter(bp <= MinMax, bp >= offset) %>%
              group_by(tf, bp) %>%
              summarize( minP = min(pValue)) %>%
              filter(!(tf %in% c("slp","cad")))

Es2MinLong <- TFminNorm %>%
              mutate(dummy=TRUE) %>%
              left_join(frame %>% mutate(dummy=TRUE)) %>%
              filter(X2<= X1.seqL, X3>= X1.seqL) %>%
              mutate(bp = X1.seqL, 
                     tf = tolower(as.character(X1)),
                     minP = 0) %>%
              distinct(tf,bp,minP)

pVals <- exp(seq(log(.001), log(.1), length.out = 100))
fimoFull <- mutate(FIMOLong, Fbin = 1) %>%
            full_join(mutate(Es2MinLong, Mbin = 1), c("tf","bp")) %>%
            mutate(minP = ifelse(!is.na(minP.x),minP.x,minP.y)) %>%
            select(-minP.y, -minP.x)

scoreMat <- vector()
for (val in pVals){
  fimo <- mutate(fimoFull,
                 Fbin = (!is.na(Fbin))*(minP <= val),
                 Mbin = !is.na(Mbin),
                 tp = Mbin*Fbin,
                 fp = Fbin*(Fbin-Mbin)) %>%
           group_by(tf) %>%
           summarize_each(funs(sum)) %>%
           mutate(pVal = val,
                  ftp = tp/Mbin,
                  ffp = fp/(fp+tp))
  scoreMat <- rbind(scoreMat,fimo)  
}

fpLim <- .3 #for minimal enhancer region I'm reasonably confident that REDFLY's
            #tfbs coordinates are close to comprehensive. Therefore, to be conservative,
            #I will find the pValue for each TF that gives highest tp % while returning
            #fewer than 30% fp.  
tpLim <- .5

scoresSorted <- arrange(scoreMat,tf,pVal) %>%
                filter(ftp >= tpLim) %>%
                group_by(tf) %>%
                mutate(rank = row_number()) %>%
                filter(((ffp <= fpLim)+(rank==1))>=1) %>%                
                summarize(max = max(ftp)) %>%
                inner_join(filter(scoreMat, (ffp <= fpLim)||(rank==1)), 
                           c("tf"="tf", "max"="ftp" )) %>%
                group_by(tf) %>%
                summarize(pChoice = min(pVal),
                          ftp = max(max),
                          ffp = min(ffp)) 
                
#take average repressor pvalue for sloppy and average activator (exclude vfl)for caudal
names <- c("slp", "cad")
SlCa <- mutate(scoresSorted, cadFlag = (tf%in%c("bcd","hb"))) %>%
        group_by(cadFlag) %>%
        summarize(pChoice = mean(pChoice)) %>%
        mutate(tf = names[cadFlag+1]) %>%
        select(tf,pChoice)

tfPval <- as.data.frame(rbind(SlCa,select(scoresSorted,tf,pChoice)))

outFIMO <- left_join(csvFIMO, mutate(
                     tfPval,tf = factor(tf)), c("pattern"="tf")) %>%
           filter(pValue <= pChoice) %>%
           select(pattern,start,stop,sequence, pValue)

TFBS <- rbind(outFIMO,TFClean)

write.csv(TFBS,file=paste0(AnalyzePath,project,'TFBSfiltered.csv'))
write.csv(tfPval,file=paste0(AnalyzePath,project,'TFBSplimits.csv'))
