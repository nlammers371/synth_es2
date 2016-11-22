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

#############################################################################
##Initial Cleaning Steps
#############################################################################

#create es2 seq char
es2Char <- as.character(es2csv[1,6])

#call tfbs annotation function
annotate_tfbs_fun <- dget(paste0(FunPath,"annotate_tfbs_fun.R")) 
es2AnnoChar <- gsub('F','-',annotate_tfbs_fun(seq=es2Char,tfbs=TFranges))

es2Anno <- DNAString(es2AnnoChar)
es2 <- DNAString(es2Char)

#############################################################################
##Search for A-Tracts (Defined as any continuous A or T series >= 3bp)
#############################################################################
dict0=PDict(c("AAA","TTT"))
mm <- matchPDict(dict0, es2Anno)
st <- as.vector(c(0,0))
end <- as.vector(c(485,0))
Atracts <- rbind(as.data.frame(mm[[1]]),as.data.frame(mm[[2]]))[,1:2]
KeepRegions <- as.data.frame(rbind(as.matrix(Atracts),as.matrix(TFranges),st,end)) %>%
               arrange(start) %>%
               mutate(interStart = end + 1,
                      interEnd = lead(start-1),
                      width = interEnd - interStart +1) %>%
               filter(interStart < interEnd) %>%
               select(interStart,interEnd, width)
              

es2KeepChar <- gsub('F*F','-',annotate_tfbs_fun(seq=es2Char,tfbs=KeepRegions[,1:2]))
wVec <- KeepRegions[,3]
es2KeepVec <- strsplit(es2KeepChar,'-')[[1]]
es2KeepVec <- es2KeepVec[es2KeepVec != ""]
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
