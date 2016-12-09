#check FIMO results against REDFLY TFBS. Just check coverage for now
rm(list = ls())
setwd(getwd())
library(Biostrings)
source('../utilities/header.R')

project <- '/ToFIMO/Results_12.01.16/'

#####OLD PATHS#####
OldPath <- "../../data/_Archive/intermediate/"
OldData <-'C:/Users/Nicholas/Documents/GitHub/synth_es2/data/'
FunPath <- 'C:/Users/Nicholas/Documents/GitHub/synth_es2/code/_Archive/functions/'
#####OLD DATA#####
TFranges <- read.csv(paste0(OldPath,'/TFranges.csv'), header = TRUE, sep = ",")[,2:3]
es2Min <- (readDNAStringSet(file=paste0(OldPath,'dmel_es2.fa')))
es2MinChar <- paste0(unlist(es2Min))

#call tfbs annotation function
annotate_tfbs_fun <- dget(paste0(FunPath,"annotate_tfbs_fun.R")) 
es2AnnoChar <- gsub('F','-',annotate_tfbs_fun(seq=es2MinChar,tfbs=TFranges))
oldCoverage <- nchar(gsub('-','',es2AnnoChar))/nchar(es2AnnoChar)##.47 use for prior in phastcons

#####READ IN NEW DATA#####
es2CASEChar <- as.character(unlist(read.csv(paste0(
              AnalyzePath,project,'dmel_ES2Case.txt'),header=FALSE)[1]))

newCoverage <-nchar(gsub('[a,t,c,g]','',es2CASEChar))/nchar(es2CASEChar)##.50--good
  
es2Full <- DNAStringSet(es2CASEChar)

combEs2 <- c(es2Full,es2Min)
es2Algn <- AlignSeqs(combEs2)

##output and view
writeXStringSet(es2Algn,file=paste0(OutPath,'es2_check.fa'))

