#given a matrix of PWMs, calculate normalized joint score
rm(list = ls())
library(dplyr)
library(Biostrings)
library(ggplot2)
library(reshape)
#library(reshape)
setwd('C:/Users/Nicholas/Documents/GitHub/synth_es2/data/')
InPath <- 'intermediate/'
WritePath <- './out/'
FunPath <- 'C:/Users/Nicholas/Documents/GitHub/synth_es2/analysis/functions/'

annotate_tfbs_fun <- dget(paste0(FunPath,"annotate_tfbs_fun.R"))
pwmScoreDriver_fun <- dget(paste0(FunPath,"pwmScoreDriver_fun.R"))
pwmScore_fun <- dget(paste0(FunPath,"pwmScore_fun.R")) #called inside pwmScoreDriver_fun
sub_seq_fun <- dget(paste0(FunPath,"sub_seq_fun.R"))
plotPWM_fun <- dget(paste0(FunPath,"plotPWM_fun.R"))

#####Read in PWM Matrix#####
MasterPWM <- read.table(file = paste0(InPath,'MasterPWM.csv'), header = TRUE, sep = ',')

TFsumm <- as.data.frame(group_by(MasterPWM,TF,ID)%>%
          summarize())
TFnames <- unique(TFsumm[,1])

#####Read in ES2 Sequence and TFBS ranges######
es2csv <- read.csv('./in/redfly/ES2MinSeq.csv', header = TRUE, sep = ",")
TFranges <- read.csv('./intermediate/TFranges.csv', header = TRUE, sep = ",")[,2:3]

#create es2 seq char
es2Char <- as.character(es2csv[1,6])

#call tfbs annotation function
es2AnnoChar <- gsub('F','-',annotate_tfbs_fun(seq=es2Char,tfbs=TFranges))

es2Anno <- DNAString(es2AnnoChar)
es2 <- DNAString(es2Char)

#####Read in Insertion Matrix######
#I want to only consider those regions that I selected for shuffling in the 
#shuffle perturbation

shuffRanges <- read.csv('./out/csv/es2_shuff_variant_insertions.csv',header=TRUE, sep=',')[1:7,1:2]
eSeq <- vector()
for (i in 1:nrow(shuffRanges)){
  eSeq[i] <- as.character(es2[shuffRanges[i,1]:shuffRanges[i,2]])  
}
eSeqDNA <- DNAStringSet(eSeq)

#####Analysis Start#####
#Calculate weight and es2 ref matrices that will remain constant throughout iteration
es2ScoresC <- pwmScoreDriver_fun(MAT=MasterPWM,SEQ=es2,DRIVER=TFsumm,COMP=1)
es2Scores <- pwmScoreDriver_fun(MAT=MasterPWM,SEQ=es2,DRIVER=TFsumm,COMP=0)
#####plot 5 known TFs as a sanity check######
PlotData <- cbind(TFnames,es2Scores)
ToPlot <- as.data.frame(PlotData[PlotData[,1]%in%c("bcd"),])
a <- plotPWM_fun(ToPlot)
PlotData1 <- cbind(TFnames,es2ScoresC)
ToPlot1 <- as.data.frame(PlotData1[PlotData1[,1]%in%c("bcd"),])
b <- plotPWM_fun(ToPlot1)

#####Generate Random Libraries for each Sequence#####

charVec <- c("A","G","C","T")
lib = 100

library <- matrix(sample(charVec,length(eSeqDNA[[1]])*lib,replace=TRUE),nrow=lib)





compSet <- DNAStringSet(c(es2,es2Anno))
sampVec <- as.vector(gregexpr('[A,T,G,C]',es2Anno)[[1]])
sampSeed <- strsplit(es2Char,NULL)[[1]]

sampMat<- matrix(ncol=4,nrow=length(sampSeed))
#make sure es2 bp is always in same position in samp mat
for (i in 1:length(sampSeed)){
  sampMat[i,] <- append(sampSeed[i],charVec[charVec!=sampSeed[i]])
}
es2Scrambled <- es2
es2Scrambled[sampVec] <- strsplit(sample(charVec,length(sampVec),replace=TRUE),NULL)[[1]]

inc <- 100 #routines
sub <- 10 #subroutines
NSUB <- 5 #substitutions per subroutine
TOL <-.5

delta <- vector()
distance <- vector()
accept <- vector()

es2Old <- es2Scrambled
ScoreOld <- pwmScoreDriver_fun(SEQ=es2Old,DRIVER=TFsumm)
deltaOld <- sum(wtsPMW*(ScoreOld-es2Scores)^2)
delta[1] <- deltaOld
distance[1] <-1-(stringDist(c(paste0(es2),paste0(es2Old)), method = "hamming",
                     diag = FALSE, upper = FALSE)[1])/length(sampVec)
iter <- 2

for (i in 1:inc){
  acc <- 0
  for(j in 1:sub){
    es2New <- sub_seq_fun(seq=es2Old,ref=es2,subVec=sampVec,nsub=NSUB,tol=TOL)
    ScoreNew <- pwmScoreDriver_fun(SEQ=es2New,DRIVER=TFsumm)
    deltaNew <- sum(wtsPMW*(ScoreNew-es2Scores)^2)
    
    if (deltaNew < deltaOld){
      
      ScoreOld <- ScoreNew
      deltaOld <- deltaNew
      es2Old <- es2New
      
      acc <- acc + 1
          
    }    
    delta[iter] <- deltaOld
    distance[iter] <- 1-(stringDist(c(paste0(es2),paste0(es2Old)), method = "hamming",
                             diag = FALSE, upper = FALSE)[1])/length(sampVec)
    iter <- iter + 1  
  }  
  accept[i] <- acc/sub  
}




