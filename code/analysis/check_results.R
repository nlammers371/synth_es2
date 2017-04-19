#Let's Try to get a feel for the results
library(RColorBrewer)
rm(list = ls())
setwd(getwd())
library(Biostrings)
library(ggplot2)
library(reshape)
source('../utilities/header.R')
source(paste0(FunPath,'utilities.R'))

#set project folder to draw from
project <- '/ToFIMO/Results_12.01.16/'

#import dmel ES2 with p<.005 binding sites in CAPS
es2CASEChar <- as.character(unlist(read.csv(paste0(
               AnalyzePath,project,'dmel_ES2Case.txt'),header=FALSE)[1]))

#import conservation scores & cons char string
RHO <- .3 #relative rate of cons substitution
conslim <- .8 #threshold for "conserved" classification 


cons_scores <- read.csv(file=paste0(AnalyzePath,project,'SegmentedES2_CONS','scores.csv'))[,2:3]
es2ConsChar <- as.character(unlist(read.csv(paste0(
                  AnalyzePath,project,'dmel_es2CASECons_rho',RHO,'cl',consLim,'.txt'),header=FALSE)[1]))

#plot tfbs footprint vs. cons_score 

tfbsVec <- rep(0,nchar(es2CASEChar))
tfbsInd <- as.vector(unlist(gregexpr('[A,T,G,C]',es2CASEChar)))
tfbsVec[tfbsInd] <- 1
ConsTfbsToPlot <- cbind(cons_scores,tfbsVec)
ConsTfbsToPlot <- melt(ConsTfbsToPlot,id='coord')

a <- ggplot(cons_scores,aes(coord,post.prob))+geom_area(color='dodgerblue3',fill='dodgerblue3',alpha=.5)
a <- ggplot(data=ConsTfbsToPlot)
a <- a + geom_line(aes(x=coord, y=value, color = variable))
a <- a + geom_ribbon(aes(ymin=0, ymax=value, x=coord,  fill = variable), alpha = 0.3)
a <- a + ggtitle('Conservation Scores and Predicted TFBS by x Position')
a
