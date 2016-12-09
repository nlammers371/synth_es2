#Let's Try to get a feel for the results

rm(list = ls())
setwd(getwd())
library(RColorBrewer)
library(Biostrings)
library(ggplot2)
library(reshape)
source('../utilities/header.R')
source(paste0(FunPath,'utilities.R'))

ITER=100
project <- 'Results_12.01.16'

type <- 'Conserved Sequences'
files = c("c0.8ms6CONS","c0.8ms36CONS","ms6TFBS", "ms36TFBS"
          )

for (file in files){
  seqResults <- read.csv(file=paste0(OutPath,'/ResultsCsv/',file,'_iter_',ITER,'_',project,'.csv'))
  mInd1 <- gregexpr('ms',file)[[1]]+2
  mInd2 <- gregexpr('\\d[A-Z]',file)
  tInd <- gregexpr('[A-Z]{4,}',file)[[1]]
  minStretch <- substr(file,mInd1,mInd2)
  type <- substr(file,tInd,tInd+3)
  plotArray <- vector()
  index <- 1
  for (i in 1:nrow(seqResults)){
    a <- rep(0,seqResults[i,4]-seqResults[i,3]+1)
    consInd <- gregexpr('[A-Z]',seqResults[i,5])[[1]]
    a[consInd] <- 1
    seq <- cbind(factor(a),seqResults[i,1],seqResults[i,3]:seqResults[i,4]
        ,index:(index+seqResults[i,4]-seqResults[i,3]))
    index <- index+ seqResults[i,4]-seqResults[i,3]+ 1
    plotArray <- rbind(plotArray,seq)
  }
  
  ToPlot <- data.frame(rbind(plotArray[,1:4],cbind(plotArray[,1:3],plotArray[,3])))
  
  ToPlot[,2] <- factor(as.character(ToPlot[,2]))
  ToPlot <- ToPlot[ToPlot[,4]%%2==1,]

  c <- ggplot(ToPlot, aes(x=X4, y=X3,color=X2, alpha = X1)) 
  c <- c + geom_point(aes())+ theme(legend.position="none")
  c <- c + labs(list(title = paste0('Original vs. Shuffled ES2:'),
                     subtitle= paste0('Minimum Permitted Stretch: '
                    ,minStretch, ' Type:',type), x = "bp (from)", y = "bp (to)"))
  
  c
  ggsave(paste0(OutPath,'/plots/',project,'_',file,'.pdf'))
}