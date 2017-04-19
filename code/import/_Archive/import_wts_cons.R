#import/generate weight matrix
rm(list = ls())
setwd(getwd())
source('../utilities/header.R')
#####Read in Phylogeny Info#####
obbardTab3 <- read.csv(paste0(ReadPath,'phylogeny/obbard_tab3_exp.csv'),sep=',',header=FALSE)

names <- obbardTab3[,1]

tempChars <- obbardTab3[,3:8]

sub <- function(x){gsub('\\(.*\\)|\\s','',x)}
tempNums <- matrix(nrow=nrow(tempChars),ncol=ncol(tempChars))
for (i in 1:ncol(tempChars)){
  vec <- as.vector((mapply(sub,x=tempChars[,i])))  
  tempNums[,i] <- as.numeric(vec)
}
tempNums <- tempNums[,c(1,3,4,5,6)] #col2 is just mutation rate (not divegence time)
tempNorm <- scale(tempNums, center=FALSE, scale=colSums(tempNums))
#relative divergence times identical for all four models put forward by obbard et al.

obTime <- cbind(as.character(names),tempNums[,4],tempNorm[,4]) #Model A1 from obbard

######calcualte position relative to Dmel#####
#Names of comparison seq in order of appearance in alignment
wtNames <- c("Dp4","DroWil","DroAna","DroEre","DroMoj","DroSim","DroYak","DroVir")

#calculations for distance using obTime
#guided by Fig. 4 in obbard
#authors note this species problematic. However, it is clear from Fig. 4 in Obbard,
#as well as Fig.2 in "Seetharam and Stuart Whole genome phylogeny for 21 Drosophila species using predicted
#2b-RAD fragments" that branched from dmel at oldest point in shared tree
DroWil <- tempNorm[8,4] #Note: col. 4 == Model A1 
AnaCheck <- obTime[5,]
DroAna <- tempNorm[5,4]
EreCheck <- obTime[4,]
DroEre <- tempNorm[4,4]
MojCheck <- obTime[,8]
DroMoj <- tempNorm[8,4]
SimCheck <- obTime[2,]
DroSim <- tempNorm[2,4]
YakCheck <- obTime[4,]
DroYak <- DroEre #should be same as Erectus
VirCheck <- obTime[8,]
DroVir <- DroMoj #same as Moj

wtVec <- c(0,DroWil,DroAna,DroEre,DroMoj,
            DroSim,DroYak,DroVir)
seq_wts_cons <- data.frame(wtNames,wtVec)

write.table(seq_wts_cons, file = paste0(OutPath,"seq_wts_cons.csv")
            ,row.names=FALSE, na="",col.names=TRUE, sep=",")

