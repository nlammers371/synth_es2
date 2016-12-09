#install.packages("genalg")
library(genalg)

left <- TFranges[1,1]
right <- TFranges[1,2]

DistLim <- floor((right-left)/2)
charMat <- matrix(c("A","T","G","C"),nrow=2)

refScore <- es2Scores[,left:right]
refSeq <- es2
seqL <- right-left+1
pad <- 10
y <- rep(1,14)
x <- y
x[2] <- 0
eval_fun <- function(x){
  readMat <- matrix(x,ncol=2)+1
  
  SeqNew <-  DNAString(paste0(charMat[readMat],collapse=''))
  
  padNew <- DNAString(paste0(refSeq[(left-pad):(left-1)],
                             SeqNew,refSeq[(right+1):(right+pad)]))
  
  ScoreNew <- pwmScoreDriver_fun(MAT=MasterPWM, SEQ=padNew, DRIVER=TFsumm)
  ScoreNew <- ScoreNew[,(pad+1):(pad+1+right-left)]
  
  current_solution_points <- sum((refScore-ScoreNew)^2)
  
  distance <- stringDist(c(paste0(refSeq[left:right]),paste0(SeqNew)), method = "hamming",
                                  diag = FALSE, upper = FALSE)[1]

  if (distance < DistLim)
    return(Inf) else return(current_solution_points)
  
}

iter = 20
POPSIZE = 20
mCHANCE <- .1
GAmodel <- rbga.bin(size = 2*seqL, popSize = POPSIZE, iters = iter, mutationChance = mCHANCE, 
                    elitism = T, evalFunc = eval_fun,verbose=TRUE)

bestSolution<-GAmodel$population[which.min(GAmodel$evaluations),]

