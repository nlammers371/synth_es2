#Import and format Pwm Matrices
rm(list = ls())
library(dplyr)
setwd('C:/Users/Nicholas/Documents/GitHub/synth_es2/data/in/pwm')
InPath <- './site_out/' #using siteout matrices for now
WritePath <- '../../intermediate/'
#import site_out namekey
siteOutKey <- read.table(file = paste0(InPath,'name_key.csv'), header = TRUE, sep = ',')

bkprob <- .25 #is this appropriate for ES2?

fnames = list.files(path = InPath, pattern="*.fm")
myfiles = lapply(paste0(InPath,fnames), read.table)

#define char lists for classifying tables by TF

getList <- function(pat){
  ind <- grepl(pattern = pat, x = fnames)
}
tfID = as.data.frame(cbind(t(mapply(getList,siteOutKey[,2]))))
iter <- 1
MasterPWM <- data.frame()
for (i in 1:nrow(tfID)){
  list <- myfiles[tfID[i,]==1]
  name <- siteOutKey[i,1]
  for (j in 1:length(list)){
    df <- as.data.frame(list[[j]])%>%
          mutate(TF = as.character(name),
                 ID = j, src = iter,
                 A = V1, C = V2, G = V3, T = V4) %>%
          select(TF,ID,A,C,G,T)
    df[3:6] <- log2(df[3:6]/bkprob)
    MasterPWM <- rbind(MasterPWM,df)
    iter <- iter + 1
  }
}

write.table(fnames, file = paste0(WritePath,"pwm_src.csv"),row.names=FALSE, na="",col.names=TRUE, sep=",")
write.table(MasterPWM, file = paste0(WritePath,"MasterPWM.csv"),row.names=FALSE, na="",col.names=TRUE, sep=",")
