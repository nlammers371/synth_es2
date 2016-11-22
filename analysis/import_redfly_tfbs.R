#Import TFBS footprints taken from redfly database
rm(list = ls())
library(dplyr)
setwd('C:/Users/Nicholas/Documents/GitHub/synth_es2/data')
InPath <- './in/redfly/tfbs/'
OutPath <- './intermediate/'
files <- dir(InPath)

tfbsRanges <- data.frame(matrix(ncol = 3, nrow = length(files)))
es2csv <- read.csv('./in/redfly/ES2MinSeq.csv', header = TRUE, sep = ",")
es2Char <- as.character(es2csv[1,6])
for (f in 1:length(files)){
  IN <- read.csv(paste0(InPath,files[f]), header = TRUE, sep = ",")
  
  tfbsRanges[f,1] <- gsub('_.*$',"",as.character(IN[1,1]))
  
  rangeChar <- as.character(IN[1,5])
  rangeChar <- sub(':','..',rangeChar)
  rangeChar  <- strsplit(rangeChar,'..',fixed = TRUE)[[1]]
  tfbsRanges[f,2] <- as.numeric(rangeChar[2])
  tfbsRanges[f,3] <- as.numeric(rangeChar[3])
}

#extract start point of ES2 
rangeChar <- as.character(es2csv[1,5])
rangeChar <- sub(':','..',rangeChar)
rangeChar  <- strsplit(rangeChar,'..',fixed = TRUE)[[1]]
frame <- vector()
Fstart = as.numeric(rangeChar[2])

TFranges <-tfbsRanges[,2:3]-Fstart + 1

freqCheck <- summarize(group_by(tfbsRanges,X1),count = n()) #compare to Depace footprint counts
ordered <- arrange(tfbsRanges,X2)


write.csv(tfbsRanges, file = paste0(OutPath,'RFtfbs.csv'))
write.csv(TFranges, file = paste0(OutPath,'TFranges.csv'))
writeXStringSet(DNAStringSet(es2Char), file=paste0(OutPath,"/dmel_es2.fa"), width=nchar(es2Char))
