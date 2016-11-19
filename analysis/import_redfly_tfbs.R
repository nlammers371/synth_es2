#Import TFBS footprints taken from redfly database
rm(list = ls())
InPath <- 'C:/Users/Nicholas/Documents/GitHub/synth_es2/data/in/redfly/tfbs/'
OutPath <- 'C:/Users/Nicholas/Documents/GitHub/synth_es2/data/intermediate/'
files <- dir(InPath)
tfbsRanges <- data.frame(matrix(ncol = 3, nrow = length(files)))

for (f in 1:length(files)){
  IN <- read.csv(paste0(InPath,files[f]), header = TRUE, sep = ",")
  
  tfbsRanges[f,1] <- gsub('_.*$',"",as.character(IN[1,1]))
  
  rangeChar <- as.character(IN[1,5])
  rangeChar <- sub(':','..',rangeChar)
  rangeChar  <- strsplit(rangeChar,'..',fixed = TRUE)[[1]]
  tfbsRanges[f,2] <- as.numeric(rangeChar[2])
  tfbsRanges[f,3] <- as.numeric(rangeChar[3])
}
write.csv(tfbsRanges, file = paste0(OutPath,'RFtfbs.csv'))
