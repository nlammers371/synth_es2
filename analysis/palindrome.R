  rm(list = ls())
  
  library(ggplot2)
  
  setwd('C:/Users/Nicholas/Documents/GitHub/synth_es2/data/')
  FunPath <- 'C:/Users/Nicholas/Documents/GitHub/synth_es2/analysis/functions/'
  WritePath <- 'C:/Users/Nicholas/Documents/GitHub/synth_es2/data/out/'
  es2csv <- read.csv('./in/redfly/ES2MinSeq.csv', header = TRUE, sep = ",")
  RFtfbs <- read.csv('./intermediate/RFtfbs.csv', header = TRUE, sep = ",")[,2:4]

  ####################################################################
  ##Search for Palindromes
  ####################################################################  

  #extract ES2 range from csv field
  rangeChar <- as.character(es2csv[1,5])
  rangeChar <- sub(':','..',rangeChar)
  rangeChar  <- strsplit(rangeChar,'..',fixed = TRUE)[[1]]
  frame <- vector()
  Fstart = as.numeric(rangeChar[2])
  frame[1] <- as.numeric(rangeChar[2])- Fstart + 1
  frame[2] <- as.numeric(rangeChar[3])- Fstart + 1
  
  TFranges <- RFtfbs[,2:3]-Fstart + 1
  #create es2 seq char
  es2 <- as.character(es2csv[1,6])
  
  #call tfbs annotation function
  annotate_tfbs_fun <- dget(paste0(FunPath,"annotate_tfbs_fun.R")) 
  es2anno <- annotate_tfbs_fun(seq=es2,tfbs=TFranges)
  
  #call palindrome_fun
  palindrome_fun <- dget(paste0(FunPath,"palindrome_fun.R"))
  palindromes = palindrome_fun(seq=es2anno,length=6)
  #note: I no palindromes >6 bp

  ####################################################################
  ##Assess Output
  ####################################################################

  #looks like there are two potentially worthwhile set of breakpoints. 
  #Lets visualize (crudely)
  StrVec <- c('A', 'T', 'G', 'C', 'F', 'P', 'Q')
  es2Pal <- es2anno
  char1 <- 'PPPPPP'
  char2 <- 'QQQQQQ'
  es2Pal <- gsub(as.character(palindromes[1,2]),char1,es2Pal)
  es2Pal <- gsub(as.character(palindromes[1,4]),char1,es2Pal)
  es2Pal <- gsub(as.character(palindromes[2,2]),char2,es2Pal)
  es2Pal <- gsub(as.character(palindromes[2,4]),char2,es2Pal)
  
  es2num <- vector()
  for (i in 1:nchar(es2Pal)){
    es2num[i] <- which(StrVec==substr(es2Pal,i,i))
  }
  break1 <- 1*(es2num==6)
  break2 <- 1*(es2num==7)
  tfbs   <- 1*(es2num==5)-.5
  sequence <- 1*(es2num < 5)-.5
  x <- 1:nchar(es2Pal)
  mat <- as.data.frame(cbind(x,break1,break2,tfbs,sequence))
  ggplot(mat, aes(x)) + 
  geom_line(aes(y = break1), colour = 'blue') + 
  geom_line(aes(y = break2)) +
  geom_point(aes(y = sequence), color = 'gray') + 
  geom_point(aes(y = tfbs), color = 'green') + 
  ylim(0,1)

  ##Both breakpoints appear to be in aceptable and potentially interesting places 
    
  ####################################################################
  ##Generate Sequence Variants
  ####################################################################

  #Variant 1:
  start1 <- as.numeric(palindromes[1,1])
  end1 <- as.numeric(palindromes[1,3]) + 6 - 1

  Seq1 <- substr(es2,start1,end1)
  RevSeq1bp <- rev(strsplit(Seq1,NULL))[[1]]
  RevSeq1 <- ""
  for (i in 1:nchar(Seq1)){
    RevSeq1 <- paste0(RevSeq1,RevSeq1bp[i])
  }
  es2pal1 <- gsub(Seq1,RevSeq1,es2)
  es2pal1out <- as.data.frame(cbind(start1,end1,es2pal1))

  write.table(es2pal1out, file = paste0(WritePath,"es2_pal_variant1.csv"),row.names=FALSE, na="",col.names=TRUE, sep=",")

  #Variant 2:
  start2 <- as.numeric(palindromes[2,1])
  end2 <- as.numeric(palindromes[2,3]) + 6 - 1
  
  Seq2 <- substr(es2,start2,end2)
  RevSeq2bp <- rev(strsplit(Seq2,NULL))[[1]]
  RevSeq2 <- ""
  for (i in 1:nchar(Seq2)){
    RevSeq2 <- paste0(RevSeq2,RevSeq2bp[i])
  }
  es2pal2 <- gsub(Seq2,RevSeq2,es2)
  es2pal2out <- as.data.frame(cbind(start2,end2,es2pal2))

  write.table(es2pal2out, file = paste0(WritePath,"es2_pal_variant2.csv"),row.names=FALSE, na="",col.names=TRUE, sep=",")
