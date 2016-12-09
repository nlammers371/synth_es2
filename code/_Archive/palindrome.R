rm(list = ls())

library(ggplot2)

setwd('C:/Users/Nicholas/Documents/GitHub/synth_es2/data/')

FunPath <- 'C:/Users/Nicholas/Documents/GitHub/synth_es2/analysis/functions/'
WritePath <- 'C:/Users/Nicholas/Documents/GitHub/synth_es2/data/out/'

es2csv <- read.csv('./in/redfly/ES2MinSeq.csv', header = TRUE, sep = ",")
TFranges <- read.csv('./intermediate/TFranges.csv', header = TRUE, sep = ",")[,2:3]

####################################################################
##Search for Palindromes
####################################################################  

#create es2 seq char
es2Char <- as.character(es2csv[1,6])

#call tfbs annotation function (replaces all tfbs bp with 'F')
annotate_tfbs_fun <- dget(paste0(FunPath,"annotate_tfbs_fun.R")) 
es2AnnoChar <- annotate_tfbs_fun(seq=es2Char,tfbs=TFranges)

#call palindrome_fun
lpal = 6#define desired length of palindrome seq
palindrome_fun <- dget(paste0(FunPath,"palindrome_fun.R"))
palindromes = palindrome_fun(seq=es2AnnoChar,len=lpal)
#note: I find 1 6bp palindrome pair. No palindrome pairs >6 bp

####################################################################
##Assess Output
####################################################################

#looks like there are two potentially worthwhile sets of breakpoints. 
#Lets visualize (crudely)
StrVec <- c('A', 'T', 'G', 'C', 'F', 'P', 'Q')
es2PalChar <- es2AnnoChar
char1 <- 'PPPPPP'
char2 <- 'QQQQQQ'

#replace palindrome pairs with string vectors
es2PalChar <- gsub(as.character(palindromes[1,2]),char1,es2PalChar)
es2PalChar <- gsub(as.character(palindromes[1,4]),char1,es2PalChar)
es2PalChar <- gsub(as.character(palindromes[2,2]),char2,es2PalChar)
es2PalChar <- gsub(as.character(palindromes[2,4]),char2,es2PalChar)

es2num <- vector()
for (i in 1:nchar(es2PalChar)){
  es2num[i] <- which(StrVec==substr(es2PalChar,i,i))
}

break1 <- 1*(es2num==6)
break2 <- 1*(es2num==7)
tfbs   <- 1*(es2num==5)-.5
sequence <- 1*(es2num < 5)-.5
x <- 1:nchar(es2PalChar)
mat <- as.data.frame(cbind(x,break1,break2,tfbs,sequence))
ggplot(mat, aes(x)) + 
geom_step(aes(y = break1), colour = 'blue') + 
geom_step(aes(y = break2)) +
geom_point(aes(y = sequence), color = 'gray') + 
geom_point(aes(y = tfbs), color = 'green') + 
ylim(0,1)

##breakpoints appear to be in acceptable and potentially interesting places 
  
####################################################################
##Generate Sequence Variants
####################################################################

#Variant 1:
start1 <- as.numeric(palindromes[1,1])
end1 <- as.numeric(palindromes[1,3]) + lpal - 1

Seq1 <- substr(es2Char,start1,end1)
RevSeq1 <- paste0(rev(strsplit(Seq1,NULL))[[1]],collapse = '')
es2pal1 <- gsub(Seq1,RevSeq1,es2Char)
es2pal1out <- as.data.frame(cbind(start1,end1,es2pal1))

write.table(es2pal1out, file = paste0(WritePath,"/csv/es2_pal_variant1.csv"),row.names=FALSE, na="",col.names=TRUE, sep=",")
writeXStringSet(DNAStringSet(es2pal1), file=paste0(WritePath,"/fasta/es2_pal_variant1.fa"), width=nchar(es2pal1))
