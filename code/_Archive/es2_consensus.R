#Look across es2 sequences from multiple species to find most highly conserved pieces
rm(list = ls())
setwd('C:/Users/Nicholas/Documents/GitHub/synth_es2/data/')

InPath <- './intermediate/'
WritePath <- './out/'
FunPath <- 'C:/Users/Nicholas/Documents/GitHub/synth_es2/analysis/functions/'

#####Import#####
annotate_tfbs_fun <- dget(paste0(FunPath,"annotate_tfbs_fun.R")) 
cons_score_fun <- dget(paste0(FunPath,"cons_score_fun.R")) 

es2csv <- read.csv('./in/redfly/ES2MinSeq.csv', header = TRUE, sep = ",")
TFranges <- read.csv('./intermediate/TFranges.csv', header = TRUE, sep = ",")[,2:3]
wtVec <- as.numeric(read.csv(paste0(InPath,'seq_wts_cons.csv'),
                                         header = TRUE, sep = ",")[,1])
Alignment <-  readDNAStringSet("./intermediate/es2_cons_align_lud.fa")

#####Cleaning#####
#create es2 seq char
es2Char <- as.character(es2csv[1,6])

#call tfbs annotation function
es2AnnoChar <- gsub('F','-',annotate_tfbs_fun(seq=es2Char,tfbs=TFranges))
es2Anno <- DNAString(es2AnnoChar)
es2 <- DNAString(es2Char)

cons_seq <- cons_score_fun(alignment=Alignment, wtVec=wtVec, seqAnno=es2AnnoChar
                           , seq=es2Char, lb=.5, scores=0)
