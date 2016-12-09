#check for HindIII sites on synthetic ES2 
rm(list = ls())
setwd(paste0("C:/Users/Nicholas/Documents/GitHub/synth_es2/code/analysis"))
library(Biostrings)
library(dplyr)
source('../utilities/header.R')

files = list.files(path = paste0(OutPath,"ResultsFasta/"), pattern = "*.fa")
flags <- vector()
iter <- 1

hind1 <- DNAString("aagctt")
hind2 <- DNAString("ttcgaa")

for (file in files){
  assign(file,unlist(readDNAStringSet(file=paste0(OutPath,"ResultsFasta/",file))))
  flags[iter] <- ((max(gregexpr(hind1,eval(as.name(file)))[[1]]) > 0) +
                  (max(gregexpr(rev(hind1),eval(as.name(file)))[[1]]) > 0) +
                  (max(gregexpr(hind2,eval(as.name(file)))[[1]]) > 0) +
                  (max(gregexpr(rev(hind2),eval(as.name(file)))[[1]]) > 0))
                            
  iter <- iter + 1
}

#all clear

