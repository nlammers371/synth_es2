#Piece together a Newick Diagram for the species I from Ciera's Compilation
#install.packages("ape")
rm(list = ls())
setwd(getwd())
library(Biostrings)
library(ape)
source('../utilities/header.R')


Ref <- "(((MEMB003D,(MEMB003F,(MEMB002D,MEMB003B))),(MEMB002F,(MEMB002A,(MEMB002C,MEMB003C)))),((DroWil1,(DroMoj3,DroVir3)),(Dp4,(DroAna3,(DroEre2,(dmel,DroSim1))))));"
RefTree <- read.tree(text = Ref)
plot(RefTree, type = "phylogram")

write(Ref, file = paste0(AnalyzePath,'RefTree.txt'))
