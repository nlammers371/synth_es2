#Scriot to align sequences and obtain conservation scores

rm(list = ls())
setwd('C:/Users/Nicholas/Documents/GitHub/synth_es2/code/analysis')
source('../utilities/header.R')
library(rphast)
library(Biostrings)
library(DECIPHER)
library(dplyr)
library(ape)
library(ggplot2)
project <- '/ToFIMO/Results_12.01.16/'

#Read in ES2 Sequences Compiled by Ciera
CombES2 <- readDNAStringSet(paste0(SeqPath,'forContructTarget_eve-striped-2_with_Montium_and_melanogaster.fa'))
RefTree <- as.character(unlist(read.table(paste0(AnalyzePath,'RefTree.txt'))))

#Read in output of TFBS segment for use as reference
es2CASEChar <- as.character(unlist(read.csv(paste0(
                    AnalyzePath,project,'dmel_ES2Case.txt'),header=FALSE)[1]))

#-------------------------Format and Clean Data-----------------------------------#
#Plot RefTree for visual Check 
RefTreePlot <- read.tree(text = RefTree)
plot(RefTreePlot, type = "phylogram")

names <- unlist(strsplit(names(CombES2),'\\|'))
ind  <- seq(3,length(names),4)
ClnNames <- names[ind] #Note: MEMB002E does not appear in phylo data. remove
ClnNames <- c(ClnNames[1:8],ClnNames[10:18])
CombES2 <- c(CombES2[1:8],CombES2[10:18])
#NOTE: dmel should be on top to ensure compatible coordinates for cons scores
ClnNames <- c(ClnNames[17],ClnNames[1:16])
CombES2 <- c(CombES2[17],CombES2[1:16])
ClnNames[1] 
CombES2[1]

#---------------------Align Sequences and Estimate Tree----------------------------#
ES2Aligned <- AlignSeqs(CombES2)

msaES2 <- msa(paste(ES2Aligned),names=ClnNames)

#estimate treeeee
neutralMod <- phyloFit(msaES2,tree=RefTree)

#Let's take a look at the fitted tree
FitTree <- "(((MEMB003D:0.041213,(MEMB003F:0.0657486,(MEMB002D:0.0487667,MEMB003B:0.116597):0.0109954):0.013095):0.0565007,(MEMB002F:0.103063,(MEMB002A:0.00887196,(MEMB002C:0.00573016,MEMB003C:0.006826):0.00556543):0.121106):0.00923898):0.0553344,((DroWil1:0.355821,(DroMoj3:0.175613,DroVir3:0.13643):0.317412):0.0964739,(Dp4:0.295486,(DroAna3:0.21584,(DroEre2:0.0592241,(dmel:0.0187653,DroSim1:0.013545):0.035115):0.163677):0.0371299):0.0152591):0.0553344);"
FitTreePlot <- read.tree(text = FitTree)
plot(FitTreePlot, type = "phylogram")

#estimate expectedlength by taking average len of CAPS in CASE. This will give a
#higher avg length than FP coordinates because many sites overlap
splitMat <- unlist(strsplit(gsub('[a,t,g,c]','-',es2CASEChar),'-'))
splitMat <- splitMat[splitMat != ""]
Tlen <- nchar(paste0(splitMat,collapse=''))/length(splitMat)

coverage <- .47 #from REDFLY TFBS footprints
#sites for minimal es2 (~.47)
RHO <- .3 #I'm not sure about this value (relative evo rate). Check sensitivity
#results are quite robust to changes: .5--.32 cons,.3--.302,.1--.26
phast <- phastCons(msaES2, neutralMod,rho=RHO, target.coverage=coverage, expected.length=Tlen)

cons_scores <- phast$post.prob.wig
plot(cons_scores) #wow. Nice distinction
hist(cons_scores$post.prob) #Will set threshold at .8
consLim <- .8
cons_scoresFilt <- filter(cons_scores, post.prob >= consLim)
consInd <- cons_scoresFilt[,1]
es2Lower <- tolower(es2CASEChar)
es2split <- unlist(strsplit(es2Lower,NULL))
es2CASECons <- es2split
es2CASECons[consInd] <- toupper(es2CASECons[consInd])
es2CASECons <- paste0(es2CASECons,collapse='')

#Export DNA String with conserved regions in CAPS
write(es2CASECons, file = paste0(AnalyzePath,project,'dmel_es2CASECons_rho',RHO,'cl',consLim,'.txt'))
#Export Cons Scores
write.csv(cons_scores,file=paste0(AnalyzePath,project,'SegmentedES2_CONS','scores.csv'))
#Export Alignment
writeXStringSet(ES2Aligned,file=paste0(OutPath,'all_seq_algn.fa'))
