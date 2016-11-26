#import/generate weight matrix
rm(list = ls())

WritePath <- 'C:/Users/Nicholas/Documents/GitHub/synth_es2/data/'

seq_wts_cons <- as.vector(rep(1,9))

write.table(seq_wts_cons, file = paste0(WritePath,"/seq_wts_cons.csv")
            ,row.names=FALSE, na="",col.names=TRUE, sep=",")
