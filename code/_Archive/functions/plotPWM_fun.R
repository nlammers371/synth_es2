function(scores){
  #simple function to plot PWM scores
  scores[is.na(scores)] <- 0  
  PlotData <- as.data.frame(melt(scores,id="TFnames")) %>%
    select(TFnames,value) %>%
    arrange(TFnames)
  PlotData[is.na(PlotData[,2]),2]<-0
  yMin <- min(PlotData[PlotData[,2]!=0,2])
  yMin <- yMin - .1*yMin
  yMax <- max(PlotData[,2])
  yMax <- yMax + .1*yMax
  
  xVec <- rep(1:(length(scores[1,])-1),length(scores[,1]))
  a <- ggplot(data=PlotData)
  a <- a + geom_line(aes(x=xVec, y=value, color = TFnames))
  a <- a + geom_ribbon(aes(ymin=0, ymax=value, x=xVec,  fill = TFnames), alpha = 0.3)
  a <- a + coord_cartesian(ylim = c(yMin,yMax)) + ggtitle('PWM Score by x Position')
  a <- a + ylab("Normalized 'Energy'") + xlab('Position Relative to Start of ES2')
  
  return(a)
}