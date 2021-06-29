plot_MDS<-function(expression,group,...){
 
  ## calculate distance for the sample
  data <- expression %>%
    t() %>%
    dist() %>%
    as.matrix()
  
  ## convert distance matrix to Classical multidimensional scaling(MDS)
  mdsData <- data.frame(cmdscale(data))
  mds <- cbind(mdsData, as.data.frame(data)) # combine with distance with mds
  
  ## plot in ggplot2
  plotmds <- ggplot(mds, aes(X1, X2,color=group)) +
    geom_point(size = 3) +
    theme_minimal() +
    theme(axis.title=element_text(size = 12,face="bold", colour = "black"),
          axis.text = element_text(size = 12),
          axis.ticks = element_line(colour='black'),
          plot.title = element_text(hjust = 0.5,size=12,face="bold"),
          legend.position = "bottom",
          legend.title = element_text(color = "Black", size = 12, face = "bold"),
          legend.text=element_text(color = "Black", size = 12, face = "bold"),
          panel.background = element_blank(),
          panel.grid.major =  element_line(colour = "grey90", size = 0.2),
          panel.grid.minor =  element_line(colour = "grey98", size = 0.5),
          panel.border = element_rect(color='black',fill=NA)) +
    labs(x = "Leading LogFC dim 1", y = "Leading LogFC dim 2", title = "MDS plot") +
    ggrepel::geom_text_repel(data = mds,aes(label = rownames(mds)))
  return(plotmds)
    
 
}
