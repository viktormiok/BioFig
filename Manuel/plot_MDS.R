plot_MDS<-function(expression,group,point.size,
                   LegendName_Color="group",
                   LegendName_Shape="shape",...){
 
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
    geom_point(size = point.size) +
    ggTheme(1) +
    labs(x = "Leading LogFC dim 1", y = "Leading LogFC dim 2", title = "MDS plot") +
    labs(shape=LegendName_Shape, col=LegendName_Color)+
    ggrepel::geom_text_repel(data = mds,aes(label = rownames(mds)))
  return(plotmds)
    
 
}
