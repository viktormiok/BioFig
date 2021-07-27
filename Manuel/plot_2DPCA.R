plot_2DPCA<-function(expression, scl,colors=NULL,shape=NULL,
                     group,samplenames,
                     ggrepelLab=TRUE,size_gglab=5,
                     size_title=1,point.size=4,
                     MahalanobisEllips=F,
                     LegendName_Color="group",LegendName_Shape="shape",LegendName="group",...){
  df_pca<-prcomp(t(expression),scale=scl)
  df_out <- as.data.frame(df_pca$x)
  df_out$group<-group
  df_out$sample_name<-samplenames
  
  percentage <- round(df_pca$sdev^2 / sum(df_pca$sdev^2)*100,2)
  percentage <- paste0(colnames(df_out)[grep("^PC",colnames(df_out))], " (", paste0(as.character(percentage), "% variance", ")") )
  print(percentage)
  
  if(is.null(shape)){
    p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group))
    if(!is.null(colors)){p<-p+scale_color_manual(values=colors,name=LegendName)}
  }else{
    p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group,shape=shape))
    if(!is.null(colors)){p<-p+scale_color_manual(values=colors,name=LegendName)}
  }
  
  
  p<-p+geom_point(size=point.size) + 
    xlab(percentage[1]) + 
    ylab(percentage[2]) +
    ggtitle(label = "PCA plot") +
    labs(shape=LegendName_Shape, col=LegendName_Color)+
    ggTheme(1)
  if(ggrepelLab){
    p<-p + geom_text_repel(aes(label=sample_name),show.legend = FALSE,
                           size=size_gglab,
                           force = 2,max.overlaps = Inf)
  }
  
  if(MahalanobisEllips){
    p<-p+stat_ellipse()
  }
  return(p)
  

}

      
