plot_2DPCA<-function(expression,group,colors=NULL,shape=NULL,samplenames,title="",LegendName_Color="group",LegendName_Shape="shape",LegendName="group",ggrepelLab=TRUE,size_gglab=5,size_title=14){
  if(is.matrix(expression)){
    df_pca<-prcomp(t(expression),scale=TRUE)
    df_out <- as.data.frame(df_pca$x)
    df_out$group<-group
    df_out$sample_name<-samplenames
    
    percentage <- round(df_pca$sdev^2 / sum(df_pca$sdev^2)*100,2)
    percentage <- paste0(colnames(df_out)[grep("^PC",colnames(df_out))], "(", paste0(as.character(percentage), "%", ")") )
    print(percentage)
    
    if(is.null(shape)){
      p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group))
      if(!is.null(colors)){p<-p+scale_color_manual(values=colors,name=LegendName)}
    }else{
      p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group,shape=shape))
      if(!is.null(colors)){p<-p+scale_color_manual(values=colors,name=LegendName)}
    }
    
    
    p<-p+geom_point(size=4)+ xlab(percentage[1]) + ylab(percentage[2])+
      #geom_text(label=row.names(df_out_raw),show.legend = FALSE,hjust=0,vjust=0.2)
      theme(plot.title = element_text(size = size_title, face = "bold",hjust = 0.5),
            legend.title=element_text(size=14), legend.text=element_text(size=12,),
            axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
            panel.background = element_blank(),
            panel.grid.major =  element_line(colour = "grey90", size = 0.2),
            panel.grid.minor =  element_line(colour = "grey98", size = 0.5),
            panel.border = element_rect(color='black',fill=NA))
    if(ggrepelLab){
      p<-p + geom_text_repel(aes(label=sample_name),show.legend = FALSE,size=size_gglab,force = 2) +
        labs(shape=LegendName_Shape, col=LegendName_Color)+
        ggtitle(title)
    }
  } else if(typeof(expression)=="S4"){
    pcaData <- plotPCA(vst(expression, fitType = "local"),
                       intgroup = intgroup,
                       returnData = TRUE
    )
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, 
           aes(PC1, PC2, color = group)) +
      geom_point(size = point.size) +
      geom_label(aes(label = colnames(dds))) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      guides(color = guide_legend(ncol = col.legend)) +
      ggtitle(title)
  
    
  }
  
  return(p)
}